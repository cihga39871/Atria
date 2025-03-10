
# f_procs(x::String) = x == "-p" || x == "--procs"

function julia_wrapper_atria_se(ARGS::Vector{String}; exit_after_help = true)

    time_program_initializing = time()

    args = parsing_args(ARGS; exit_after_help = exit_after_help)

    if args === nothing  # ARGS is ["-h"]
        return 0
    end
    args_range_test(args)

    outdir = args["output-dir"]

    nfile = length(args["read1"])
    file_range = 1:nfile

    #================== Arguments ====================#

    command                  = `$ARGS`
    nthread                  =  args["threads"               ]
    max_chunk_size           =  2 ^ args["log2-chunk-size"   ]
    force                    = args["force"]
    # polyX
    min_poly_length          = args["poly-length"            ]
    poly_mismatch_per_16mer  = args["poly-mismatch-per-16mer"]
    # N
    max_N                    =  args["max-n"                 ]
    # adapter
    adapter1s                = args["adapter1"]
    adapter1_seqheadsets     = SeqHeadSet.(adapter1s)
    # NOTE: TruncSeq has some unknown accuracy problems.
    kmer_tolerance           = args["kmer-tolerance"          ]
    kmer_tolerance_consensus = args["kmer-tolerance-consensus"]
    trim_score               = args["trim-score-se"           ]
    tail_length              = args["tail-length"             ]
    # consensus
    # hard clip
    nclip_end          = args["clip3-r1"     ]
    nclip_after        = args["clip-after-r1"]
    nclip_front        = args["clip5-r1"     ]
    # quality
    quality_offset     = Trimmer.get_quality_offset(args["quality-format"])
    quality_kmer       = args["quality-kmer"]
    quality_single     = convert(UInt8, args["quality-score"] + quality_offset)
    quality_threshold  = convert(UInt64, quality_single * quality_kmer)
    # length
    length_range = Trimmer.get_length_range(args["length-range"]::String)
    # complexity
    min_complexity = args["min-complexity"]

    # feature enable/disable
    do_check_identifier      =  false
    do_pcr_dedup             =  args["pcr-dedup"           ]
    do_pcr_dedup_count       =  args["pcr-dedup-count"     ]
    do_hash_collision        =  get(ENV, "ATRIA_HASH_COLLISION", "") in ("1", "true", "True", "T") ? true : false
    do_polyG                 =  args["polyG"               ]
    do_polyT                 =  args["polyT"               ]
    do_polyA                 =  args["polyA"               ]
    do_polyC                 =  args["polyC"               ]
    do_length_filtration     = !args["no-length-filtration"]
    do_adapter_trimming      = !args["no-adapter-trim"     ]
    do_consensus_calling     =  false
    do_hard_clip_end         =  nclip_end > 0
    do_hard_clip_after       =  nclip_after > 0
    do_hard_clip_front       =  nclip_front > 0
    do_quality_trimming      = !args["no-quality-trim"     ]
    do_tail_n_trimming       = !args["no-tail-n-trim"      ]
    do_tail_low_qual_trimming=  do_polyG || do_polyT || do_polyA || do_polyC || do_adapter_trimming
    do_max_n_filtration      =  max_N > 0
    do_read_stats            =  args["stats"]
    do_complexity_filtration =  args["enable-complexity-filtration"]

    mkpath(outdir)


    #================== Main function and common variables ====================#
    dup_dict = Dict{Vector{UInt64}, DupCount}()
    dup_dict_lock = ReentrantLock()
    
    stats_r1 = TrimStats()

    if do_hash_collision
        error("Dev: Hash Collision not supported for single read.")
    end

    #======= Check identifier (not applicable in single end)=======#
    CheckIdentifier = nothing

    #======= PolyX =======#
    PolyG = do_polyG ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_G, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            @atomic $stats_r1.polyG += 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyT = do_polyT ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_T, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            @atomic $stats_r1.polyT += 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyA = do_polyA ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_A, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            @atomic $stats_r1.polyA += 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyC = do_polyC ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_C, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            @atomic $stats_r1.polyC += 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    #======= Complexity =======#
    ComplexityFilter = do_complexity_filtration ? quote
        if seq_complexity(r1) < $min_complexity
            @atomic $stats_r1.complexity_filtered += 1
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing

    #======= Hard clips =======#
    HardClipEndR1 = do_hard_clip_end ? quote
        len_clip = max(0, length(r1.seq) - $nclip_end_r1)
        tail_trim!(r1::FqRecord, len_clip)
    end : nothing
    HardClipEndR2 = nothing

    HardClipAfterR1 = do_hard_clip_after ? quote
        if length(r1.seq) > $nclip_after    
            @atomic $stats_r1.hard_clip_after += 1
            tail_trim!(r1::FqRecord, $nclip_after)
        end
    end : nothing
    HardClipAfterR2 = nothing

    HardClipFrontR1 = do_hard_clip_front ? quote
        front_trim!(r1::FqRecord, $nclip_front::Int64)
    end : nothing
    HardClipFrontR2 = nothing

    #======= Trim N from the tail =======#
    TailNTrim = do_tail_low_qual_trimming ? quote
        tail_low_qual_trim!(r1::FqRecord, $stats_r1)
    end : do_tail_n_trimming ? quote
        tail_N_trim!(r1::FqRecord, $stats_r1)
    end : nothing


    #======= If reads are too short, remove =======#
    LengthFilter = do_length_filtration ? quote
        if !isinreadlength!(r1::FqRecord, $length_range)
            @atomic $stats_r1.length_filtered += 1            
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing


    #======= If reads have so many N, remove =======#
    MaxNFilter = do_max_n_filtration ? quote
        if !isnotmuchN!(r1::FqRecord, $max_N)
            @atomic $stats_r1.max_n_filtered += 1            
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing

    #======= PCR dedup =======#
    PCRDedup = if do_pcr_dedup && do_pcr_dedup_count
        quote
            hash_key = hash_dna(r1.seq)
            lock($dup_dict_lock)
            dup_count = get!($dup_dict, hash_key) do
                DupCount(0, String(copy(r1.id)))
            end
            unlock($dup_dict_lock)
            
            @atomic dup_count.count += 1
            if dup_count.count > 1  # dup, not unique
                is_good = false; @goto stop_read_processing # return false
            end
        end
    elseif do_pcr_dedup
        quote
            hash_key = hash_dna(r1.seq)
            lock($dup_dict_lock)
            dup_count = get!($dup_dict, hash_key) do
                DupCount(0)
            end
            unlock($dup_dict_lock)
            
            @atomic dup_count.count += 1
            if dup_count.count > 1  # dup, not unique
                is_good = false; @goto stop_read_processing # return false
            end
        end
    else
        nothing
    end

    #======= Quality trimming =======#
    QualityTrim = do_quality_trimming ? quote
        quality_index1 = qualitymatch(r1::FqRecord, $quality_single::UInt8, $quality_threshold::UInt64, $quality_kmer::Int64)
        if quality_index1 != -1  # -1 means no need to do quality trimming
            @atomic $stats_r1.quality_trim += 1
            tail_trim!(r1::FqRecord, quality_index1::Int64)
        end
    end : nothing


    #======= adapter trimming =======#
    AdapterTrim = do_adapter_trimming ? quote
        nremain = adapter_match_se($adapter1_seqheadsets, r1, $kmer_tolerance, $trim_score)
        if nremain < length(r1.seq)
            @atomic $stats_r1.adapter_trim += 1
            tail_trim!(r1::FqRecord, nremain)
        end
    end : nothing

    ReadProcess = quote end
    steps = OrderedDict{String, Vector{Union{Nothing,Expr}}}(
        "DefaultOrder" => [CheckIdentifier, PolyG, PolyT, PolyA, PolyC, LengthFilter, AdapterTrim, HardClipEndR1, HardClipEndR2, HardClipAfterR1, HardClipAfterR2, HardClipFrontR1, HardClipFrontR2, QualityTrim, TailNTrim, MaxNFilter, LengthFilter, ComplexityFilter, PCRDedup],
        "CheckIdentifier" => [CheckIdentifier],
        "PolyX" => [PolyG, PolyT, PolyA, PolyC],
        "PolyG" => [PolyG],
        "PolyT" => [PolyT],
        "PolyA" => [PolyA],
        "PolyC" => [PolyC],
        "AdapterTrim" => [AdapterTrim],
        "HardClipEnd"     => [HardClipEndR1, HardClipEndR2],
        "HardClipEndR1"   => [HardClipEndR1],
        "HardClipEndR2"   => [HardClipEndR2],
        "HardClipAfter" => [HardClipAfterR1, HardClipAfterR2],
        "HardClipAfterR1" => [HardClipAfterR1],
        "HardClipAfterR2" => [HardClipAfterR2],
        "HardClipFront" => [HardClipFrontR1, HardClipFrontR2],
        "HardClipFrontR1" => [HardClipFrontR1],
        "HardClipFrontR2" => [HardClipFrontR2],
        "QualityTrim" => [QualityTrim],
        "TailNTrim" => [TailNTrim],
        "MaxNFilter" => [MaxNFilter],
        "LengthFilter" => [LengthFilter],
        "ComplexityFilter" => [ComplexityFilter],
        "PCRDedup" => [PCRDedup]
    )

    for step in args["order"]
        step_exprs = get(steps, step, step)
        if step_exprs isa AbstractString
            name_list = join(collect(keys(steps)), ", ")
            error("Invalid process name in -O or --order: $step. Possible names are $name_list.")
        end
        for step_expr in step_exprs
            if isnothing(step_expr)
                continue
            end
            append!(ReadProcess.args, step_expr.args)
        end
    end

    # moved from thread_trim.jl to fix world age error.
    function processing_reads_range!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, reads_range::UnitRange{Int64})
        this_threadid = Threads.threadid()
        for i in reads_range
            @inbounds isgoods[i] = read_processing!(r1s[i], this_threadid)
        end
    end

    # processing_reads_threads! without dedup
    @eval function processing_reads_threads!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
        if length(isgoods) != n_reads
            resize!(isgoods, n_reads)
        end
        # split reads to N reads per batch
        Threads.@threads :static for reads_start in 1:512:n_reads
            reads_end = min(reads_start + 511, n_reads)
            reads_range = reads_start:reads_end
            thread_id = Threads.threadid()
    
            for i in reads_range
                r1 = r1s[i]
                is_good = true
                $ReadProcess
                @label stop_read_processing
                @inbounds isgoods[i] = is_good
            end
        end
    end
    
    in1bytes = Vector{UInt8}(undef, max_chunk_size)

    # number of jobs to boxing FqRecord from UInt8 Vector
    njobs = nthread * 5
    vr1s = ntuple(_ -> Vector{FqRecord}(), njobs)

    r1s = Vector{FqRecord}()

    isgoods_even = Vector{Bool}()
    isgoods_odd = Vector{Bool}()

    outr1s = Vector{Vector{UInt8}}()

    time_program_initializing = time() - time_program_initializing

    #================== Iteration for files ====================#
    for filenum in file_range
        # filenum = 1
        time_file_initializing = time()


        #===== file names =====#

        file1 = args["read1"][filenum]

        # check whether this sample is processed before
        outjson = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.log.json", count=1))
        if force
            rm(outjson, force=true)
        else
            if isfile(outjson) && filesize(outjson) > 0
                @warn "Skip completed analysis: $outjson (use --force to disable the feature)" _module=nothing _group=nothing _id=nothing _file=nothing
                continue
            end
        end

        isingzip = occursin(r"\.gz$"i, file1)
        isinbzip2 = occursin(r"\.bz2$"i, file1)
        outcompress = uppercase(args["compress"])
        if outcompress == "AUTO"
            if isingzip
                outcompress = "GZIP"
            elseif isinbzip2
                outcompress = "BZIP2"
            else
                outcompress = "NO"
            end
        elseif outcompress == "GZ"
            outcompress = "GZIP"
        elseif outcompress == "BZ2"
            outcompress = "BZIP2"
        end

        outfile1 = joinpath(outdir, replace(basename(file1), r"(fastq$|fq$|[^.]*)(\.gz|\.bz2)?$"i => s"atria.\1", count=1))
        if outcompress == "GZIP"
            outfile1 *= ".gz"
        elseif outcompress == "BZIP2"
            outfile1 *= ".bz2"
        end

        outlog = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.log", count=1))


        #===== file IO =====#

        if isingzip
            io1 = open(`pigz -p$nthread -cd $file1`, write=false)
        elseif isinbzip2
            io1 = open(`pbzip2 -p$nthread -cd $file1`, write=false)
        else
            io1 = open(file1, "r")
        end

        if outcompress == "GZIP"
            io1out = open(`pigz -p$nthread -c`, open(outfile1, "w"), write=true, read=false)
        elseif outcompress == "BZIP2"
            io1out = open(`pbzip2 -p$nthread -c`, open(outfile1, "w"), write=true, read=false)
        else
            io1out = open(outfile1, "w")
        end

        iolog = open(outlog, "w+")
        logger = Logging.SimpleLogger(iolog)

        @info "ATRIA VERSIONS" atria=atria_version julia=string("v", VERSION)
        @info "ATRIA ARGUMENTS" command
        @info "ATRIA OUTPUT FILES" read=outfile1
        @info("ATRIA TRIMMERS AND FILTERS",
            tail_polyG_trimming = do_polyG,
            tail_polyT_trimming = do_polyT,
            tail_polyA_trimming = do_polyA,
            tail_polyC_trimming = do_polyC,
            adapter_trimming = do_adapter_trimming,
            consensus_calling = do_consensus_calling,
            hard_clip_3_end = do_hard_clip_after,
            hard_clip_5_end = do_hard_clip_front,
            quality_trimming = do_quality_trimming,
            tail_N_trimming = do_tail_n_trimming,
            max_N_filtering = do_max_n_filtration,
            length_filtering = do_length_filtration,
            complexity_filtering = do_complexity_filtration)

        with_logger(logger) do
            @info "ATRIA VERSIONS" atria=atria_version julia=string("v", VERSION)
            @info "ATRIA ARGUMENTS" command
            @info "ATRIA OUTPUT FILES" read=outfile1
            @info("ATRIA TRIMMERS AND FILTERS",
                tail_polyG_trimming = do_polyG,
                tail_polyT_trimming = do_polyT,
                tail_polyA_trimming = do_polyA,
                tail_polyC_trimming = do_polyC,
                adapter_trimming = do_adapter_trimming,
                consensus_calling = do_consensus_calling,
                hard_clip_3_end = do_hard_clip_after,
                hard_clip_5_end = do_hard_clip_front,
                quality_trimming = do_quality_trimming,
                tail_N_trimming = do_tail_n_trimming,
                max_N_filtering = do_max_n_filtration,
                length_filtering = do_length_filtration,
                complexity_filtering = do_complexity_filtration)
        end


        #================== Renew variables for read processing ====================#


        # setting chunk size for file 1 and file2: skip for single end

        # clear common variables
        empty!(r1s)
        empty!(dup_dict)
        empty!(stats_r1)

        #=
        global n_reads
        global n_r1
        global nbatch
        global total_read_copied_in_loading
        global total_n_goods
        global total_n_reads
        global total_n_bytes_read1
        global chunk_size1
        =#

        n_reads = 0
        n_r1 = 0
        nbatch = 0
        total_read_copied_in_loading = 0
        total_n_goods = 0
        total_n_reads = 0
        total_n_pcr_dup = 0
        total_n_bytes_read1 = 0
        in1bytes_nremain = 0

        #================== File processing ====================#
        time_file_initializing = time() - time_file_initializing
        time_read_processing = time()

        # the first cycle to generate compiled code?
        function cycle_wrapper_single_end(task_r1s_unbox, task_write1)
            nbatch += 1
            n_r1_before = length(r1s) - n_reads

            if typeof(io1) <: IOStream  # not compressed
                (n_r1, r1s, ncopied) = load_fqs_threads!(io1, in1bytes, vr1s, r1s, task_r1s_unbox; remove_first_n = n_reads, njobs=njobs)
            else  # gziped
                total_n_bytes_read1 += length(in1bytes)  # will read INT in this batch
                (n_r1, r1s, in1bytes_nremain, ncopied) = load_fqs_threads!(
                    io1,
                    in1bytes,
                    in1bytes_nremain,
                    vr1s,
                    r1s,
                    task_r1s_unbox;
                    remove_first_n = n_reads,
                    njobs = njobs
                )
            end

            n_reads = n_r1
            total_read_copied_in_loading += ncopied

            isgoods = nbatch % 2 == 0 ? isgoods_even : isgoods_odd  # write task is still using the other one! 
            
            if do_pcr_dedup && nbatch == 1
                # size hint estimate
                chunk_size1, chunk_size2, uncompressed_size1, uncompressed_size2 = chunk_sizes(file1, file1, max_chunk_size)

                if uncompressed_size1 == -1
                    if isingzip
                        uncompressed_size = filesize(file1) / 0.14
                    elseif isinbzip2
                        uncompressed_size = filesize(file1) / 0.12
                    else
                        uncompressed_size = filesize(file1)
                    end
                else
                    uncompressed_size = uncompressed_size1
                end

                processed_rate = length(in1bytes) / uncompressed_size
                estimated_dict_size = round(Int, n_reads / processed_rate * 1.2)
                sizehint!(dup_dict, estimated_dict_size)
            end

            Base.invokelatest(processing_reads_threads!, r1s, isgoods, n_reads)

            isgoods_in_range = view(isgoods, 1:n_reads)
            task_sum = Threads.@spawn sum(isgoods_in_range)

            task_r1s_unbox, task_write1 = write_fqs_threads!(
                io1out::IO,
                outr1s::Vector{Vector{UInt8}},
                r1s::Vector{FqRecord},
                n_reads::Int, isgoods_in_range, task_write1)

            n_goods = fetch(task_sum)
            total_n_goods += n_goods
            total_n_reads += n_reads


            @info "Cycle $nbatch: read $n_reads/$total_n_reads; wrote $n_goods/$total_n_goods; (copied $ncopied/$total_read_copied_in_loading)"

            return task_r1s_unbox, task_write1
        end

        task_r1s_unbox = Threads.@spawn 1
        task_write1 = Threads.@spawn 1

        while !eof(io1::IO)
            task_r1s_unbox, task_write1 = cycle_wrapper_single_end(task_r1s_unbox, task_write1)
        end

        @info "ATRIA COMPLETE" read=outfile1
        with_logger(logger) do
            @info "ATRIA COMPLETE" read=outfile1
        end

        wait(task_write1)

        time_read_processing = time() - time_read_processing


        #================== Close files ====================#
        time_post_processing = time()

        close(io1)
        close(io1out)
        close(iolog)

        #================== Dedup stats ====================#
        if do_pcr_dedup 
            if do_pcr_dedup_count
                out_pcr_dedup_count = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.pcr_dedup_count.tsv", count=1))
                @atomic stats_r1.pcr_dedup_removed = write_pcr_dedup_count(out_pcr_dedup_count, dup_dict)
            else
                @atomic stats_r1.pcr_dedup_removed = get_dup_count(dup_dict)
            end
        end

        #================== JSON log ====================#

        logjson = OrderedDict()

        logjson["arguments"] = sort!(OrderedDict(args))
        logjson["version"] = OrderedDict(
            "julia-version" => VERSION,
            "atria-version" => atria_version
        )
        logjson["input-files"] = OrderedDict(
            "input-read1" => file1,
        )
        logjson["output-files"] = OrderedDict(
            "good-read1" => outfile1,
            "log" => outlog,
            "json-log" => outjson
        )
        time_post_processing = time() - time_post_processing
        logjson["runtime"] = OrderedDict(
            "program-initializing-time" => time_program_initializing,
            "file-initializing-time"=> time_file_initializing,
            "read-processing-time" => time_read_processing,
            "post-processing-time"=> time_post_processing
        )
        logjson["trimming-details"] = OrderedDict(
            "good-reads" => total_n_goods,
            "total-reads" => total_n_reads,
        )
        logjson["stats"] = stats_r1

        iologjson = open(outjson, "w+")
        JSON.print(iologjson, sort!(logjson), 4)
        close(iologjson)

    end

    return 0
end # func

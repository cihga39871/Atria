
# f_procs(x::String) = x == "-p" || x == "--procs"

function julia_wrapper_atria_pe(ARGS::Vector{String}; exit_after_help = true)

    time_program_initializing = time()

    args = parsing_args(ARGS; exit_after_help = exit_after_help)

    if args === nothing  # ARGS is ["-h"]
        return 0
    end
    args_range_test(args)
    
    nthread = args["threads"]
    outdir = args["output-dir"]

    nfile = length(args["read1"])
    file_range = 1:nfile

    #================== Arguments ====================#

    command                  = `$ARGS`
    max_chunk_size           =  2 ^ args["log2-chunk-size"]
    force                    = args["force"]
    # polyX
    min_poly_length          = args["poly-length"            ]
    poly_mismatch_per_16mer  = args["poly-mismatch-per-16mer"]
    # N
    max_N                    = args["max-n"                  ]
    # adapter
    adapter1s                 = LongDNA{4}.(args["adapter1"])
    adapter2s                 = LongDNA{4}.(args["adapter2"])
    adapter1_seqheadsets      = SeqHeadSet.(adapter1s)
    adapter2_seqheadsets      = SeqHeadSet.(adapter2s)
    # first adapter
    adapter1                 = adapter1s[1]
    adapter2                 = adapter2s[1]
    adapter1_seqheadset      = adapter1_seqheadsets[1]
    adapter2_seqheadset      = adapter2_seqheadsets[1]
    # adapter and primer match options
    op                        = PEOptions(args)

    # hard clip
    nclip_after_r1        = args["clip-after-r1"]
    nclip_after_r2        = args["clip-after-r2"]
    nclip_front_r1        = args["clip5-r1"     ]
    nclip_front_r2        = args["clip5-r2"     ]
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
    do_check_identifier      =  args["check-identifier"    ]
    do_pcr_dedup             =  args["pcr-dedup"           ]
    do_pcr_dedup_count       =  args["pcr-dedup-count"     ]
    do_polyG                 =  args["polyG"               ]
    do_polyT                 =  args["polyT"               ]
    do_polyA                 =  args["polyA"               ]
    do_polyC                 =  args["polyC"               ]
    do_length_filtration     = !args["no-length-filtration"]
    do_adapter_trimming      = !args["no-adapter-trim"     ]
    do_consensus_calling     = !args["no-consensus"]
    do_hard_clip_3_end_r1    =  nclip_after_r1 > 0
    do_hard_clip_3_end_r2    =  nclip_after_r2 > 0
    do_hard_clip_5_end_r1    =  nclip_front_r1 > 0
    do_hard_clip_5_end_r2    =  nclip_front_r2 > 0
    do_quality_trimming      = !args["no-quality-trim"     ]
    do_tail_n_trimming       = !args["no-tail-n-trim"      ]
    do_tail_low_qual_trimming=  do_polyG || do_polyT || do_polyA || do_polyC || do_adapter_trimming
    do_max_n_filtration      =  max_N > 0
    # do_read_stats            =  args["stats"                       ]
    do_complexity_filtration =  args["enable-complexity-filtration"]

    mkpath(outdir)


    #================== Main function and common variables ====================#
    r1_seq_rc_threads = [LongDNA{4}() for _ in 1:nthread]
    r2_seq_rc_threads = [LongDNA{4}() for _ in 1:nthread]
    
    dup_dict = Dict{Vector{UInt64}, DupCount}()
    dup_dict_lock = ReentrantLock()

    #======= Check identifier =======#
    CheckIdentifier = do_check_identifier ? quote
        if !check_identifier(r1, r2)
            throw_identifier_error(r1, r2)
        end
    end : nothing

    #======= PolyX =======#
    PolyG = do_polyG ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_G, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
        polyX_idx, polyX_length = polyX_tail_scan(DNA_G, r2, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r2, polyX_idx)
        end
    end : nothing

    PolyT = do_polyT ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_T, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
        polyX_idx, polyX_length = polyX_tail_scan(DNA_T, r2, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r2, polyX_idx)
        end
    end : nothing

    PolyA = do_polyA ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_A, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
        polyX_idx, polyX_length = polyX_tail_scan(DNA_A, r2, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r2, polyX_idx)
        end
    end : nothing

    PolyC = do_polyC ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_C, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
        polyX_idx, polyX_length = polyX_tail_scan(DNA_C, r2, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r2, polyX_idx)
        end
    end : nothing

    #======= Complexity =======#
    ComplexityFilter = do_complexity_filtration ? quote
        if seq_complexity(r1) < $min_complexity
            is_good = false; @goto stop_read_processing # return false
        elseif seq_complexity(r2) < $min_complexity
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing

    #======= Hard clips =======#
    HardClip3EndR1 = do_hard_clip_3_end_r1 ? quote
        tail_trim!(r1::FqRecord, $nclip_after_r1)
    end : nothing
    HardClip3EndR2 = do_hard_clip_3_end_r2 ? quote
        tail_trim!(r2::FqRecord, $nclip_after_r2)
    end : nothing

    HardClip5EndR1 = do_hard_clip_5_end_r1 ? quote
        front_trim!(r1::FqRecord, $nclip_front_r1::Int64)
    end : nothing
    HardClip5EndR2 = do_hard_clip_5_end_r2 ? quote
        front_trim!(r2::FqRecord, $nclip_front_r2::Int64)
    end : nothing


    #======= Trim N from the tail =======#
    TailNTrim = do_tail_low_qual_trimming ? quote
        tail_low_qual_trim!(r1::FqRecord)
        tail_low_qual_trim!(r2::FqRecord)
    end : do_tail_n_trimming ? quote
        tail_N_trim!(r1::FqRecord)
        tail_N_trim!(r2::FqRecord)
    end : nothing


    #======= If reads are too short, remove =======#
    LengthFilter = do_length_filtration ? quote
        if !isinreadlength!(r1::FqRecord, r2::FqRecord, $length_range)
            # global count_reads_pair_invalid_length[thread_id] += 1
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing


    #======= If reads have so many N, remove =======#
    MaxNFilter = do_max_n_filtration ? quote
        if !isnotmuchN!(r1::FqRecord, r2::FqRecord, $max_N)
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing

    #======= PCR dedup =======#
    PCRDedup = do_pcr_dedup ? quote
        hash_key = hash_dna(r1.seq, r2.seq)
        lock($dup_dict_lock)
        dup_count = get($dup_dict, hash_key, nothing)
        
        if isnothing(dup_count)  # unique
            new_val = DupCount(1, String(r1.seq) * "__" * String(r2.seq))
            # new_val = DupCount(1, String(copy(r1.id)))
            $dup_dict[hash_key] = new_val
            unlock($dup_dict_lock)
        else # dup
            unlock($dup_dict_lock)
            @atomic dup_count.count += 1
            this = String(r1.seq) * "__" * String(r2.seq) 
            if this != dup_count.id
                s1, s2 = split(dup_count.id, "__")
                if s1[1:end-1] != String(r1.seq)[1:end-1] || s2[1:end-1] != String(r2.seq)[1:end-1]
                    @warn "conflict" this that = dup_count.id
                end
            end
            is_good = false; @goto stop_read_processing # return false
        end
    end : nothing

    #======= Quality trimming =======#
    QualityTrim = do_quality_trimming ? quote
        quality_index1 = qualitymatch(r1::FqRecord, $quality_single::UInt8, $quality_threshold::UInt64, $quality_kmer::Int64)
        if quality_index1 != -1  # -1 means no need to do quality trimming
            # global count_r1_quality_trimming[thread_id] += 1
            tail_trim!(r1::FqRecord, quality_index1::Int64)
        end

        quality_index2 = qualitymatch(r2::FqRecord, $quality_single::UInt8, $quality_threshold::UInt64, $quality_kmer::Int64)
        if quality_index2 != -1  # -1 means no need to do quality trimming
            # global count_r2_quality_trimming[thread_id] += 1
            tail_trim!(r2::FqRecord, quality_index2::Int64)
        end
    end : nothing


    #======= adapter trimming =======#
    AdapterTrim = if do_adapter_trimming
        if length(adapter1_seqheadsets) > 1
            quote
                r1_seq_rc = $r1_seq_rc_threads[thread_id]
                r2_seq_rc = $r2_seq_rc_threads[thread_id]
                adapter_match_and_trim_pe!(
                    $adapter1_seqheadsets, $adapter2_seqheadsets,
                    $adapter1s, $adapter2s, r1, r2,
                    true, r1_seq_rc, r2_seq_rc, $op
                )
            end
        else
            quote
                r1_seq_rc = $r1_seq_rc_threads[thread_id]
                r2_seq_rc = $r2_seq_rc_threads[thread_id]
                adapter_match_and_trim_pe!(
                    $adapter1_seqheadset, $adapter2_seqheadset,
                    $adapter1, $adapter2, r1, r2,
                    r1_seq_rc, r2_seq_rc, $op
                )
            end
        end
    else
        nothing
    end


    ReadProcess = quote end
    steps = OrderedDict{String, Vector{Union{Nothing,Expr}}}(
        "DefaultOrder" => [CheckIdentifier, PolyG, PolyT, PolyA, PolyC, LengthFilter, AdapterTrim, HardClip3EndR1, HardClip3EndR2, HardClip5EndR1, HardClip5EndR2, QualityTrim, TailNTrim, MaxNFilter, LengthFilter, ComplexityFilter, PCRDedup],
        "CheckIdentifier" => [CheckIdentifier],
        "PolyX" => [PolyG, PolyT, PolyA, PolyC],
        "PolyG" => [PolyG],
        "PolyT" => [PolyT],
        "PolyA" => [PolyA],
        "PolyC" => [PolyC],
        "AdapterTrim" => [AdapterTrim],
        "HardClip3End" => [HardClip3EndR1, HardClip3EndR2],
        "HardClip3EndR1" => [HardClip3EndR1],
        "HardClip3EndR2" => [HardClip3EndR2],
        "HardClip5End" => [HardClip5EndR1, HardClip5EndR2],
        "HardClip5EndR1" => [HardClip5EndR1],
        "HardClip5EndR2" => [HardClip5EndR2],
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

    # processing_reads_threads! without dedup
    @eval function processing_reads_threads!(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
        if length(isgoods) != n_reads
            resize!(isgoods, n_reads)
        end

        # split reads to N reads per batch
        Threads.@threads :static for reads_start in 1:1024:n_reads
            reads_end = min(reads_start + 1023, n_reads)
            reads_range = reads_start:reads_end
            thread_id = Threads.threadid()
            @inbounds for i in reads_range
                r1 = r1s[i]
                r2 = r2s[i]
                is_good = true
                $ReadProcess
                @label stop_read_processing
                @inbounds isgoods[i] = is_good
            end
        end
    end

    in1bytes = Vector{UInt8}(undef, max_chunk_size)
    in2bytes = Vector{UInt8}(undef, max_chunk_size)

    # number of jobs to boxing FqRecord from UInt8 Vector
    njobs = nthread * 6 - 1
    vr1s = ntuple(_ -> Vector{FqRecord}(), njobs)
    vr2s = ntuple(_ -> Vector{FqRecord}(), njobs)

    r1s = Vector{FqRecord}()
    r2s = Vector{FqRecord}()

    isgoods_odd = Vector{Bool}()
    isgoods_even = Vector{Bool}()

    outr1s = Vector{Vector{UInt8}}()
    outr2s = Vector{Vector{UInt8}}()

    time_program_initializing = time() - time_program_initializing

    #================== Iteration for paired files ====================#
    for filenum in file_range
        # filenum = 1
        time_file_initializing = time()


        #===== file names =====#

        file1 = args["read1"][filenum]
        file2 = args["read2"][filenum]

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
        outfile2 = joinpath(outdir, replace(basename(file2), r"(fastq$|fq$|[^.]*)(\.gz|\.bz2)?$"i => s"atria.\1", count=1))
        if outcompress == "GZIP"
            outfile1 *= ".gz"
            outfile2 *= ".gz"
        elseif outcompress == "BZIP2"
            outfile1 *= ".bz2"
            outfile2 *= ".bz2"
        end

        outlog = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.log", count=1))


        #===== file IO =====#
        halfthread = cld(nthread, 2)
        if isingzip
            io1 = open(`pigz -p$halfthread -cd $file1`, write=false)
            io2 = open(`pigz -p$halfthread -cd $file2`, write=false)
        elseif isinbzip2
            io1 = open(`pbzip2 -p$halfthread -cd $file1`, write=false)
            io2 = open(`pbzip2 -p$halfthread -cd $file2`, write=false)
        else
            io1 = open(file1, "r")
            io2 = open(file2, "r")
        end

        if outcompress == "GZIP"
            io1out = open(`pigz -p$halfthread -c`, open(outfile1, "w"), write=true, read=false)
            io2out = open(`pigz -p$halfthread -c`, open(outfile2, "w"), write=true, read=false)
        elseif outcompress == "BZIP2"
            io1out = open(`pbzip2 -p$halfthread -c`, open(outfile1, "w"), write=true, read=false)
            io2out = open(`pbzip2 -p$halfthread -c`, open(outfile2, "w"), write=true, read=false)
        else
            io1out = open(outfile1, "w")
            io2out = open(outfile2, "w")
        end

        iolog = open(outlog, "w+")
        logger = Logging.SimpleLogger(iolog)

        @info "ATRIA VERSIONS" atria=atria_version julia=string("v", VERSION)
        @info "ATRIA ARGUMENTS" command
        @info "ATRIA OUTPUT FILES" read1=outfile1 read2=outfile2
        @info("ATRIA TRIMMERS AND FILTERS",
            pcr_dedup = do_pcr_dedup,
            adapter_trimming = do_adapter_trimming,
            consensus_calling = do_consensus_calling,
            hard_clip_3_end_r1 = do_hard_clip_3_end_r1,
            hard_clip_3_end_r2 = do_hard_clip_3_end_r2,
            hard_clip_5_end_r1 = do_hard_clip_5_end_r1,
            hard_clip_5_end_r2 = do_hard_clip_5_end_r2,
            quality_trimming = do_quality_trimming,
            tail_N_trimming = do_tail_n_trimming,
            max_N_filtering = do_max_n_filtration,
            length_filtering = do_length_filtration)

        with_logger(logger) do
            @info "ATRIA VERSIONS" atria=atria_version julia=string("v", VERSION)
            @info "ATRIA ARGUMENTS" command
            @info "ATRIA OUTPUT FILES" read1=outfile1 read2=outfile2
            @info("ATRIA TRIMMERS AND FILTERS",
                pcr_dedup = do_pcr_dedup,
                adapter_trimming = do_adapter_trimming,
                consensus_calling = do_consensus_calling,
                hard_clip_3_end_r1 = do_hard_clip_3_end_r1,
                hard_clip_3_end_r2 = do_hard_clip_3_end_r2,
                hard_clip_5_end_r1 = do_hard_clip_5_end_r1,
                hard_clip_5_end_r2 = do_hard_clip_5_end_r2,
                quality_trimming = do_quality_trimming,
                tail_N_trimming = do_tail_n_trimming,
                max_N_filtering = do_max_n_filtration,
                length_filtering = do_length_filtration)
        end


        #================== Renew variables for read processing ====================#


        # setting chunk size for file 1 and file2
        chunk_size1, chunk_size2, uncompressed_size1, uncompressed_size2 = chunk_sizes(file1, file2, max_chunk_size)
        if (uncompressed_size1 == -1 || uncompressed_size2 == -1) && (isingzip || isinbzip2)
            # file is gzip but uncompressed size not known.
            # do not resize. just assume R1/2 is the original data, which means insert size is evenly-distributed.
            chunk_size1 = length(in1bytes)
            chunk_size2 = length(in2bytes)
        else
            resize!(in1bytes, chunk_size1)
            resize!(in2bytes, chunk_size2)
        end

        # clear common variables
        empty!(r1s)
        empty!(r2s)
        empty!(dup_dict)

        n_reads = 0
        n_r1 = 0
        n_r2 = 0
        nbatch = 0
        total_read_copied_in_loading = 0
        total_n_goods = 0
        total_n_reads = 0
        total_n_pcr_dup = 0
        total_n_bytes_read1 = 0
        total_n_bytes_read2 = 0
        in1bytes_nremain = 0
        in2bytes_nremain = 0

        #================== File processing ====================#
        time_file_initializing = time() - time_file_initializing
        time_read_processing = time()

        # the first cycle to generate compiled code?
        function cycle_wrapper(task_r1s_unbox, task_r2s_unbox, task_write1, task_write2)
            nbatch += 1
            n_r1_before = length(r1s) - n_reads
            n_r2_before = length(r2s) - n_reads

            if typeof(io1) <: IOStream  # not compressed
                length(in1bytes) == chunk_size1 || resize!(in1bytes, chunk_size1)
                length(in2bytes) == chunk_size2 || resize!(in2bytes, chunk_size2)
                (n_r1, n_r2, r1s, r2s, ncopied) = load_fqs_threads!(io1, io2, in1bytes, in2bytes, vr1s, vr2s, r1s, r2s, task_r1s_unbox, task_r2s_unbox; remove_first_n = n_reads, njobs = njobs, quality_offset = quality_offset)

                # it only get the sizes, did not change the sizes. Size changing is done in the "Read" part.
                chunk_size1, chunk_size2 = get_ideal_inbyte_sizes(in1bytes, in2bytes, n_r1, n_r2, n_r1_before, n_r2_before, max_chunk_size, chunk_size1, chunk_size2)

            else  # gziped
                total_n_bytes_read1 += length(in1bytes)  # will read INT in this batch
                total_n_bytes_read2 += length(in2bytes)  # will read INT in this batch
                will_eof1 = total_n_bytes_read1 >= uncompressed_size1
                will_eof2 = total_n_bytes_read2 >= uncompressed_size2
                (n_r1, n_r2, r1s, r2s, in1bytes_nremain, in2bytes_nremain, ncopied) = load_fqs_threads!(
                    io1, io2,
                    in1bytes, in2bytes, in1bytes_nremain, in2bytes_nremain,
                    vr1s, vr2s, r1s, r2s, task_r1s_unbox, task_r2s_unbox;
                    will_eof1 = will_eof1, will_eof2 = will_eof2,
                    in1bytes_resize_before_read = chunk_size1,
                    in2bytes_resize_before_read = chunk_size2,
                    remove_first_n = n_reads, quality_offset = quality_offset,
                    njobs = njobs
                );
            end

            n_reads = min(n_r1, n_r2)
            total_read_copied_in_loading += ncopied

            # processing reads
            isgoods = nbatch % 2 == 0 ? isgoods_even : isgoods_odd  # write task is still using the other one! 

            if do_pcr_dedup && nbatch == 1
                # size hint estimate
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

            Base.invokelatest(processing_reads_threads!, r1s, r2s, isgoods, n_reads)

            isgoods_in_range = view(isgoods, 1:n_reads)
            task_sum = Threads.@spawn sum(isgoods_in_range)

            task_r1s_unbox, task_r2s_unbox, task_write1, task_write2 = write_fqs_threads!(
                io1out::IO, io2out::IO,
                outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
                r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
                n_reads, isgoods_in_range, task_write1, task_write2)

            n_goods = fetch(task_sum)
            total_n_goods += n_goods
            total_n_reads += n_reads

            @info "Cycle $nbatch: read $n_reads/$total_n_reads pairs; wrote $n_goods/$total_n_goods pairs; (copied $ncopied/$total_read_copied_in_loading reads)"

            return task_r1s_unbox, task_r2s_unbox, task_write1, task_write2
        end

        task_r1s_unbox = Threads.@spawn 1
        task_r2s_unbox = Threads.@spawn 1
        task_write1 = Threads.@spawn 1
        task_write2 = Threads.@spawn 1

        # this has Precompiling
        if !eof(io1::IO) || !eof(io2::IO)
            task_r1s_unbox, task_r2s_unbox, task_write1, task_write2 = cycle_wrapper(task_r1s_unbox, task_r2s_unbox, task_write1, task_write2)
        end

        while !eof(io1::IO) || !eof(io2::IO)
            task_r1s_unbox, task_r2s_unbox, task_write1, task_write2 = cycle_wrapper(task_r1s_unbox, task_r2s_unbox, task_write1, task_write2)
        end
        
        @info "ATRIA COMPLETE" read1=outfile1 read2=outfile2
        with_logger(logger) do
            @info "ATRIA COMPLETE" read1=outfile1 read2=outfile2
        end

        wait(task_write1)
        wait(task_write2)

        time_read_processing = time() - time_read_processing


        #================== Close files ====================#
        time_post_processing = time()

        close(io1)
        close(io2)
        close(io1out)
        close(io2out)
        close(iolog)

        #================== Dedup stats ====================#
        if do_pcr_dedup 
            if do_pcr_dedup_count
                out_pcr_dedup_count = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.pcr_dedup_count.tsv", count=1))
                dup_count = write_pcr_dedup_count(out_pcr_dedup_count, dup_dict)
            else
                dup_count = get_dup_count(dup_dict)
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
            "input-read2" => file2
        )
        logjson["output-files"] = OrderedDict(
            "good-read1" => outfile1,
            "good-read2" => outfile2,
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
            "good-read-pairs" => total_n_goods,
            "total-read-pairs" => total_n_reads,
        )
        if do_pcr_dedup
            logjson["trimming-details"]["pcr-duplicate-pairs"] = dup_count
        end

        iologjson = open(outjson, "w+")
        JSON.print(iologjson, sort!(logjson), 4)
        close(iologjson)

    end

    return 0
end # func

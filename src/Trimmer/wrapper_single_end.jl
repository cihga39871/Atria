
# f_procs(x::String) = x == "-p" || x == "--procs"

function julia_wrapper_atria_single_end(ARGS::Vector{String}; exit_after_help = true)

    time_program_initializing = time()

    atria_version = "v3.0.2"

    args = parsing_args(ARGS; ver = atria_version, exit_after_help = exit_after_help)

    if args === nothing  # ARGS is ["-h"]
        return 0
    end
    args_range_test(args)

    outdir = args["output-dir"]

    #================== Parallel control ====================
    --procs 4onlyrun2 is a hidden feature for parallel computing...
    onlyrun2 means only run the second read fastqs.
    MECHANISM:
    when -p | --procs is specified, Atria will launch additional worker processes.
    the controller process (main process) will only manage worker processes,
    and will not run samples.
    the new workers will run additional Atria with a modified ARGS: -p XonlyrunN
    onlyrunN refresh `file_range` to N:N,
    so one worker will only run one sample each time. =#

    nfile = length(args["read1"])
    file_range = 1:nfile

    index_procs = findfirst(f_procs, ARGS)
    if index_procs != nothing  # argument --procs is specified
        nparallel = tryparse(Int64, args["procs"])
        if nparallel == nothing  # --procs might be 4onlyrun2
            if occursin(r"\d+onlyrun\d+", args["procs"])
                _filenum = parse(Int64, match(r"\d+onlyrun(\d+)$", args["procs"]).captures[1])
                file_range = _filenum:_filenum
            else
                @error "--procs INT must be positive integer" _module=:. _group=:. _id=:. _file="."
                exit(3)
            end
        else
            # start run this code as parallel.
            addprocs(nparallel)
            julia_command = Base.julia_cmd()

            function sub_procs_single_end(filenum::Int64, r1_filename::String, nfile::Int64, julia_command::Cmd)
                println("Atria: ($filenum/$nfile) $r1_filename")
                new_args = ARGS[1:end]
                # presumption: --procs is specified in ARGS
                new_args[index_procs+1] = new_args[index_procs+1] * "onlyrun$filenum"

                logs = joinpath(outdir, replace(basename(r1_filename), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.stdlog", count=1))
                try
                    run(pipeline(`$julia_command -e 'Atria = Base.loaded_modules[Base.PkgId(Base.UUID("226cbef3-b485-431c-85c2-d8bd8da14025"), "Atria")]; Atria.julia_main()' -- $new_args`, stdout=logs, stderr=logs))
                catch e
                    rethrow(e)
                end

                return 0
            end
            @eval @everywhere ARGS = $ARGS
            @eval @everywhere sub_procs_single_end = $sub_procs_single_end
            @eval @everywhere outdir = $outdir
            @info "Parallel mode: logs saved to $outdir/*.atria.(std)log(.json)"
            pmap(sub_procs_single_end,
                file_range,
                args["read1"],
                repeat(Int64[nfile], nfile),
                repeat(Cmd[julia_command], nfile)
            )
            file_range = 1:0  # no file will be processed in the main thread.
            return 0  # do not run Atria in the main thread. end of julia_wrapper_atria()
        end
    end


    #================== Arguments ====================#

    command                  = `$ARGS`
    nthread                  =  args["threads"               ]
    max_chunk_size           =  2 ^ args["log2-chunk-size"   ]
    # polyX
    min_poly_length          = args["poly-length"            ]
    poly_mismatch_per_16mer  = args["poly-mismatch-per-16mer"]
    # N
    max_N                    =  args["max-n"                 ]
    # adapter
    adapter1                 = LongDNASeq(args["adapter1"]) |> bitsafe!
    adapter1_seqheadset      = SeqHeadSet(adapter1)
    # NOTE: TruncSeq has some unknown accuracy problems.
    kmer_tolerance           = args["kmer-tolerance"          ]
    kmer_tolerance_consensus = args["kmer-tolerance-consensus"]
    trim_score               = args["trim-score-se"           ]
    tail_length              = args["tail-length"             ]
    # consensus
    # hard clip
    nclip_after        = args["clip-after"]
    nclip_front        = args["clip5"     ]
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
    do_polyG                 =  args["polyG"               ]
    do_polyT                 =  args["polyT"               ]
    do_polyA                 =  args["polyA"               ]
    do_polyC                 =  args["polyC"               ]
    do_length_filtration     = !args["no-length-filtration"]
    do_adapter_trimming      = !args["no-adapter-trim"     ]
    do_consensus_calling     =  false
    do_hard_clip_3_end       =  nclip_after > 0
    do_hard_clip_5_end       =  nclip_front > 0
    do_quality_trimming      = !args["no-quality-trim"     ]
    do_tail_n_trimming       = !args["no-tail-n-trim"      ]
    do_tail_low_qual_trimming=  do_polyG || do_polyT || do_polyA || do_polyC || do_adapter_trimming
    do_max_n_filtration      =  max_N > 0
    do_read_stats            =  args["stats"]
    do_complexity_filtration =  args["enable-complexity-filtration"]

    mkpath(outdir)


    #================== Main function and common variables ====================#


    #======= Check identifier (not applicable in single end)=======#
    CheckIdentifier = nothing

    #======= PolyX =======#
    PolyG = do_polyG ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_G, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyT = do_polyT ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_T, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyA = do_polyA ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_A, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    PolyC = do_polyC ? quote
        polyX_idx, polyX_length = polyX_tail_scan(DNA_C, r1, $poly_mismatch_per_16mer)
        if polyX_length >= $min_poly_length
            polyX_idx -= 1
            tail_trim!(r1, polyX_idx)
        end
    end : nothing

    #======= Complexity =======#
    ComplexityFilter = do_complexity_filtration ? quote
        if seq_complexity(r1) < $min_complexity
            return false
        end
    end : nothing

    #======= Hard clips =======#
    HardClip3End = do_hard_clip_3_end ? quote
        tail_trim!(r1::FqRecord, $nclip_after)
    end : nothing

    HardClip5End = do_hard_clip_5_end ? quote
        front_trim!(r1::FqRecord, $nclip_front::Int64)
    end : nothing


    #======= Trim N from the tail =======#
    TailNTrim = do_tail_low_qual_trimming ? quote
        tail_low_qual_trim!(r1::FqRecord)
    end : do_tail_n_trimming ? quote
        tail_N_trim!(r1::FqRecord)
    end : nothing


    #======= If reads are too short, remove =======#
    LengthFilter = do_length_filtration ? quote
        if !isinreadlength!(r1::FqRecord, $length_range)
            # global count_reads_pair_invalid_length[thread_id] += 1
            return false
        end
    end : nothing


    #======= If reads have so many N, remove =======#
    MaxNFilter = do_max_n_filtration ? quote
        if !isnotmuchN!(r1::FqRecord, $max_N)
            return false
        end
    end : nothing


    #======= Quality trimming =======#
    QualityTrim = do_quality_trimming ? quote
        quality_index1 = qualitymatch(r1::FqRecord, $quality_single::UInt8, $quality_threshold::UInt64, $quality_kmer::Int64)
        if quality_index1 != -1  # -1 means no need to do quality trimming
            # global count_r1_quality_trimming[thread_id] += 1
            tail_trim!(r1::FqRecord, quality_index1::Int64)
        end
    end : nothing


    #======= adapter trimming =======#
    AdapterTrim = do_adapter_trimming ? quote

        r1_adapter_pos, r1_adapter_nmatch = bitwise_scan($adapter1_seqheadset, r1.seq, 1, $kmer_tolerance)

        r1_insert_size = r1_adapter_pos - 1

        r1_adapter_prob = probmean(r1, r1_adapter_pos, r1_adapter_pos + 15)

        match_length = length(r1.seq) - r1_insert_size
        if match_length >= 16
            r1_adapter_score = @fastmath r1_adapter_nmatch * r1_adapter_prob
        else
            r1_adapter_score = @fastmath(
                16 * r1_adapter_nmatch * r1_adapter_prob / match_length
            )
        end

        if r1_adapter_score > $trim_score
            # r1_insert_size can be -1
            tail_trim!(r1::FqRecord, r1_insert_size < 0 ? 0 : r1_insert_size)
        end
    end : nothing


    @eval function read_processing!(r1::FqRecord, thread_id::Int)::Bool
        $CheckIdentifier
        $PolyG
        $PolyT
        $PolyA
        $PolyC
        $LengthFilter
        $AdapterTrim
        $HardClip3End
        $HardClip5End
        $QualityTrim
        $TailNTrim
        $MaxNFilter
        $LengthFilter
        $ComplexityFilter
        return true
    end


    in1bytes = Vector{UInt8}(undef, max_chunk_size)

    # number of jobs to boxing FqRecord from UInt8 Vector
    njobs = nthread * 5
    vr1s = ntuple(_ -> Vector{FqRecord}(), njobs)

    r1s = Vector{FqRecord}()

    isgoods = Vector{Bool}()

    outr1s = Vector{Vector{UInt8}}()



    time_program_initializing = time() - time_program_initializing

    #================== Iteration for files ====================#
    for filenum in file_range
        # filenum = 1
        time_file_initializing = time()


        #===== file names =====#

        file1 = args["read1"][filenum]

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

        outfile1 = joinpath(outdir, replace(basename(file1), r"(fastq$|fq$|[^.]*)(\.gz)?$"i => s"atria.\1", count=1))
        if outcompress == "GZIP"
            outfile1 *= ".gz"
        elseif outcompress == "BZIP2"
            outfile1 *= ".bz2"
        end

        outlog = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.log", count=1))
        outjson = joinpath(outdir, replace(basename(file1), r"fastq$|fq$|[^.]*(\.gz)?$"i => "atria.log.json", count=1))


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

        @info "ATRIA ARGUMENTS" command
        @info "ATRIA OUTPUT FILES" read1=outfile1
        @info("ATRIA TRIMMERS AND FILTERS",
            tail_polyG_trimming = do_polyG,
            tail_polyT_trimming = do_polyT,
            tail_polyA_trimming = do_polyA,
            tail_polyC_trimming = do_polyC,
            adapter_trimming = do_adapter_trimming,
            consensus_calling = do_consensus_calling,
            hard_clip_3_end = do_hard_clip_3_end,
            hard_clip_5_end = do_hard_clip_5_end,
            quality_trimming = do_quality_trimming,
            tail_N_trimming = do_tail_n_trimming,
            max_N_filtering = do_max_n_filtration,
            length_filtering = do_length_filtration,
            complexity_filtering = do_complexity_filtration)

        with_logger(logger) do
            @info "ATRIA ARGUMENTS" command
            @info "ATRIA OUTPUT FILES" read1=outfile1
            @info("ATRIA TRIMMERS AND FILTERS",
                tail_polyG_trimming = do_polyG,
                tail_polyT_trimming = do_polyT,
                tail_polyA_trimming = do_polyA,
                tail_polyC_trimming = do_polyC,
                adapter_trimming = do_adapter_trimming,
                consensus_calling = do_consensus_calling,
                hard_clip_3_end = do_hard_clip_3_end,
                hard_clip_5_end = do_hard_clip_5_end,
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
        total_n_bytes_read1 = 0
        in1bytes_nremain = 0

        #================== File processing ====================#
        time_file_initializing = time() - time_file_initializing
        time_read_processing = time()

        # the first cycle to generate compiled code?
        function cycle_wrapper_single_end()
            nbatch += 1
            n_r1_before = length(r1s) - n_reads

            if typeof(io1) <: IOStream  # not compressed
                (n_r1, r1s, ncopied) = load_fqs_threads!(io1, in1bytes, vr1s, r1s; remove_first_n = n_reads, njobs=njobs)
            else  # gziped
                total_n_bytes_read1 += length(in1bytes)  # will read INT in this batch
                will_eof1 = false
                (n_r1, r1s, in1bytes_nremain, ncopied) = load_fqs_threads!(
                    io1,
                    in1bytes,
                    in1bytes_nremain,
                    vr1s,
                    r1s;
                    will_eof1 = will_eof1,
                    remove_first_n = n_reads,
                    njobs = njobs
                )
            end

            n_reads = n_r1
            total_read_copied_in_loading += ncopied

            processing_reads_threads!(r1s, isgoods, n_reads)

            #= debug
            for i in 1:length(r1s)
                @info i
                read_processing!(r1s[i], 1)
            end
            continue
            =#

            isgoods_in_range = view(isgoods, 1:n_reads)
            task_sum = Threads.@spawn sum(isgoods_in_range)
            write_fqs_threads!(
                io1out::IO,
                outr1s::Vector{Vector{UInt8}},
                r1s::Vector{FqRecord},
                n_reads::Int, isgoods_in_range)

            n_goods = fetch(task_sum)
            total_n_goods += n_goods
            total_n_reads += n_reads

            # @info "Cycle $nbatch: processed $n_reads reads ($total_n_reads in total), in which $n_goods passed filtration ($total_n_goods in total). ($ncopied/$total_read_copied_in_loading reads copied)"
            @info "Cycle $nbatch: read $n_reads/$total_n_reads pairs; wrote $n_goods/$total_n_goods; (copied $ncopied/$total_read_copied_in_loading)"
        end

        while !eof(io1::IO)
            cycle_wrapper_single_end()
        end

        time_read_processing = time() - time_read_processing


        #================== Close files ====================#
        time_post_processing = time()

        close(io1)
        close(io1out)
        close(iolog)

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
            "good-read" => total_n_goods,
            "total-read" => total_n_reads,
        )

        iologjson = open(outjson, "w+")
        JSON.print(iologjson, sort!(logjson), 4)
        close(iologjson)

    end

    return 0
end # func

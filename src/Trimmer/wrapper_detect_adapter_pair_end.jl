
# f_procs(x::String) = x == "-p" || x == "--procs"

function julia_wrapper_detect_adapter_pe(ARGS::Vector{String}; exit_after_help = true)

    time_program_initializing = time()

    args = parsing_args(ARGS; exit_after_help = exit_after_help)

    if args === nothing  # ARGS is ["-h"]
        return 0
    end
    args_range_test(args)
    
    nthread = args["threads"]

    nfile = length(args["read1"])
    file_range = 1:nfile

    #================== Arguments ====================#

    max_chunk_size           =  2 ^ args["log2-chunk-size"]

    # NOTE: TruncSeq has some unknown accuracy problems.
    kmer_tolerance           = args["kmer-tolerance"          ]
    kmer_n_match             = args["kmer-n-match"            ]

    # quality
    quality_offset     = Trimmer.get_quality_offset(args["quality-format"])


    mkpath(outdir)


    #================== Main function and common variables ====================#

    in1bytes = Vector{UInt8}(undef, max_chunk_size)
    in2bytes = Vector{UInt8}(undef, max_chunk_size)

    # number of jobs to boxing FqRecord from UInt8 Vector
    njobs = nthread * 10
    vr1s = ntuple(_ -> Vector{FqRecord}(), njobs)
    vr2s = ntuple(_ -> Vector{FqRecord}(), njobs)

    r1s = Vector{FqRecord}()
    r2s = Vector{FqRecord}()


    time_program_initializing = time() - time_program_initializing

    adapter_detection_summary = init_adapter_detection_summary()
    #================== Iteration for paired files ====================#
    for filenum in file_range
        # filenum = 1
        time_file_initializing = time()


        #===== file names =====#

        file1 = args["read1"][filenum]
        file2 = args["read2"][filenum]

        # check whether this sample is processed before

        isingzip = occursin(r"\.gz$"i, file1)
        isinbzip2 = occursin(r"\.bz2$"i, file1)


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

        n_reads = 0
        n_r1 = 0
        n_r2 = 0
        nbatch = 0
        total_read_copied_in_loading = 0
        total_n_bytes_read1 = 0
        total_n_bytes_read2 = 0
        in1bytes_nremain = 0
        in2bytes_nremain = 0

        #================== File processing ====================#


        # the first cycle to generate compiled code?
        function cycle_wrapper()
            nbatch += 1
            n_r1_before = length(r1s) - n_reads
            n_r2_before = length(r2s) - n_reads

            if typeof(io1) <: IOStream  # not compressed
                length(in1bytes) == chunk_size1 || resize!(in1bytes, chunk_size1)
                length(in2bytes) == chunk_size2 || resize!(in2bytes, chunk_size2)
                (n_r1, n_r2, r1s, r2s, ncopied) = load_fqs_threads!(io1, io2, in1bytes, in2bytes, vr1s, vr2s, r1s, r2s; remove_first_n = n_reads, njobs = njobs, quality_offset = quality_offset)
            else  # gziped
                total_n_bytes_read1 += length(in1bytes)  # will read INT in this batch
                total_n_bytes_read2 += length(in2bytes)  # will read INT in this batch
                will_eof1 = total_n_bytes_read1 >= uncompressed_size1
                will_eof2 = total_n_bytes_read2 >= uncompressed_size2
                (n_r1, n_r2, r1s, r2s, in1bytes_nremain, in2bytes_nremain, ncopied) = load_fqs_threads!(
                    io1, io2,
                    in1bytes, in2bytes, in1bytes_nremain, in2bytes_nremain,
                    vr1s, vr2s, r1s, r2s;
                    will_eof1 = will_eof1, will_eof2 = will_eof2,
                    in1bytes_resize_before_read = chunk_size1,
                    in2bytes_resize_before_read = chunk_size2,
                    remove_first_n = n_reads, quality_offset = quality_offset,
                    njobs = njobs
                )
            end

            n_reads = min(n_r1, n_r2)
            total_read_copied_in_loading += ncopied

            # it only get the sizes, did not change the sizes. Size changing is done in the "Read" part.
            chunk_size1, chunk_size2 = adjust_inbyte_sizes(in1bytes, in2bytes, n_r1, n_r2, n_r1_before, n_r2_before, max_chunk_size, chunk_size1, chunk_size2)

            # check_fq_ids(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, n_reads::Int)::nothing

            # processing reads
            r1_stats, r2_stats = check_pe_match(r1s, r2s; kmer_tolerance = kmer_tolerance + 1, kmer_n_match = kmer_n_match, occurance = 0.0004)

            show_paired_adapter_result(file1, r1_stats, n_reads)
            show_paired_adapter_result(file2, r2_stats, n_reads)
            push_adapter_detection_summary!(adapter_detection_summary, file1, r1_stats, file2, r2_stats)
        end

        cycle_wrapper()

        #================== Close files ====================#

        close(io1)
        close(io2)
    end

    timestamp = replace(string(now()), r"[T:\.]" => "-")
    adapter_detection_summary_file = joinpath(outdir, "atria_adapter_detect_summary.$timestamp.txt")
    CSV.write(adapter_detection_summary_file, adapter_detection_summary, delim = '\t')
    println("""
    _________________________________

    Summary of detected adapters is saved to $adapter_detection_summary_file
    
    _________________________________

    Paired-end Adapter Detection Note: 
    
    Atria detects adapter sequences using paired-end information. Adapter sequences are truncated to 16-bp, which are accurate enough for trimming. From experiments of many popular trimmers, increasing adapter length from 16 to 33 does not increase accuracy (Figure 4C of https://doi.org/10.46471/gigabyte.31).

    Adapter detection is the last choice because its accuracy is highly based on your data. If your data have been trimmed, the remaining adapters may not enough for accurate guessing. We suggest to use adapter detection only when you cannot find the actual adapter sequence.

    Besides, Atria does not automatically trim auto-detected adapters. It is your responsibility to check whether the detected adapters are real.
    
    Those rules can be used to check the adapter results: 
    
    (1) An Illumina sequence file only has ONE adapter sequence. 
    
    (2) In the same batch of NGS experiments, all R1 samples should have the SAME adapter sequence, so do R2 samples. The most prevalent adapters of R1 and R2 might be true for all your data.
    _________________________________

    Summary of detected adapters is saved to $adapter_detection_summary_file
    _________________________________

    """)

    return 0
end # func

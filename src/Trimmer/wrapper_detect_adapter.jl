
# f_procs(x::String) = x == "-p" || x == "--procs"

function julia_wrapper_detect_adapter(ARGS::Vector{String}; exit_after_help = true)

    args = parsing_args(ARGS; exit_after_help = exit_after_help)

    if args === nothing  # ARGS is ["-h"]
        return 0
    end
    args_range_test(args)


    #================== Arguments ====================#

    nthread                  =  args["threads"            ]
    max_chunk_size           =  2 ^ args["log2-chunk-size"]

    #================== Main function and common variables ====================#

    in1bytes = Vector{UInt8}(undef, max_chunk_size)

    # number of jobs to boxing FqRecord from UInt8 Vector
    njobs = nthread * 5
    vr1s = ntuple(_ -> Vector{FqRecord}(), njobs)

    r1s = Vector{FqRecord}()

    #================== Iteration for files ====================#

    append!(args["read1"], args["read2"])

    for file1 in args["read1"]
        # filenum = 1
        time_file_initializing = time()

        #===== file names =====#

        isingzip = occursin(r"\.gz$"i, file1)
        isinbzip2 = occursin(r"\.bz2$"i, file1)

        #===== file IO =====#

        if isingzip
            io1 = open(`pigz -p$nthread -cd $file1`, write=false)
        elseif isinbzip2
            io1 = open(`pbzip2 -p$nthread -cd $file1`, write=false)
        else
            io1 = open(file1, "r")
        end


        #================== Renew variables for read processing ====================#

        # clear common variables
        empty!(r1s)

        n_reads = 0
        n_r1 = 0
        in1bytes_nremain = 0

        #================== File processing ====================#
        time_file_initializing = time() - time_file_initializing
        time_read_processing = time()

        # the first cycle to generate compiled code?
        function cycle_wrapper_detect_adapter()
            n_r1_before = length(r1s) - n_reads

            if typeof(io1) <: IOStream  # not compressed
                (n_r1, r1s, ncopied) = load_fqs_threads!(io1, in1bytes, vr1s, r1s; remove_first_n = n_reads, njobs=njobs)
            else  # gziped
                (n_r1, r1s, in1bytes_nremain, ncopied) = load_fqs_threads!(
                    io1,
                    in1bytes,
                    in1bytes_nremain,
                    vr1s,
                    r1s;
                    remove_first_n = n_reads,
                    njobs = njobs
                )
            end

            top3 = detect_adapter_threads!(r1s)

            @info "adapter detection in $file1 ($n_r1 reads)" top_3_adapters=top3
        end

        cycle_wrapper_detect_adapter()

        #================== Close files ====================#

        close(io1)
    end

    return 0
end # func

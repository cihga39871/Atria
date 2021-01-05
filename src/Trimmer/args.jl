function parsing_args(args::Vector; ver::String="x.x.x", exit_after_help = true)
    settings = ArgParseSettings(
        prog = "atria",
        description = "Atria $ver",
        version = ver, add_version = true, exit_after_help = exit_after_help,
        epilog = """Author: Jiacheng (Eric) Chuan\n"""
    )
    @add_arg_table! settings begin
        "--threads", "-t"
            help = "use INT threads to process one sample (multi-threading parallel)."
            metavar = "INT"
            default = Threads.nthreads()
            arg_type = Int
        "--procs", "-p"
            # hidden feature: --procs 8onlyrun2 will only run the second pair-end FASTQ!
            help = "process at most INT samples at the same time (multi-core parallel)"
            metavar = "INT"
            default = "1"
            arg_type = String
        "--log2-chunk-size"
            help = "read at most 2^INDEX bits each time. Suggest to process 200,000 reads each time. Reduce INDEX to lower the memory usage."
            metavar = "INDEX"
            default = 26
            arg_type = Int
    end
    add_arg_group!(settings, "input/output: input read 1 and read 2 should be in the same order")
    @add_arg_table! settings begin
        "--read1", "-r"
            help = "input read 1 fastq file(s) (if unix & gzip, require `pigz` in environment)"
            nargs = '+'
            metavar = "R1-FASTQ"
            required = true
        "--read2", "-R"
            help = "input read 2 fastq file(s) (paired with R1-FASTQ)"
            nargs = '+'
            metavar = "R2-FASTQ"
            required = true
        "--output-dir", "-o"
            help = "store output files and stats to PATH"
            metavar = "PATH"
            default = pwd()
        "--gzip", "-g"
            help = "gzip the output fastq files (if unix & gzip, require `pigz` in environment)"
            metavar = "AUTO|YES|NO"
            default = "AUTO"
        "--check-identifier"
            help = "check whether the identifiers of r1 and r2 are the same"
            action = :store_true
    end

    add_arg_group!(settings, "adapter trimming")
    @add_arg_table! settings begin
        "--no-adapter-trim"
            help = "disable adapter and pair-end trimming"
            action = :store_true
        "--adapter1", "-a"
            help = "read 1 adapter"
            metavar = "SEQ"
            default = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        "--adapter2", "-A"
            help = "read 2 adapter"
            metavar = "SEQ"
            default = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        "--kmer-tolerance", "-T"
            help = "# of mismatch allowed in 16-mers adapter and pair-end matching"
            default = 2
            metavar = "INT"
            arg_type = Int64
        "--pe-adapter-diff", "-d"
            help = "number of bases allowed when disconcordance found between adapter and pair-end search"
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--r1-r2-diff", "-D"
            help = "number of bases allowed when the insert sizes of r1 and r2 are different"
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--r1-r2-score-diff", "-s"
            help = "if the score difference between the 16-mers of r1 and r2 is greater than FLOAT, the insert size of the low-score read will be adapted to the high-score read"
            default = 3.0
            metavar = "FLOAT"
            arg_type = Float64
        "--trim-score", "-S"
            help = "if the final score of the read pair is greater than FLOAT, the reads will be trimmed."
            default = 8.0
            metavar = "FLOAT"
            arg_type = Float64
        "--tail-length", "-l"
            help = "if the adapter is in the tail region, and insert size of pe match is smaller than this region, do not trim the read."
            default = 8
            metavar = "INT"
            arg_type = Int64
        "--stats"
            help = "write stats to description lines of r2 reads."
            action = :store_true
    end

    add_arg_group!(settings, "consensus/merging in adapter trimming")
    @add_arg_table! settings begin
        "--no-consensus"
            help = "disable generating consensus paired reads. If adapter trimming is disabled, consensus calling is not performed even the flag is not set."
            action = :store_true
        "--kmer-tolerance-consensus"
            help = "# of mismatch allowed in 16-mers matching in consensus calling"
            default = 10
            metavar = "INT"
            arg_type = Int64
        "--min-ratio-mismatch"
            help = "if the ratio of mismatch of the overlapped region is less than FLOAT, skip consensus calling."
            default = 0.28
            metavar = "FLOAT"
            arg_type = Float64
        "--overlap-score"
            help = "if no adapter was found, scan the tails of the paired reads. Then, if the maximum score of the overlapped 16-mers are less than FLOAT, skip consensus calling for the read pair. If adapters were found, this step is ignored."
            default = 0.0
            metavar = "FLOAT"
            arg_type = Float64
        "--prob-diff"
            help = "when doing consensus calling, if the bases were not complementary, the base with the higher quality probability is selected unless the quality probability difference are less than FLOAT"
            default = 0.0
            metavar = "FLOAT"
            arg_type = Float64
    end

    add_arg_group!(settings, "hard clipping: trim a fixed length (after adapter trimming)")
    @add_arg_table! settings begin
        "--clip-after", "-C"
            help = "hard clip the 3' tails to contain only INT bases. 0 to disable."
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--clip5", "-c"
            help = "remove the first INT bases from 5' end."
            default = 0
            metavar = "INT"
            arg_type = Int64
    end

    add_arg_group!(settings, "quality trimming: trim the tail when the average quality of bases in\na sliding window is low (after hard clipping)")
    @add_arg_table! settings begin
        "--no-quality-trim"
            help = "skip quality trimming"
            action = :store_true
        "--quality-score", "-q"
            help ="threshold of quality score; 0 means turn off quality trimming"
            default = 20
            metavar = "INT"
            arg_type = Int64
        "--quality-kmer"
            help = "trim the tail once found the average quality of bases in a sliding window is low"
            default = 5
            metavar = "INT"
            arg_type = Int64
        "--quality-format"
            help = "the format of the quality score (Illumina1.3, Illumina1.8, Sanger, Illumina1.5, Solexa); or the ASCII number when quality score == 0"
            default = "33"
            metavar = "FORMAT"
            arg_type = String
    end

    add_arg_group!(settings, "N trimming (after quality trimming)")
    @add_arg_table! settings begin
        "--no-tail-n-trim"
            help = "disable removing NNNNN tail."
            action = :store_true
        "--max-n", "-n"
            help = "# N allowed in each read; N tails not included if --no-tail-n-trim; INT<0 to disable"
            metavar = "INT"
            default = 15
            arg_type = Int64
    end

    add_arg_group!(settings, "length filter (after N trimming)")
    @add_arg_table! settings begin
        "--no-length-filtration"
            help = "disable length filtration"
            action = :store_true
        "--length-range"
            help = "length range of good reads; format is min:max"
            default = "50:500"
            metavar = "INT:INT"
            arg_type = String
    end
    return parse_args(args, settings)
end

function isdna(x::String)
    all(map(d -> uppercase(d) in ['A', 'G', 'C', 'T', 'N'], collect(x)))
end

function args_range_test(args::Dict{String,Any}; test_only::Bool=false)
    ispass = true

    @static if !Sys.iswindows()
        try
            x = readchomp(`pigz --version`)
        catch
            @error "Dependency missing: pigz not found." _module=nothing _group=nothing _id=nothing _file=nothing
            ispass = false
        end
    end

    if args["threads"] != Threads.nthreads()
        if args["threads"] > Sys.CPU_THREADS
            @error "Failed to set --threads INT (INT=$(args["threads"])): greater than available CPU threads." _module=nothing _group=nothing _id=nothing _file=nothing
        else
            @static if Sys.iswindows()
                @error "Cannot set number of threads using --threads INT on Windows. Please add JULIA_NUM_THREADS=INT to your environment variables to enable multi-threading." _module=nothing _group=nothing _id=nothing _file=nothing
            else
                @error "Failed to set --threads INT (INT=$(args["threads"])). Was Atria bash script changed? If you use IOS, you may post an issue at Github Atria." _module=nothing _group=nothing _id=nothing _file=nothing
            end
        end
        ispass = false
    end

    if !occursin(r"^[1-9][0-9]*$|^[1-9][0-9]*onlyrun\d+$", args["procs"])
        @error "--procs INT must be positive integer" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if args["log2-chunk-size"] < 20 || args["log2-chunk-size"] > 30
        @error "--log2-chunk-size INDEX is not in a suggested range." _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # input
    if length(args["read1"]) != length(args["read2"])
        @error "--read1 R1-FASTQ... and --read2 R2-FASTQ... must be paired respectively" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    for file in [args["read1"]; args["read2"]]
        if !isfile(file)
            @error "input fastq file not valid" FILE=file _module=nothing _group=nothing _id=nothing _file=nothing
            ispass = false
        end
    end
    if !(uppercase(args["gzip"]) in ["AUTO", "YES", "NO"])
        @error "--gzip -g OPTION invalid; options should be one of AUTO, YES, NO." _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    try
        mkpath(args["output-dir"])
        tmp_file = joinpath(args["output-dir"], basename(tempname()))
        touch(tmp_file)
        rm(tmp_file)
    catch
        @error "--output-dir -o PATH is not writable. Do you have write permission?" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # N trimming
    # if args["max-n"] < 0
    #     @warn "--max-n -n INT less than 0: disable max-n filtration" INT=args["max-n"] _module=nothing _group=nothing _id=nothing _file=nothing
    # end

    # adapter triming
    if !isdna(args["adapter1"])
        @error "--adapter1 -a SEQ contains invalid bases (not ACGTN)" SEQ=args["adapter1"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if !isdna(args["adapter2"])
        @error "--adapter2 -A SEQ contains invalid bases (not ACGTN)" SEQ=args["adapter2"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if args["kmer-tolerance"] >= 16 || args["kmer-tolerance"] < 0
        @error "--kmer-tolerance -T INT not in range of 0:15" INT=args["kmer-tolerance"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["pe-adapter-diff"] < 0
        @error "--pe-adapter-diff -d INT < 0 is not valid" INT=args["pe-adapter-diff"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["r1-r2-diff"] < 0
        @error "--r1-r2-diff -D INT < 0 is not valid" INT=args["r1-r2-diff"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["r1-r2-score-diff"] < 0.0
        @error "--r1-r2-score-diff -s FLOAT < 0 is not valid" INT=args["r1-r2-score-diff"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["trim-score"] < 0.0
        @error "--trim-score -S FLOAT < 0 is not valid" INT=args["trim-score"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["tail-length"] < 1
        @warn "--tail-length -l INT < 1: the result might be unexpected" INT=args["tail-length"] _module=nothing _group=nothing _id=nothing _file=nothing
    end

    # consensus
    if args["kmer-tolerance-consensus"] >= 16 || args["kmer-tolerance-consensus"] < 0
        @error "--kmer-tolerance-consensus -T INT not in range of 0:15" INT=args["kmer-tolerance"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["overlap-score"] < 0.0 || args["overlap-score"] > 16.0
        @error "--overlap-score FLOAT not in 0-16" INT=args["overlap-score"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["min-ratio-mismatch"] < 0.0 || args["min-ratio-mismatch"] > 1.0
        @error "--min-ratio-mismatch FLOAT not in 0-1" INT=args["min-ratio-mismatch"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["prob-diff"] < 0.0 || args["prob-diff"] > 1.0
        @error "--prob-diff FLOAT not in 0-1" INT=args["prob-diff"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # hard clip
    if 0 < args["clip-after"] < 30
        @warn "--clip-after -C INT is very small. Did you mis-understand this option? It does not clip INT bases from 3' end, but only keep the first INT bases of each read. For example, if you want to clip 5 bases from 3' end, and your read length is 101, you should use `--clip-after 96` or `-C 96`." INT=args["clip-after"] _module=nothing _group=nothing _id=nothing _file=nothing
    end
    if args["clip5"] > 70
        @warn "--clip5 -c INT is very large. Did you mis-understand this option? It does not keep the last INT bases of each read, but clip the first INT bases from 5' end." INT=args["clip5"] _module=nothing _group=nothing _id=nothing _file=nothing
    end

    # quality trimming
    # quality-score can be anything
    if args["quality-kmer"] < 0
        @error "--quality-kmer INT invalid" INT=args["quality-kmer"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if !occursin(r"^\d+$|^ILLUMINA1\.[358]$|^SOLEXA$|^SANGER$", uppercase(args["quality-format"]))
        @error "--quality-format FORMAT invalid; FORMAT shoule be the format of the quality score (Illumina1.3, Illumina1.8, Sanger, Illumina1.5, Solexa) or the ASCII number when quality score == 0" FORMAT=args["quality-format"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # filter
    if !occursin(r"^\d+\:\d+$", args["length-range"])
        @error "--length-range RANGE invalid; RANGE should be MIN:MAX where MIN and MAX are positive integer" RANGE=args["length-range"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if test_only
        return nothing
    end
    # output in the last
    try
        mkpath(args["output-dir"])
        testfile = joinpath(args["output-dir"], basename(tempname()))
        touch(testfile)
        rm(testfile)
    catch
        @error "output directory not writable" PATH=args["output-dir"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if !ispass
        exit(3)
    end
    nothing
end

function get_quality_offset(quality_format_string::String)
    if occursin(r"^\d+$", quality_format_string)
        parse(Int64, quality_format_string)
    else
        fmt = uppercase(quality_format_string)
        if fmt in ["ILLUMINA1.3", "ILLUMINA1.5", "SOLEXA"]
            64
        elseif fmt in ["SANGER", "ILLUMINA1.8"]
            33
        else
            @error "--quality-format FORMAT invalid; FORMAT shoule be the format of the quality score (Illumina1.3, Illumina1.8, Sanger, Illumina1.5, Solexa) or the ASCII number when quality score == 0" FORMAT=quality_format_string  _module=nothing _group=nothing _id=nothing _file=nothing
            exit(1)
        end
    end
end

function get_length_range(length_range_string::String)
    r = parse.(Int, split(length_range_string, ':'))
    r[1]:r[2]
end

function parsing_args(args::Vector; exit_after_help = true)
    settings = ArgParseSettings(
        prog = "atria",
        description = "Atria $atria_version",
        version = atria_version, add_version = true, exit_after_help = exit_after_help,
        preformatted_epilog = true,
        epilog = """

        ----------
        Process names for --order | -O:
            DefaultOrder    = [CheckIdentifier, PolyG, PolyT, PolyA, PolyC, LengthFilter, AdapterTrim, HardClip3EndR1, HardClip3EndR2, HardClip5EndR1, HardClip5EndR2, QualityTrim, TailNTrim, MaxNFilter, LengthFilter, ComplexityFilter],
            CheckIdentifier ,
            PolyX           = [PolyG, PolyT, PolyA, PolyC],
            PolyG           ,
            PolyT           ,
            PolyA           ,
            PolyC           ,
            AdapterTrim     ,
            HardClip3End    = [HardClip3EndR1, HardClip3EndR2],
            HardClip3EndR1  ,
            HardClip3EndR2  ,
            HardClip5End    = [HardClip5EndR1, HardClip5EndR2],
            HardClip5EndR1  ,
            HardClip5EndR2  ,
            QualityTrim     ,
            TailNTrim       ,
            MaxNFilter      ,
            LengthFilter    ,
            ComplexityFilter.

        ----------

        Jiacheng Chuan, Aiguo Zhou, Lawrence Richard Hale, Miao He, Xiang Li, Atria: an ultra-fast and accurate trimmer for adapter and quality trimming, Gigabyte, 1, 2021 https://doi.org/10.46471/gigabyte.31\n
        """
    )
    @add_arg_table! settings begin
        "--threads", "-t"
            help = "use INT threads to process one sample (multi-threading parallel)."
            metavar = "INT"
            default = Threads.nthreads()
            arg_type = Int
        "--log2-chunk-size"
            help = "read at most 2^INDEX bits each time. Suggest to process 200,000 reads each time. Reduce INDEX to lower the memory usage."
            metavar = "INDEX"
            default = 26
            arg_type = Int
        "--force", "-f"
            help = "force to analyze all samples; not skip completed ones"
            action = :store_true
    end
    add_arg_group!(settings, "input/output: input read 1 and read 2 should be in the same order")
    @add_arg_table! settings begin
        "--read1", "-r"
            help = "input read 1 fastq file(s), or single-end fastq files"
            nargs = '+'
            metavar = "R1-FASTQ"
            required = true
        "--read2", "-R"
            help = "input read 2 fastq file(s) (paired with R1-FASTQ)"
            nargs = '*'
            metavar = "R2-FASTQ"
        "--output-dir", "-o"
            help = "store output files and stats to PATH"
            metavar = "PATH"
            default = pwd()
        "--compress", "-g"
            help = "compression methods for output files (AUTO: same as input, NO: no compression, GZ|GZIP: gzip with `pigz`, BZ2|BZIP2: bzip2 with `pbzip2`)"
            metavar = "AUTO|NO|GZ|GZIP|BZ2|BZIP2"
            default = "AUTO"
        "--check-identifier"
            help = "check whether the identifiers of r1 and r2 are the same"
            action = :store_true
        "--detect-adapter"
            help = "detect possible adapters for each sample only"
            action = :store_true
    end

    add_arg_group!(settings, "processing order")
    @add_arg_table! settings begin
        "--order", "-O"
            help = "order of trimming and filtration processing methods. Unlisted process will not be done. See epilog for process names"
            metavar = "PROCESS"
            default = String["DefaultOrder"]
    end

    add_arg_group!(settings, "poly X tail trimming")
    @add_arg_table! settings begin
        "--polyG"
            help = "enable trimming poly G tails"
            action = :store_true
        "--polyT"
            help = "enable trimming poly T tails"
            action = :store_true
        "--polyA"
            help = "enable trimming poly A tails"
            action = :store_true
        "--polyC"
            help = "enable trimming poly C tails"
            action = :store_true
        "--poly-length"
            help = "the minimum length of poly X"
            default = 10
            arg_type = Int64
        "--poly-mismatch-per-16mer"
            help = "the number of mismatch allowed in 16 mer poly X"
            default = 2
            metavar = "INT"
            arg_type = Int64
    end

    add_arg_group!(settings, "adapter trimming")
    @add_arg_table! settings begin
        "--no-adapter-trim"
            help = "disable adapter and pair-end trimming"
            action = :store_true
        "--adapter1", "-a"
            help = "read 1 adapter(s) appended to insert DNA at 3' end"
            nargs = '+'
            metavar = "SEQ"
            default = String["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"]
        "--adapter2", "-A"
            help = "read 2 adapter(s) appended to insert DNA at 3' end"
            nargs = '+'
            metavar = "SEQ"
            default = String["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]
        "--kmer-tolerance", "-T"
            help = "# of mismatch allowed in 16-mers adapter and pair-end matching"
            default = 2
            metavar = "INT"
            arg_type = Int64
        "--pe-adapter-diff", "-d"
            help = "(FOR PAIRED END) number of bases allowed when disconcordance found between adapter and pair-end search"
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--r1-r2-diff", "-D"
            help = "(FOR PAIRED END) number of bases allowed when the insert sizes of r1 and r2 are different"
            default = 0 #1
            metavar = "INT"
            arg_type = Int64
        "--kmer-n-match", "-s"
            help = "(FOR PAIRED END) if n base matched [0-16] is less than INT, loosen matches will be made based on the match with the highest n base match"
            default = 9
            metavar = "INT"
            arg_type = Int64
        "--trim-score-pe"
            help = "(FOR PAIRED END) if final score [0-32] of read pair is greater than FLOAT, the reads will be trimmed."
            default = 10.0
            metavar = "FLOAT"
            arg_type = Float64
        "--trim-score-se"
            help = "(FOR SINGLE END) if final score [0-16] of read is greater than FLOAT, the reads will be trimmed."
            default = 10.0
            metavar = "FLOAT"
            arg_type = Float64
        "--tail-length", "-l"
            help = "(FOR PAIRED END) if the adapter is in the tail region, and insert size of pe match is smaller than this region, do not trim the read."
            default = 12
            metavar = "INT"
            arg_type = Int64
        "--stats"
            help = "(DEV ONLY) write stats to description lines of r2 reads."
            action = :store_true
    end

    add_arg_group!(settings, "consensus/merging in adapter trimming (FOR PAIRED END)")
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

    # add_arg_group!(settings, "primer trimming (PRIMER)")
    # @add_arg_table! settings begin
    #     "--primer1", "-m"
    #         help = "Primers(s) at 5' end of read 1, and their reverse complement appended to 3' end of read 2"
    #         nargs = '+'
    #         metavar = "SEQ"
    #         default = String[]
    #     "--primer2", "-M"
    #         help = "Primers(s) at 5' end of read 2, and their reverse complement appended to 3' end of read 1"
    #         nargs = '+'
    #         metavar = "SEQ"
    #         default = String[]
    #     "--primers", "-P"
    #         help = "Primer table file. Each line is a primer set. Columns are primer1, primer2, primer name and delimited by TAB. Lines starts with `#` are ignored."
    #         metavar = "FILE"
    #         default = ""
    # end

    add_arg_group!(settings, "hard clipping: trim a fixed length")
    @add_arg_table! settings begin
        "--clip-after-r1", "-b"
            help = "hard clip the 3' tails of R1 to contain only INT bases. 0 to disable."
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--clip-after-r2", "-B"
            help = "hard clip the 3' tails of R2 to contain only INT bases. 0 to disable."
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--clip5-r1", "-e"
            help = "remove the first INT bases from 5' front of R1."
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--clip5-r2", "-E"
            help = "remove the first INT bases from 5' front of R2."
            default = 0
            metavar = "INT"
            arg_type = Int64
    end

    add_arg_group!(settings, "quality trimming: trim the tail when the average quality of bases in\na sliding window is low")
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

    add_arg_group!(settings, "N trimming")
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

    add_arg_group!(settings, "length filtration")
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

    add_arg_group!(settings, "complexity filtration")
    @add_arg_table! settings begin
        "--enable-complexity-filtration"
            help = "enable complexity filtration"
            action = :store_true
        "--min-complexity"
            help = "complexity threshold"
            default = 0.3
            metavar = "FLOAT"
            arg_type = Float64
    end

    add_arg_group!(settings, "legacy arguments")
    @add_arg_table! settings begin
        "--procs", "-p"
            help = "ignored (multi-proc is disabled)"
            metavar = "INT"
            default = "1"
            arg_type = String
        "--clip-after", "-C"
            help = "removed (use --clip-after-r1 and --clip-after-r2 instead)"
            default = 0
            metavar = "INT"
            arg_type = Int64
        "--clip5", "-c"
            help = "removed (use --clip5-r1 and --clip5-r2 instead)"
            default = 0
            metavar = "INT"
            arg_type = Int64
    end

    return parse_args(args, settings)
end

function isdna(x::String)
    try
        foreach(x -> convert(DNA, x), collect(x))
        true
    catch
        false
    end
end

function isdna(xs::Vector)
    for x in xs
        if isdna(x)
            continue
        else
            return false
        end
    end
    true
end

function args_range_test(args::Dict{String,Any}; test_only::Bool=false)
    ispass = true

    try
        run(pipeline(`pigz --version`, stderr=devnull))
    catch
        @error "Dependency missing: pigz not found." _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    try
        run(pipeline(`pbzip2 --version`, stderr=devnull))
    catch
        @error "Dependency missing: pbzip2 not found." _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
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

    # if !occursin(r"^[1-9][0-9]*$|^[1-9][0-9]*onlyrun\d+$", args["procs"])
    #     @error "--procs INT must be positive integer" _module=nothing _group=nothing _id=nothing _file=nothing
    #     ispass = false
    # end
    if args["procs"] != "1"
        @warn "--procs (-p) is removed from Atria v3.2.0" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = true
    end

    if args["log2-chunk-size"] < 20 || args["log2-chunk-size"] > 30
        @error "--log2-chunk-size INDEX is not in a suggested range." _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # input
    if length(args["read2"]) != 0 && length(args["read1"]) != length(args["read2"])
        @error "paired-end input: --read1 R1-FASTQ... and --read2 R2-FASTQ... must be paired respectively" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    for file in [args["read1"]; args["read2"]]
        if !isfile(file)
            @error "input fastq file not valid" FILE=file _module=nothing _group=nothing _id=nothing _file=nothing
            ispass = false
        end
    end
    compress = uppercase(args["compress"])
    if !(compress in ["AUTO", "NO", "GZIP", "GZ", "BZIP2", "BZ2"])
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

    # poly X tail trimming
    if args["poly-mismatch-per-16mer"] >= 16 || args["poly-mismatch-per-16mer"] < 0
        @error "--poly-mismatch-per-16mer not in range of 0:15" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if args["poly-length"] <= 0
        @error "--poly-length < 1 is not allowed" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    # adapter triming
    if length(args["read2"]) != 0 && length(args["adapter1"]) != length(args["adapter2"])
        @error "--adapter1 -a SEQ... and --adapter2 -A SEQ...: length of adapter 1 and adapter 2 should be the same. If the same sequence is used in adapter 1 and 2, you need to explicitly specify the sequence to both adapter 1 and 2." ADAPTER1=args["adapter1"] ADAPTER2=args["adapter2"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if !isdna(args["adapter1"])
        @error "--adapter1 -a SEQ... contain(s) invalid bases" SEQ=args["adapter1"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if !isdna(args["adapter2"])
        @error "--adapter2 -A SEQ... contain(s) invalid bases" SEQ=args["adapter2"] _module=nothing _group=nothing _id=nothing _file=nothing
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
    # if args["r1-r2-score-diff"] < 0.0
    #     @error "--r1-r2-score-diff -s FLOAT < 0 is not valid" INT=args["r1-r2-score-diff"] _module=nothing _group=nothing _id=nothing _file=nothing
    #     ispass = false
    # end

    if args["kmer-n-match"] < 0 || args["kmer-n-match"] > 16
        @error "--kmer-n-match -s INT < 0 or > 16 are not valid" INT=args["kmer-n-match"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["trim-score-se"] < 0.0 || args["trim-score-se"] > 16.0
        @error "--trim-score-se FLOAT < 0 or > 16 are not valid" INT=args["trim-score-se"] _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["trim-score-pe"] < 0.0 || args["trim-score-pe"] > 32.0
        @error "--trim-score-pe FLOAT < 0 or > 32 are not valid" INT=args["trim-score-pe"] _module=nothing _group=nothing _id=nothing _file=nothing
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
    if args["clip-after"] != 0
        @error "--clip-after -C is removed. Please use --clip-after-r1 and --clip-after-r2 instead" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end
    if args["clip5"] != 0
        @error "--clip5 -c is removed. Please use --clip5-r1 and --clip5-r2 instead" _module=nothing _group=nothing _id=nothing _file=nothing
        ispass = false
    end

    if 0 < args["clip-after-r1"] < 30
        @warn "--clip-after-r1 -b INT is very small. Did you mis-understand this option? It does not clip INT bases from 3' end, but only keep the first INT bases of each read. For example, if you want to clip 5 bases from 3' end, and your read length is 101, you should use `--clip-after-r1 96` or `-b 96`." INT=args["clip-after-r1"] _module=nothing _group=nothing _id=nothing _file=nothing
    end
    if 0 < args["clip-after-r2"] < 30
        @warn "--clip-after-r2 -B INT is very small. Did you mis-understand this option? It does not clip INT bases from 3' end, but only keep the first INT bases of each read. For example, if you want to clip 5 bases from 3' end, and your read length is 101, you should use `--clip-after-r2 96` or `-B 96`." INT=args["clip-after-r2"] _module=nothing _group=nothing _id=nothing _file=nothing
    end

    if args["clip5-r1"] > 70
        @warn "--clip5-r1 -e INT is very large. Did you mis-understand this option? It does not keep the last INT bases of each read, but clip the first INT bases from 5' end." INT=args["clip5-r1"] _module=nothing _group=nothing _id=nothing _file=nothing
    end
    if args["clip5-r2"] > 70
        @warn "--clip5-r2 -E INT is very large. Did you mis-understand this option? It does not keep the last INT bases of each read, but clip the first INT bases from 5' end." INT=args["clip5-r1"] _module=nothing _group=nothing _id=nothing _file=nothing
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

    if args["min-complexity"] >= 0.8 || args["min-complexity"] <= 0
        @error "--min-complexity is invalid" _module=nothing _group=nothing _id=nothing _file=nothing
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

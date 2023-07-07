@noinline function test_trimmer_and_benchmark()
@testset "trimmer_and_benchmark" begin

    pwd_backup = pwd()

    tmp_path = tempname()
    mkpath(tmp_path)
    cd(tmp_path)

    try
        args = Trimmer.parsing_args(["-r", "peReadSimulated.R1.fastq", "-R", "peReadSimulated.R2.fastq"])
        logjson = Trimmer.OrderedDict()
        logjson["version"] = Trimmer.OrderedDict(
            "julia-version" => VERSION,
            "atria-version" => "9.9.9"
        )
        logjson["arguments"] = sort!(Trimmer.OrderedDict(args))
        fio = open("json","w+")
        Trimmer.JSON.print(fio, sort!(logjson), 4)
        close(fio)

        Trimmer.Distributed.addprocs(1)
        var = ["var"]
        @eval Trimmer.Distributed.@everywhere var = $var
        Trimmer.Distributed.pmap(+,[1,2],[4,5])


        Benchmark.julia_wrapper_simulate(["-o" ,"peReadSimulated", "-x", "2000"])
        Benchmark.julia_wrapper_simulate(["-h"], exit_after_help=false)

        Benchmark.julia_wrapper_randtrim(["peReadSimulated.R1.fastq", "peReadSimulated.R2.fastq"])
        Benchmark.julia_wrapper_randtrim(["-h"])

        if Sys.iswindows()
            julia_wrapper_atria(["-r", "peReadSimulated.R1.randtrim.fastq", "-R", "peReadSimulated.R2.randtrim.fastq", "-c", "8", "--compress", "gz", "-f"])
        else
            run(`pigz --keep peReadSimulated.R1.randtrim.fastq`)
            run(`pigz --keep peReadSimulated.R2.randtrim.fastq`)
            Trimmer.julia_wrapper_atria(["-r", "peReadSimulated.R1.randtrim.fastq.gz", "-R", "peReadSimulated.R2.randtrim.fastq.gz", "-c", "8", "--compress", "gz", "--check-identifier", "-f"])
            julia_wrapper_atria_single_end(["-r", "peReadSimulated.R1.randtrim.fastq.gz", "-c", "8", "--compress", "gz", "-f"])
            julia_wrapper_detect_adapter(["-r", "peReadSimulated.R1.randtrim.fastq.gz", "-c", "8", "--compress", "gz"])

            run(`pbzip2 peReadSimulated.R1.randtrim.fastq`)
            run(`pbzip2 peReadSimulated.R2.randtrim.fastq`)
            julia_wrapper_atria(["-r", "peReadSimulated.R1.randtrim.fastq.gz", "-R", "peReadSimulated.R2.randtrim.fastq.gz", "-c", "8", "--compress", "bz2", "--check-identifier", "-f"])
            julia_wrapper_atria_single_end(["-r", "peReadSimulated.R1.randtrim.fastq.gz", "-c", "8", "--compress", "bz2", "-f"])
        end


        julia_wrapper_atria(["-r", "peReadSimulated.R1.fastq", "-R", "peReadSimulated.R2.fastq", "--polyG", "--enable-complexity-filtration", "-f", "--stats"])
        julia_wrapper_atria_single_end(["-r", "peReadSimulated.R1.fastq",  "--polyG", "--enable-complexity-filtration", "-f", "--stats"])

        # skip finished
        julia_wrapper_atria(["-r", "peReadSimulated.R1.fastq", "-R", "peReadSimulated.R2.fastq", "--polyG", "--enable-complexity-filtration"])
        julia_wrapper_atria_single_end(["-r", "peReadSimulated.R1.fastq",  "--polyG", "--enable-complexity-filtration"])

        # detect adapter
        julia_wrapper_detect_adapter(["-r", "peReadSimulated.R1.fastq",  "--polyG", "--enable-complexity-filtration"])

        julia_wrapper_atria(["-h"], exit_after_help=false)
        julia_wrapper_atria_single_end(["-h"], exit_after_help=false)

        julia_wrapper_readstat(["peReadSimulated.R1.atria.fastq", "peReadSimulated.R2.atria.fastq"])
        julia_wrapper_readstat(["-h"])

        Rscript_check_package = """
        if (is.na(packageDescription("argparse")[1])) writeLines("R package 'argparse' not found. Please run `install.packages('argparse')` in R session.")
        if (is.na(packageDescription("plotly")[1])) writeLines("R package 'plotly' not found. Please run `install.packages('plotly')` in R session.")
        if (is.na(packageDescription("ggsci")[1])) writeLines("R package 'ggsci' not found. Please run `install.packages('ggsci')` in R session.")
        """

        julia_wrapper_rscript(Rscript_check_package, ["-h"])

        atria_markdown_help()

        ARGS_old = deepcopy(ARGS)
        empty!(ARGS)
        push!(ARGS, "prog")
        Atria.julia_main()

        empty!(ARGS)
        append!(ARGS, ARGS_old)

        @info "Precompiling/test passed without errors."

    catch e
        @error "Precompiling/test failed!" exception=e
        cd(pwd_backup)
        rm(tmp_path, recursive=true, force=true)
        rethrow(e)
    finally
        cd(pwd_backup)
        rm(tmp_path, recursive=true, force=true)
    end
end
end
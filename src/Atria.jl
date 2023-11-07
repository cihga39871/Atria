
module Atria

# add ArgParse BioSymbols BioSequences Printf JSON Statistics DelimitedFiles Distributed Logging DataStructures Markdown PrettyTables

# using ArgParse
# using BioSymbols
# using BioSequences
# using Printf
# using JSON
# using Statistics
# using DelimitedFiles
# using Distributed
# using Base.Threads
# using Logging
# using DataStructures
# using Markdown
# using PrettyTables

using Reexport

include(joinpath("BioBits", "BioBits.jl"))
@reexport using .BioBits

include(joinpath("FqRecords", "FqRecords.jl"))
@reexport using .FqRecords

include(joinpath("Trimmer", "Trimmer.jl"))
@reexport using .Trimmer

include(joinpath("Benchmark", "Benchmark.jl"))
@reexport using .Benchmark

include(joinpath("AtriaTest", "AtriaTest.jl"))
@reexport using .AtriaTest


function julia_main()::Cint

    help_programs = """
    Available programs:
        atria       Pair-end trimming software (default)
        simulate    Generate artificial pair-end reads
        randtrim    Randomly trim R1 or R2 at a random position
        readstat    Collect trimming statistics
                        (reads should be generated by `atria simulate`)
        statplot    Plot trimming statistics
                        (`Rscript` in PATH required)
        test        Test Atria program
        p | prog    Show this program list
    """

    if length(ARGS)::Int64 >= 1
        if ARGS[1] in ["prog", "p"]
            println(help_programs)
        elseif ARGS[1] in ("atria", "Atria")
            if "--detect-adapter" in ARGS
                if "-R" in ARGS || "--read2" in ARGS
                    julia_wrapper_detect_adapter_pe(ARGS[2:end])
                else
                    julia_wrapper_detect_adapter_se(ARGS[2:end])
                end
            elseif "-R" in ARGS || "--read2" in ARGS
                # paired-end
                julia_wrapper_atria_pe(ARGS[2:end]::Vector{String})
            else
                julia_wrapper_atria_se(ARGS[2:end]::Vector{String})
            end
        elseif ARGS[1] == "simulate"
            julia_wrapper_simulate(ARGS[2:end]::Vector{String})
        elseif ARGS[1] == "randtrim"
            julia_wrapper_randtrim(ARGS[2:end]::Vector{String})
        elseif ARGS[1] == "readstat"
            julia_wrapper_readstat(ARGS[2:end]::Vector{String})
        elseif ARGS[1] == "statplot"
            julia_wrapper_rscript(statplot_code, ARGS[2:end]::Vector{String})
        elseif ARGS[1] == "test"
            test_atria()
        else
            if "--detect-adapter" in ARGS
                if "-R" in ARGS || "--read2" in ARGS
                    julia_wrapper_detect_adapter_pe(ARGS)
                else
                    julia_wrapper_detect_adapter_se(ARGS)
                end
            elseif "-R" in ARGS || "--read2" in ARGS
                # paired-end
                julia_wrapper_atria_pe(ARGS::Vector{String})
            else
                julia_wrapper_atria_se(ARGS::Vector{String})
            end
        end
    else
        atria_markdown_help()
    end
    return 0
end



end  # module end

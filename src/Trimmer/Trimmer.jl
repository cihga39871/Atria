
module Trimmer

export julia_wrapper_atria_pe,
julia_wrapper_atria_se,
julia_wrapper_detect_adapter_se,
julia_wrapper_detect_adapter_pe,
sub_procs, sub_procs_single_end,
atria_markdown_help,
processing_reads!,
processing_reads_range!,
processing_reads_threads!,
parsing_args,
args_range_test,
get_quality_offset,
get_length_range,
f_procs

using Reexport

@reexport using ArgParse
@reexport using BioSymbols
@reexport using BioSequences
@reexport using Distributed
@reexport using Logging
@reexport using JSON
@reexport using DataStructures
@reexport using Printf
@reexport using Markdown
@reexport using PrettyTables
@reexport using DataFrames
@reexport using CSV
@reexport using Dates
@reexport using Statistics

@reexport using ..BioBits
@reexport using ..FqRecords

using Pkg
const atria_version = @eval($(string("v", Pkg.project().version)))

include("markdown_help.jl")
include("args.jl")
include("thread_trim.jl")
include("wrapper_pe.jl")
include("wrapper_se.jl")
include("detect_adapter.jl")
include("wrapper_detect_adapter_se.jl")
include("wrapper_detect_adapter_pe.jl")
end

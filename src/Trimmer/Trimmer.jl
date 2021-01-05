
module Trimmer

export julia_wrapper_atria,
sub_procs,
atria_markdown_help,
processing_reads!,
processing_reads_range!,
processing_reads_threads!,
parsing_args,
args_range_test,
get_quality_offset,
get_length_range,
f_procs

using ArgParse
using BioSymbols
using BioSequences
using Distributed
using TimerOutputs
using Logging
using JSON
using DataStructures
using Printf
using Markdown

using ..BioBits
using ..BioBits.BioSymbols
using ..BioBits.BioSequences
using ..FqRecords

include("args.jl")
include("thread_trim.jl")
include("wrapper.jl")
include("markdown_help.jl")

end

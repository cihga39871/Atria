
module Benchmark

export julia_wrapper_simulate,
julia_wrapper_randtrim,
julia_wrapper_readstat,
julia_wrapper_rscript,
statplot_code

using ArgParse
using Statistics

using ..BioBits.BioSymbols
using ..BioBits.BioSequences
using ..FqRecords

include("read_simulation.jl")
include("rand_trim.jl")
include("read_stats.jl")
include("external_code.jl")

end


module AtriaTest
export test_atria

using Test

using ..BioBits
using ..BioBits.BioSymbols
using ..BioBits.BioSequences
using ..FqRecords
using ..Trimmer
using ..Benchmark
using ...Atria

#=
using Test
using .Atria
using .Atria.BioBits
using .Atria.BioBits.BioSymbols
using .Atria.BioBits.BioSequences
using .Atria.FqRecords
using .Atria.Trimmer
using .Atria.Benchmark
=#



@noinline function test_atria()
    @testset "Atria" begin
        include(joinpath("BioBits", "runtests.jl"))
        include(joinpath("FqRecords", "runtests.jl"))
        include("trimmer_and_benchmark.jl")
    end
    true
end

end

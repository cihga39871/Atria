
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

include(joinpath("BioBits", "runtests.jl"))
include(joinpath("FqRecords", "runtests.jl"))
include("trimmer_and_benchmark.jl")

@noinline function test_atria()
    @testset "Atria" begin
        @test test_biobits()
        @test test_fqrecords()
        @test test_trimmer_and_benchmark()
    end
    true
end

end

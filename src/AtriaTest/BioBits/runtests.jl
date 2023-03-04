

# include("insert_size_decision.jl")
include("algorithm_basis.jl")
include("biosequences_safety.jl")
include("get_seq.jl")
include("bit_match.jl")

@noinline function test_bio_bits()
    @testset "BioBits" begin
        test_algorithm_basis()
        test_biosequences_safety()
        test_get_seq()
        test_bit_match()
    end
end

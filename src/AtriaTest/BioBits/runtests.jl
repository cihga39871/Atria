
include("algorithm_basis.jl")
include("biosequences_safety.jl")
include("get_seq.jl")
include("bit_match.jl")
# include("insert_size_decision.jl")

@noinline function test_biobits()
    test_algorithm_basis()
    test_biosequences_safety()
    test_get_seq()
    test_bit_match()
    true
end

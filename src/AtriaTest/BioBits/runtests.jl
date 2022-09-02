

# include("insert_size_decision.jl")

@testset "BioBits" begin
    include("algorithm_basis.jl")
    include("biosequences_safety.jl")
    include("get_seq.jl")
    include("bit_match.jl")
end

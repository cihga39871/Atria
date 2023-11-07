include("fq_records.jl")
include("primer_match.jl")

@noinline function test_fq_records()
    @testset "BioBits" begin
        test_fq_records_basis()
        # test_primer_match()
    end
end

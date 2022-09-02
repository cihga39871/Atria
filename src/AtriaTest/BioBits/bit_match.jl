
@testset "bit match" begin
    a = dna"ACCCGGTCAGTACGTCAGTACGCAGTAGTGTA" |> bitsafe!
    b = dna"NNNACCCGGTCAGTACGTCAGTACGCAGTAGTGTA" |> bitsafe!
    c = dna"NNNNACCCGGTCAGTACGTCAGTACGCAGTAGTGTA" |> bitsafe!
    d = dna"GGTCAGTACGTCAGTACGCAGTAGTGTANNNNACCC" |> bitsafe!
    e = dna"GGTCAGTACGTCAGTACGCAGTAGTGTANNNNCCC" |> bitsafe!
    f = dna"GGTCAGTACGTCAGTACGCAGTAGTGTATTTTCCC" |> bitsafe!
    g = dna"GGTCAGTACGTCAGTACGCAGTAGTGTATTTTACCC" |> bitsafe!
    h = dna"GGTCAGTACGTCAGTACGCAGTAGTGTATTTTAC" |> bitsafe!
    i = dna"GGTCAGTACGTCAGTACGCAGTAGTGTATTTAC" |> bitsafe!
    j = dna"GGTCAGTACGTCAGTACGCAGTAGTGTATTTACC" |> bitsafe!

    # to speed up, bitwise_scan does not handle tail match well.
    @test bitwise_scan(a, b, 1, 0) == (4, 16)
    @test bitwise_scan(a, c, 1, 0) == (5, 16)
    @test bitwise_scan(a, d, 1, 0) == (33, 4)
    @test bitwise_scan(a, e, 1, 0) == (31, 4) # actually best is 32,4
    @test bitwise_scan(a, f, 1, 5) == (32, 3)
    @test bitwise_scan(a, g, 1, 5) == (33, 4)
    @test bitwise_scan(a, h, 1, 5) == (33, 2)
    @test bitwise_scan(a, i, 1, 5) == (32, 1) # actually best is 32,2
    @test bitwise_scan(a, j, 1, 5) == (31, 2) # actually best is 32,3
end

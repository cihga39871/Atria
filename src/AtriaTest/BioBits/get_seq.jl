
@testset "get seq" begin
    a = dna"NNNNATCGNNSANNNNNNNNNNNN" |> bitsafe!

    a.data = N2gap.(a.data)
    @test a == dna"----ATCG--SA------------"

    a = dna"ATCGACTGCGTACGTACGTAC" |> bitsafe!
    SeqHeadSet(a)

    b = dna"" |> bitsafe!
    SeqHeadSet(b)

    pa = get_pointer(0x00, a)
    @test unsafe_load(pa) == 0x81

    pa = get_pointer(0x0000, a)
    @test unsafe_load(pa) == 0x4281

    pa = get_pointer(0x00000000, a)
    @test unsafe_load(pa) == 0x48214281

    pa = get_pointer(0x0000000000000000, a)
    @test unsafe_load(pa) == 0x1842184248214281

    @test unsafe_bitseq(pa, 1) == 0x1842184248214281
    @test unsafe_bitseq(pa, 2) == 0x0184218424821428
    @test unsafe_bitseq(pa, 21, 21) == (0x0000000000000002, 1)
end


@noinline function test_algorithm_basis()

    @test UInt === UInt64
    @test Int === Int64
    @test sizeof(UInt64) == 8
    @test sizeof(UInt32) == 4
    @test sizeof(UInt16) == 2
    @test sizeof(UInt8) == 1


    seq = dna"ANATATATATATATGGANNNNATATATNNNGGGG"

    @test typeof(seq) === LongDNA{4}
    @test typeof(seq) === LongSequence{DNAAlphabet{4}}

    @test typeof(seq.data) === Array{UInt64,1}
    @test typeof(seq.len) === UInt

    @test seq.data == UInt64[0x44818181818181f1,
                             0x44fff818181ffff1,
                             0x0000000000000044]

    p_seq = pointer(seq.data)
    @test unsafe_load(p_seq, 2) == 0x44fff818181ffff1
    @test unsafe_load(p_seq + 1) == 0xf144818181818181

end

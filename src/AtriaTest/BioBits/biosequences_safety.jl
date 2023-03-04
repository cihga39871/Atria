@noinline function test_biosequences_safety()
        
    @testset "biosequences safety" begin
        @testset "bitsafe" begin
            s1 = dna""
            s2 = dna"NASTTGGTTATCNNNN"
            s3 = LongDNA{4}([0x4214824181422181, 0x0ff0000084128142, 0x0000000084128142], UInt(40))

            bitsafe!(s1)
            @test length(s1) == 0
            @test length(s1.data) == 1
            # @test s1.data[1] == 0x0

            bitsafe!(s2)
            @test s2.data[1] == 0xffff28188448861f
            @test length(s2.data) == 2

            @test !isbitsafe(s3)
            @test isbitsafe(s1)
            @test isbitsafe(s2)
        end

        @testset "bitsafe resize" begin
            s3 = LongDNA{4}([0x4214824181422181, 0x0ff0000084128142, 0x0000000084128142], UInt(40))
            resize!(s3, 5)
            @test isbitsafe(s3)

            s3 = LongDNA{4}([0x4214824181422181, 0x0ff0000084128142, 0x0000000084128142], UInt(40))
            resize!(s3, 100)
            @test length(s3.data) == 8
            @test isbitsafe(s3)

            s3 = LongDNA{4}([0x4214824181422181, 0x0ff0000084128142, 0x0000000084128142], UInt(40))
            s3 = s3[4:36]
            resize!(s3, 40)
            @test s3.data[1] == 0x1424214824181422
            @test s3.data[2] == 0x1420ff0000084128
            @test s3.data[3] |  0xfffffffffffffff0 == 0xfffffffffffffff8
            @test length(s3.data) == 4
        end

        @testset "bitsafe reverse complement" begin
            s3 = LongDNA{4}([0x4214824181422181, 0x0ff0000084128142, 0x0000000084128142], UInt(40))
            s3_rc = reverse_complement(s3)
            true_s3_rc_data = [0x00000ff042814821,
                            0x8241284242814821,
                            0x0000000081844281,
                            0x0000000000000000]
            @test s3_rc.data[1] == 0x00000ff042814821
            @test s3_rc.data[2] == 0x8241284242814821
            @test s3_rc.data[3] |  0xffffffff00000000 == 0xffffffff81844281

            bitsafe!(s3)
            s3_rc = reverse_complement(s3)
            @test s3_rc.data[1] == 0x00000ff042814821
            @test s3_rc.data[2] == 0x8241284242814821
            @test s3_rc.data[3] |  0xffffffff00000000 == 0xffffffff81844281
            @test length(s3.data) == 4

            s4 = s3[17:25]
            s4_rc = reverse_complement!(s4)
            @test s4_rc.data[1] |  0xfffffff000000000 == 0xfffffff428148210
            @test length(s4.data) == 2
        end
    end
end

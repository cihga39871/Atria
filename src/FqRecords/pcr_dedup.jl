
mutable struct DupCount
    @atomic count::Int
    id::String
end


function write_pcr_dedup_count(out_pcr_dedup_count::AbstractString, dup_dict::Dict{Vector{UInt64}, DupCount})
    dup_count = 0
    open(out_pcr_dedup_count, "w+") do io
        println(io, "count\tid")
        for v in values(dup_dict)
            if v.count > 1
                @inbounds println(io, "$(v.count)\t$(v.id)")
                dup_count += v.count
            end
        end
    end
    dup_count
end

function get_dup_count(dup_dict::Dict)
    dup_count = 0
    for v in values(dup_dict)
        if v.count > 1
            dup_count += v.count
        end
    end
    dup_count
end

function dna_2bit(a::DNA)
    if a == DNA_A
        return 0x00
    elseif a == T
        return 0x11
    elseif a == C
        return 0x01
    elseif a == G
        return 0x10
    else
        return 0x01  # unknown to C
    end
end
function alphabet_dna_2bit()
    ab = Vector{UInt8}(undef, 16)
    fill!(ab, 0x01)  # unknown to C
    ab[reinterpret(UInt8, DNA_A)+1] = 0b00
    ab[reinterpret(UInt8, DNA_T)+1] = 0b11
    ab[reinterpret(UInt8, DNA_G)+1] = 0b10
    ab
end
const ALPHABET_DNA_2BIT = alphabet_dna_2bit()

function alphabet_2dna()
    ab = Vector{UInt8}(undef, 256)
    fill!(ab, 0b01 << 2 | 0b01)  # unknown to CC
    for x in (0b0001, 0b0010, 0b0100, 0b1000)
        x2 = ALPHABET_DNA_2BIT[x+1]
        for y in (0b0001, 0b0010, 0b0100, 0b1000)
            y2 = ALPHABET_DNA_2BIT[y+1]
            double_dna_8bit = x << 4 | y
            ab[double_dna_8bit+1] = x2 << 2 | y2
        end
    end
    ab
end
const ALPHABET_2DNA = alphabet_2dna()


"""
    hash_dna(s1::LongDNA{4})

Hash DNA to `[num_bits_in_it; ::LongDNA{2}.data]`. Ambiguous/Gap DNA converts to C.
"""
function hash_dna(s1::LongDNA{4})
    global ALPHABET_2DNA

    len = length(s1)
    data = zeros(UInt64, 1 + BioSequences.seq_data_len(DNAAlphabet{2}, len))
    
    count_bits = sum(count_ones(x) for x in s1.data)
    @inbounds data[1] = count_bits

    dt_re = reinterpret(reshape, UInt8, data)

    s1_re = reinterpret(reshape, UInt8, s1.data)
    dt_re_offset = 8

    double_dna_size = cld(len, 4)
    @inbounds @simd for i in 1:double_dna_size
        double_dna1 = ALPHABET_2DNA[s1_re[2i-1] + 1]
        double_dna2 = ALPHABET_2DNA[s1_re[2i] + 1]

        dt_re[dt_re_offset + i] = double_dna2 << 4 | double_dna1
    end
    data
end

function hash_dna(s1::LongDNA{4}, s2::LongDNA{4})
    global ALPHABET_2DNA

    len1 = length(s1)
    len2 = length(s2)
    data_len1 = BioSequences.seq_data_len(DNAAlphabet{2}, len1)
    data_len2 = BioSequences.seq_data_len(DNAAlphabet{2}, len2)
    data = zeros(UInt64, 1 + data_len1 + data_len2)
    
    count_bits = sum(count_ones(x) for x in s1.data)
    count_bits += sum(count_ones(x) for x in s2.data)
    @inbounds data[1] = UInt64(count_bits)

    dt_re = reinterpret(reshape, UInt8, data)

    s1_re = reinterpret(reshape, UInt8, s1.data)
    dt_re_offset = 8

    double_dna_size1 = cld(len1, 4)
    @inbounds for i in 1:double_dna_size1
        double_dna1 = ALPHABET_2DNA[s1_re[2i-1] + 1]
        double_dna2 = ALPHABET_2DNA[s1_re[2i] + 1]

        dt_re[dt_re_offset + i] = double_dna2 << 4 | double_dna1
    end

    s2_re = reinterpret(reshape, UInt8, s2.data)
    dt_re_offset += data_len1 * 8

    double_dna_size2 = cld(len2, 4)
    @inbounds for i in 1:double_dna_size2
        double_dna1 = ALPHABET_2DNA[s2_re[2i-1] + 1]
        double_dna2 = ALPHABET_2DNA[s2_re[2i] + 1]

        dt_re[dt_re_offset + i] = double_dna2 << 4 | double_dna1
    end
    data
end

"""
    N2gap(bit::T) where T <: Union{UInt8, UInt16, UInt32, UInt64, UInt128}

Convert N (1111) to gap (0000) in biological `bit`.
"""
function N2gap end

for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
    @eval @inline function N2gap(bit::$T)
        nbase_1 = $(sizeof(T) * 2 - 1)

        N_bit = $(convert(T, 0b1111))
        if bit & N_bit == N_bit
            bit &= ~N_bit
        end
        for i in 1:nbase_1
            N_bit = $(convert(T, 0b1111)) << (4*i)
            if bit & N_bit == N_bit
                bit &= ~N_bit
            end
        end
        bit
    end
end


struct SeqHead{T}
    a::T
    b::T
    function SeqHead{T}(a::T, b::T) where T <: Union{UInt8, UInt16, UInt32, UInt64}
        new(a, b)
    end
end

@inline function SeqHead(a::T, b::T) where T <: Union{UInt8, UInt16, UInt32, UInt64}
    SeqHead{T}(a,b)
end

"""
    SeqHead(::T, seq::LongDNASeq) where T <: Union{UInt8, UInt16, UInt32, UInt64}

# Fields

 - `a::T`: the bits of sequence from index 1.

 - `b::T`: the bits of sequence from index 2.

# Argument

 - `seq::LongDNASeq`: the seq has to be `bitsafe!`.
"""
function SeqHead(::T, seq::LongDNASeq) where T <: Union{UInt8, UInt16, UInt32}
    bit = seq.data[1]
    a = unsafe_trunc(T, bit)
    b = unsafe_trunc(T, bit >> 4)
    SeqHead{T}(a, b)
end
function SeqHead(::UInt64, seq::LongDNASeq)
    p = pointer(seq.data)
    a = unsafe_load(p)
    if length(seq) > 0
        c = unsafe_load(p+1)
        b = (a >> 4) | (c << 4)
    else
        b = a >> 4
    end
    SeqHead{UInt64}(a, b)
end
for T in (UInt8, UInt16, UInt32, UInt64)
    @eval @inline SeqHead{$T}(seq::LongDNASeq) = SeqHead($(typemin(T)), seq)
end


struct SeqHeadSet
    s64::SeqHead{UInt64}
    s32::SeqHead{UInt32}
    s16::SeqHead{UInt16}
    s8::SeqHead{UInt8}
    function SeqHeadSet(seq::LongDNASeq)
        s64 = SeqHead{UInt64}(seq)
        s32 = SeqHead{UInt32}(seq)
        s16 = SeqHead{UInt16}(seq)
        s8 = SeqHead{UInt8}(seq)
        new(s64, s32, s16, s8)
    end
end

function BioSequences.LongDNASeq(s::SeqHeadSet)
    LongDNASeq([s.s64.a], 1:16, false)
end

"""
    TruncSeq(::T, seq::LongDNASeq) where T <: Union{UInt8, UInt16, UInt32, UInt64}

# Fields

 - `a::T`: the bits of sequence from index 1.

 - `b::T`: the bits of sequence from index 2.

 - `a1::T`: the bits of `seq[1]`.

# Argument

 - `seq::LongDNASeq`: the seq has to be `bitsafe!`.
"""
struct TruncSeq{T}
    a::T
    b::T
    a1::T
    function TruncSeq{T}(a::T, b::T, a1::T) where T <: Union{UInt8, UInt16, UInt32, UInt64}
        new(a, b, a1)
    end
end

for T in (UInt8, UInt16, UInt32, UInt64)
    @eval @inline TruncSeq(a::$T, b::$T, a1::$T) = TruncSeq{$T}(a,b,a1)
end

for T in (UInt8, UInt16, UInt32)
@eval @inline function TruncSeq(::$T, seq::LongDNASeq)
    bit = seq.data[1] #|> N2gap
    a = unsafe_trunc($T, bit)
    b = unsafe_trunc($T, bit >> 4)
    a1 = $(T(0b1111)) & a
    TruncSeq{T}(a, b, a1)
end
end
function TruncSeq(::UInt64, seq::LongDNASeq)
    p = pointer(seq.data)
    a = unsafe_load(p) #|> N2gap
    c = unsafe_load(p+1) #|> N2gap
    if length(seq) > 0
        c = unsafe_load(p+1)
        b = (a >> 4) | (c << 4)
    else
        b = a >> 4
    end
    a1 = 0x000000000000000f & a
    TruncSeq{UInt64}(a, b, a1)
end


for T in (UInt8, UInt16, UInt32, UInt64)
    @eval @inline TruncSeq{$T}(seq::LongDNASeq) = TruncSeq($(typemin(T)), seq)
end

"""
    get_pointer(::T, seq::LongDNASeq) where T <: {UInt8, UInt16, UInt32, UInt64}
"""
@inline get_pointer(::UInt64, seq::LongDNASeq) = pointer(seq.data)
for T in (UInt8, UInt16, UInt32)
    @eval @inline get_pointer(::$T, seq::LongDNASeq) =
        Core.bitcast($(Ptr{T}), pointer(seq.data))
end


"""
    get_unsafe_index_of_last_bitseq(::T, seq::LongDNASeq)
    get_unsafe_index_of_last_bitseq(::T, seq.part::UnitRange{Int64})
    get_unsafe_index_of_last_bitseq(::T, seq.part.stop::Int64)

 - `::T` is one of UInt8, UInt16, UInt32, UInt64.

Get the index of the last full-long bitseq. It is unsafe because the returned index can be negative.
"""
function get_unsafe_index_of_last_bitseq end

for T in (UInt8, UInt16, UInt32, UInt64)
    @eval @inline get_unsafe_index_of_last_bitseq(::$T, seq::LongDNASeq) =
        seq.part.stop - $(sizeof(T) * 2 - 2)
    @eval @inline get_unsafe_index_of_last_bitseq(::$T, seq_part::UnitRange{Int64}) =
        seq_part.stop - $(sizeof(T) * 2 - 2)
    @eval @inline get_unsafe_index_of_last_bitseq(::$T, seq_part_stop::Int64) =
        seq_part_stop - $(sizeof(T) * 2 - 2)
end

"""
    unsafe_bitseq(seq_data_ptr::Ptr{T}, idx::Int) => bitseq
    unsafe_bitseq(seq_data_ptr::Ptr{T}, idx::Int, max_idx::Int) => bitseq, num_base_extracted

 - `seq_data_ptr::Ptr{T}`: the pointer to `(seq::LongDNASeq).data`. `Ptr{T}` can be converted to `Ptr` of `UInt8`, `UInt16`, `UInt32`, or `UInt64`.

 - `idx`: nucleotide index of `(seq::LongDNASeq).data`. It is not adjusted from `(seq::LongDNASeq).part.start`.

 - `max_idx`: should be equal to `(seq::LongDNASeq).part.stop`. Change bits after it to 0. It does not mask bits if `max_idx` < `(seq::LongDNASeq).part.stop`, but affects num_base_extracted.

# Caution

When `idx` is even, the bitseq will always start from 0b0000, because it simply shift 4 bits from `idx - 1`.
"""
function unsafe_bitseq end

for T in (UInt8, UInt16, UInt32, UInt64)
@eval @inline function unsafe_bitseq(seq_data_ptr::Ptr{$T}, idx::Int)
    idx_c = idx - 1
    bitseq = unsafe_load(seq_data_ptr + idx_c ÷ 2)
    access_by_shift = idx_c % 2 == 1
    if access_by_shift
        # cannot accee to this index directly
        # INFO:
        bitseq >>= 0x04
    end
    return bitseq
end
end


for T in (UInt8, UInt16, UInt32, UInt64)
@eval @inline function unsafe_bitseq(seq_data_ptr::Ptr{$T}, idx::Int, max_idx::Int)
    idx_c = idx - 1
    bitseq = unsafe_load(seq_data_ptr + idx_c ÷ 2)
    access_by_shift = idx_c % 2 == 1
    if access_by_shift
        # cannot accee to this index directly
        bitseq >>= 0x04
    end

    nbase = $(sizeof(T) * 2)
    idx_stop = idx_c + nbase
    nbase_overflow = idx_stop - max_idx
    if nbase_overflow > 0
        # mask bases after idx_stop. but it is assumed masked by bitsafe!
        # bitseq &= ($(typemax(T)) >> (nbase_overflow * 4))
        num_base_extracted = nbase - nbase_overflow
    else
        num_base_extracted = access_by_shift ? nbase - 1 : nbase
    end
    bitseq, num_base_extracted
end
end


function bin(x)
    replace(bitstring(x), r"(....)" => s"\1 ")
end

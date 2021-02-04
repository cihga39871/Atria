


"""
    bitwise_scan(seq_head_set::SeqHeadSet, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

    bitwise_scan(a::LongDNASeq, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

 - `a::LongDNASeq`: it will be convert to `SeqHeadSet`, and then match to `b`.

 - `b::LongDNASeq`: the sequence to match from `from` to `until`.

 - `until::Int64`: scan until this index of `b`. It does not guarantee the bases after it are encoded as 0000. To guarantee this, subset `b` and then call `bitsafe!(b_subset)`.
"""
@inline function bitwise_scan(seq_head_set::SeqHeadSet, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    bp = get_pointer(0x0000000000000000, b)
    stop = b.part.stop
    stop_fullseq = get_unsafe_index_of_last_bitseq(0x0000000000000000, stop)
    ptr_from = (from - 1) ÷ 2
    ptr_stop = (stop - 1) ÷ 2
    ptr_until = (until - 1) ÷ 2
    ptr_stop_fullseq = (stop_fullseq - 1) ÷ 2

    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s64, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, 0, 0)

    if best_ncompatible >= 8 || ptr_until <= ptr_stop_fullseq
        return (best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    ptr_stop_fullseq += 4
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s32, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, best_idx, best_ncompatible)

    if best_ncompatible >= 4 || ptr_until <= ptr_stop_fullseq
        return (best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    ptr_stop_fullseq += 2
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s16, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, best_idx, best_ncompatible)

    if best_ncompatible >= 2 || ptr_until <= ptr_stop_fullseq
        return (best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s8, bp, ptr_from, ptr_stop, allowed_mismatch, best_idx, best_ncompatible)

    # happens when length seq is odd, and some match found at the end of base.
    if best_idx > stop
        return 0, 0
    end

    return (best_idx, best_ncompatible)
end

@inline bitwise_scan(a::LongDNASeq, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807) = bitwise_scan(SeqHeadSet(a), b, from, allowed_mismatch; until = until)

# internal for bitwise_scan!
for T in (UInt8, UInt16, UInt32, UInt64)
@eval @inline function _bitwise_scan_fullseq(seq_head::SeqHead{$T}, bp::Ptr, ptr_from::Int64, ptr_stop_fullseq::Int64, allowed_mismatch::Int64, best_idx::Int64, best_ncompatible::Int64)
    bp = Core.bitcast($(Ptr{T}), bp)
    min_compatible_fullseq = $(sizeof(T) * 2) - allowed_mismatch

    @inbounds while ptr_from <= ptr_stop_fullseq
        # b_bit = unsafe_load(bp + ptr_from)  # even slower defining it?
        ncompatible = count_ones(unsafe_load(bp + ptr_from) & seq_head.b)
        if ncompatible >= min_compatible_fullseq && ncompatible > best_ncompatible
            best_idx = ptr_from * 2
            best_ncompatible = ncompatible
        end

        ncompatible = count_ones(unsafe_load(bp + ptr_from) & seq_head.a)
        if ncompatible >= min_compatible_fullseq && ncompatible > best_ncompatible
            best_idx = ptr_from * 2 + 1
            best_ncompatible = ncompatible
        end

        ptr_from += 1
    end

    best_idx, best_ncompatible
end
end

"""
    bitwise_scan_rc!(rc_dest::LongDNASeq, a, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

 - `a::Union{LongDNASeq,SeqHeadSet}`: it will be convert to `SeqHeadSet`, and then match to the reverse-complementary `b`.

 - `from::Int64`: after computing the reverse complement of `b` (`b_rc`), the index of `b_rc` to search from.

 - `until::Int64`: scan until this index of `b_rc`.

Note: the reverse complement process is `bitsafe!`.
"""
@inline function bitwise_scan_rc!(rc_dest::LongDNASeq, a, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    copy!(rc_dest, b)
    b_rc = reverse_complement!(rc_dest)
    best_idx_rc, best_ncompatible = bitwise_scan(a, b_rc, from, allowed_mismatch; until = until)
    best_idx_rc == 0 && return 0,0
    best_idx = length(b) - best_idx_rc + 1
    best_idx, best_ncompatible
end

"""
    bitwise_scan_rc(a, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

- `a::Union{LongDNASeq,SeqHeadSet}`: it will be convert to `SeqHeadSet`, and then match to the reverse-complementary `b`.

- `from::Int64`: after computing the reverse complement of `b` (`b_rc`), the index of `b_rc` to search from.

- `until::Int64`: scan until this index of `b_rc`.

Note: the reverse complement process is `bitsafe!`.
"""
@inline function bitwise_scan_rc(a, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    rc_dest = LongDNASeq()
    bitwise_scan_rc!(rc_dest, a, b, from, allowed_mismatch; until = until)
end


"""
    bitwise_scan(trunc_seq::TruncSeq, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

    bitwise_scan(a::LongDNASeq, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

> Caution: `TruncSeq` version of `bitwise_scan` is still under development and it is not mature.

The fastest scaning method. It is robust when matching tail of b.

 - `a::LongDNASeq`: it will be convert to `TruncSeq{UInt64}`, and then match to `b`.

 - `b::LongDNASeq`: the sequence to match from `from` to `until`.

 - `until::Int64`: scan until this index of `b`. It does not guarantee the bases after it are encoded as 0000. To guarantee this, subset `b` and then call `bitsafe!(b_subset)`.
"""
function bitwise_scan end

for T in (UInt8, UInt16, UInt32, UInt64)
@eval @inline function bitwise_scan(trunc_seq::TruncSeq{$T}, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    bp = get_pointer(0x0000000000000000, b)
    stop = b.part.stop
    ptr_from = (from - 1) ÷ 2
    ptr_stop = (stop - 1) ÷ 2

    best_idx = -1
    best_ncompatible = 0

    if ptr_from == 0
        nbase = min(stop - ptr_from * 2, $(sizeof(T) * 2))
        this_allowed_mismatch = cld(allowed_mismatch * nbase, $(sizeof(T) * 2))
        min_compatible = max(nbase - this_allowed_mismatch, 1)

        ncompatible = count_ones(unsafe_load(bp + ptr_from) & trunc_seq.a)
        if ncompatible >= min_compatible && ncompatible > best_ncompatible
            best_idx = ptr_from * 2 + 1
            best_ncompatible = ncompatible
        end

        ptr_from += 1
    end

    # ptr_from > 0
    @inbounds while ptr_from <= ptr_stop
        # compare to trunc_seq.b
        nbase = min(stop - ptr_from * 2, $(sizeof(T) * 2))
        this_allowed_mismatch = cld(allowed_mismatch * nbase, $(sizeof(T) * 2))
        min_compatible = max(nbase - this_allowed_mismatch, 1)

        if ptr_from != 0
            ncompatible = count_ones(unsafe_load(bp + ptr_from) & trunc_seq.b)
            # check compatible of ignored base (idx = 2ptr_from)
            ncompatible += (unsafe_load(bp + ptr_from - 1) >> 4 ) & trunc_seq.a1 == trunc_seq.a1
            if ncompatible >= min_compatible && ncompatible > best_ncompatible
                best_idx = ptr_from * 2
                best_ncompatible = ncompatible
            end
        end

        # compare to trunc_seq.a
        ncompatible = count_ones(unsafe_load(bp + ptr_from) & trunc_seq.a)
        if ncompatible >= min_compatible && ncompatible > best_ncompatible
            best_idx = ptr_from * 2 + 1
            best_ncompatible = ncompatible
        end

        ptr_from += 1
    end

    return (best_idx, best_ncompatible)
end
end

# @inline bitwise_scan(a::LongDNASeq, b::LongDNASeq, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807) = bitwise_scan(TruncSeq{UInt64}(a), b, from, allowed_mismatch; until = until)

# trunc_seq = BioBits.TruncSeq(typemin(UInt64), dna"GGGGGGGGGGGGGGGGGGGGGGGGG")

# @benchmark bitwise_scan(seq_head_set, b, 1, 2)
# @benchmark bitwise_scan(trunc_seq, b, 1, 2)

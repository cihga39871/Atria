
mutable struct MatchRes
    idx::Int64
    ncompatible::Int64
    prob::Float64
    score::Float64
end
function MatchRes(idx::Int64, ncompatible::Int64)
    MatchRes(idx, ncompatible, NaN, NaN)
end
function Base.:(==)(a::MatchRes, b::MatchRes)
    (a.idx == b.idx) &&
    (a.ncompatible == b.ncompatible) &&
    (a.prob == b.prob || (isnan(a.prob) && isnan(b.prob))) &&
    (a.score == b.score || (isnan(a.score) && isnan(b.score)))
end

function Base.isless(a::MatchRes, b::MatchRes)
    a.ncompatible < b.ncompatible
end

"""
    bitwise_scan(seq_head_set::SeqHeadSet, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

    bitwise_scan(a::LongDNA{4}, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

 - `a::LongDNA{4}`: it will be convert to `SeqHeadSet`, and then match to `b`.

 - `b::LongDNA{4}`: the sequence to match from `from` to `until`.

 - `until::Int64`: scan until this index of `b`. It does not guarantee the bases after it are encoded as 0000. To guarantee this, subset `b` and then call `bitsafe!(b_subset)`.
"""
@inline function bitwise_scan(seq_head_set::SeqHeadSet, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    bp = get_pointer(0x0000000000000000, b)
    stop = b.len % Int64
    stop_fullseq = get_unsafe_index_of_last_bitseq(0x0000000000000000, stop)
    ptr_from = (from - 1) ÷ 2
    ptr_stop = (stop - 1) ÷ 2
    ptr_until = (until - 1) ÷ 2
    ptr_stop_fullseq = (stop_fullseq - 1) ÷ 2

    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s64, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, 0, 0)

    if best_ncompatible >= 8 || ptr_until <= ptr_stop_fullseq
        return MatchRes(best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    ptr_stop_fullseq += 4
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s32, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, best_idx, best_ncompatible)

    if best_ncompatible >= 4 || ptr_until <= ptr_stop_fullseq
        return MatchRes(best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    ptr_stop_fullseq += 2
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s16, bp, ptr_from, min(ptr_stop_fullseq, ptr_until), allowed_mismatch, best_idx, best_ncompatible)

    if best_ncompatible >= 2 || ptr_until <= ptr_stop_fullseq
        return MatchRes(best_idx, best_ncompatible)
    end

    ptr_from = ptr_stop_fullseq + 1
    allowed_mismatch = cld(allowed_mismatch, 2)
    best_idx, best_ncompatible = _bitwise_scan_fullseq(seq_head_set.s8, bp, ptr_from, ptr_stop, allowed_mismatch, best_idx, best_ncompatible)

    # happens when length seq is odd, and some match found at the end of base.
    if best_idx > stop
        return MatchRes(0, 0)
    end

    return MatchRes(best_idx, best_ncompatible)
end

@inline bitwise_scan(a::LongDNA{4}, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807) = bitwise_scan(SeqHeadSet(a), b, from, allowed_mismatch; until = until)

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
    bitwise_scan_rc!(rc_dest::LongDNA{4}, a, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

 - `a::Union{LongDNA{4},SeqHeadSet}`: it will be convert to `SeqHeadSet`, and then match to the reverse-complementary `b`.

 - `from::Int64`: after computing the reverse complement of `b` (`b_rc`), the index of `b_rc` to search from.

 - `until::Int64`: scan until this index of `b_rc`.

 - `init_rc_dest::Bool`: whether `rc_dest` need to be computed to be the reverse complement of `b`. If false, `rc_dest` has been initialized.

Note: the reverse complement process is `bitsafe!`.
"""
@inline function bitwise_scan_rc!(rc_dest::LongDNA{4}, a, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807, init_rc_dest::Bool = true)
    if init_rc_dest
        copy!(rc_dest, b)
        reverse_complement!(rc_dest)
    end
    match_res = bitwise_scan(a, rc_dest, from, allowed_mismatch; until = until)

    if match_res.idx == 0
        match_res.ncompatible = 0
        return match_res
    end
    match_res.idx = length(b) - match_res.idx + 1  # idx compute from rc_dest to b
    match_res
end

"""
    bitwise_scan_rc(a, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

The fastest scaning method. It is robust when matching tail of b.

- `a::Union{LongDNA{4},SeqHeadSet}`: it will be convert to `SeqHeadSet`, and then match to the reverse-complementary `b`.

- `from::Int64`: after computing the reverse complement of `b` (`b_rc`), the index of `b_rc` to search from.

- `until::Int64`: scan until this index of `b_rc`.

Note: the reverse complement process is `bitsafe!`.
"""
@inline function bitwise_scan_rc(a, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    rc_dest = LongDNA{4}()
    bitwise_scan_rc!(rc_dest, a, b, from, allowed_mismatch; until = until)
end


"""
    bitwise_scan(trunc_seq::TruncSeq, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

    bitwise_scan(a::LongDNA{4}, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = typemax(Int64))

> Caution: `TruncSeq` version of `bitwise_scan` is still under development and it is not mature.

The fastest scaning method. It is robust when matching tail of b.

 - `a::LongDNA{4}`: it will be convert to `TruncSeq{UInt64}`, and then match to `b`.

 - `b::LongDNA{4}`: the sequence to match from `from` to `until`.

 - `until::Int64`: scan until this index of `b`. It does not guarantee the bases after it are encoded as 0000. To guarantee this, subset `b` and then call `bitsafe!(b_subset)`.
"""
function bitwise_scan end

for T in (UInt8, UInt16, UInt32, UInt64)
@eval @inline function bitwise_scan(trunc_seq::TruncSeq{$T}, b::LongDNA{4}, from::Int64, allowed_mismatch::Int64; until::Int64 = 9223372036854775807)
    bp = get_pointer(0x0000000000000000, b)
    stop = b.len % Int64
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

    return MatchRes(best_idx, best_ncompatible)
end
end

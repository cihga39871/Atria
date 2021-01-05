
@inline function isinreadlength!(r::FqRecord, length_range::UnitRange{Int64})::Bool
    (length(r.seq)::Int64 in length_range)::Bool
end
@inline function isinreadlength!(r1::FqRecord, r2::FqRecord, length_range::UnitRange{Int64})::Bool
    res1 = (length(r1.seq)::Int64 in length_range)::Bool
    res2 = (length(r2.seq)::Int64 in length_range)::Bool
    (res1 && res2)
end

@inline function count_N(r::FqRecord)::Float64
    # A/T/G/C: 0001, 0010, 0100, 1000: count_ones == 1
    # N      : 1111                  : count_ones == 4
    n_1s = 0
    for b in r.seq.data
        n_1s += count_ones(b::UInt64)::Int64
    end
    @fastmath((n_1s - length(r.seq)::Int64)::Int64 / 3.0)::Float64
end

@inline function isnotmuchN!(r1::FqRecord, r2::FqRecord, max_N::Int64)::Bool
    c1 = count_N(r1)::Float64
    c2 = count_N(r2)::Float64
    res1 = (c1 <= max_N::Int64)
    res2 = (c2 <= max_N::Int64)
    res1 && res2
end

@inline function front_trim!(r::FqRecord, ntrim::Int64)::Nothing
    if ntrim <= 0
    elseif ntrim < length(r.seq)
        delete_range = 1:ntrim
        deleteat!(r.seq, delete_range)
        deleteat!(r.qual, delete_range)
        deleteat!(r.prob, delete_range)
    else  # ntrim >= length(r.seq)
        resize!(r.seq, 0)
        resize!(r.qual, 0)
        resize!(r.prob, 0)
    end
    return
end

# @inline function tail_trim!(r::FqRecord, m::AlignMatch)::Nothing
#     resize!(r.seq::LongDNASeq, m.insert_size::Int64)
#     resize!(r.qual, m.insert_size)
#     resize!(r.prob, m.insert_size)
#     return
# end

@inline function tail_trim!(r::FqRecord, nremain::Int64)::Nothing
    if nremain < length(r.seq::LongDNASeq)
        resize!(r.seq::LongDNASeq, nremain::Int64)
        resize!(r.qual, nremain)
        resize!(r.prob, nremain)
    end
    return
end

@inline function tail_N_trim!(r::FqRecord)::Nothing
    nbase = length(r.seq::LongDNASeq)::Int64
    # trim end
    n = nbase::Int64
    @inbounds while n::Int64 >= 1
        (r.seq::LongDNASeq)[n]::DNA == DNA_N ? n -= 1 : break
    end
    if n::Int64 != nbase::Int64
        resize!(r.seq::LongDNASeq, n::Int64)
        resize!(r.qual, n)
        resize!(r.prob, n)
    end
    return
end


"""
    qualitymatch(r::FqRecord, q0::UInt8, qn::UInt64, n::Int64)::Int64

# ARGUMENTS
1. `r::FqRecord` is FastQ record.
2. `q0::UInt8` is the adjusted quality score. (Eg: +33 if Illumina 1.9+ version).
3. `qn::UInt64` is the adjusted quality score * n. (Eg: +33 if Illumina 1.9+ version).
4. `n::Int64` is the length of sliding window to iterate the reads.

Return the length `n` of reads to keep. `-1` means no need for quality trimming.
"""
@inline function qualitymatch(r::FqRecord, q0::UInt8, qn::UInt64, n::Int64)::Int64
    quals = r.qual
    nqual = length(quals)
    N = n - 1
    nbase = nqual - N
    i = 1

    ### check any qual less than q0
    @inbounds while i <= nqual
        if quals[i] < q0
            break  # start matching sliding window
        end
        i += 1
    end

    (i > nqual) && return -1  # no qual less than q0: not trim
    (i > nbase) && @goto tail_qual_match  # i in the last n bases, go to tail_qual_match

    ### check sliding window
    qual_sum = UInt64(@inbounds quals[i])::UInt64
    start = i + 1
    stop = i + N
    @inbounds for m in start:stop
        qual_sum += quals[m]
    end

    (qual_sum < qn) && @goto tail_qual_match  # ith failed quality match

    i += 1
    while i <= nbase
        @inbounds qual_sum += quals[i+N]
        @inbounds qual_sum -= quals[i-1]
        (qual_sum < qn) && @goto tail_qual_match  # ith failed quality match
        i += 1
    end

    @label tail_qual_match
    while i <= nqual
        (quals[i] < q0) && return i-1  # ith failed quality match
        i += 1
    end

    return -1  # no trim
end

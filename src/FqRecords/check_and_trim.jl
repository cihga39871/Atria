
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

@inline function isnotmuchN!(r::FqRecord, max_N::Int64)::Bool
    c1 = count_N(r)::Float64
    c1 <= max_N
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

"""
    seq_complexity(r::FqRecord)
    seq_complexity(seq::LongDNASeq)

The complexity is defined as the percentage of bases that are different from their next bases (base[i] != base[i+1]). However, here we use an approximation algorithm.

The performance of the algorithm:
```
# Test Sequence            True  Computed Complexity
NNNNNNNNNNNNNNNNNNNNNNNN: (0.0  -2.8260869565217392)
------------------------: (0.0  1.0)
AAAAAAAAAAAAAAAAAAAAAAAA: (0.0  0.04347826086956519)
ATATATATATATATATATATATAT: (1.0  1.0)
ATTATTATTATTATTATTATTATT: (0.65 0.6521739130434783)
ATATATATGGGGGGGG        : (0.5  0.5333333333333333)
NANANANANANANANA        : (NaN  0.0)
```
"""
@inline function seq_complexity(seq::LongDNASeq)
    nbase = seq.part.stop  # cannot use length(r.seq) because seq may start from mid, which is not compatible with the algorithm
    seq_data = seq.data
    n_valid_seq_data = length(seq_data) - 1  # -1 because of bitsafe
    n_ones = 0
    for i in 1:n_valid_seq_data
        b = seq_data[i]
        n_ones += count_ones(b & (b << 4))
        # Test Sequence                            True   Computed Complexity (1 - x/15)
        # NNNNNNNNNNNNNNNN: 60 ones, 4 zeros      (0.0    -3.0)
        # ----------------: 0 ones, 64 zeros      (0.0    1.0)
        # AAAAAAAAAAAAAAAA: 15 ones, 49 zeros     (0.0    0.0)
        # ATATATATATATATAT: 0 ones, 64 zeros      (1.0  1.0)
        # ATTATTATTATTATTA: 5 ones, 59 zeros      (0.65   0.6666666666666667)
        # ATATATATGGGGGGGG: 7 ones, 57 zeros      (0.50   0.5333333333333333)
        # NANANANANANANANA: 15 ones, 49 zeros     (NN    0.0)
    end
    n_compensate = nbase % 16
    if n_compensate == 0
        complexity = @fastmath(1 - n_ones / (15*n_valid_seq_data))
    else
        complexity = @fastmath(1 - n_ones / (15*(n_valid_seq_data - 1) + n_compensate))
    end
end

@inline seq_complexity(r::FqRecord) = seq_complexity(r.seq)


@inline function polyX_tail_scan(a::DNA, b::LongDNASeq, allowed_mismatch_per_16mer::Int64; until::Int64 = 1)
    best_idx = 0
    n = length(b)
    n_mismatch = 0
    allowed_mismatch = allowed_mismatch_per_16mer
    n_polyX_length = 0
    while n >= until
        if b[n] === a
            best_idx = n
        else
            n_mismatch += 1
            n_mismatch <= allowed_mismatch && break
        end
        n -= 1
        n_polyX_length += 1
        if n_polyX_length % 16 == 0
            allowed_mismatch += allowed_mismatch_per_16mer
        end
    end
    best_idx, n_polyX_length
end

@inline polyX_tail_scan(a::DNA, b::FqRecord, allowed_mismatch_per_16mer::Int64; until::Int64 = 1) = polyX_tail_scan(a, b.seq, allowed_mismatch_per_16mer; until = until)
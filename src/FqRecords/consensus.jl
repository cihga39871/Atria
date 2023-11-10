
"""
    pe_consensus!(r1::FqRecord, r2::FqRecord, r2_seq_rc::LongDNA{4}, insert_size::Int64; min_ratio_mismatch::Float64 = 0.2, prob_diff::Float64 = 0.0)

Paired-end consensus calling for read pairs with adapters trimmed. Return `is_consensused::Bool`.
"""
function pe_consensus!(r1::FqRecord, r2::FqRecord, r2_seq_rc::LongDNA{4}, insert_size::Int64; min_ratio_mismatch::Float64 = 0.2, prob_diff::Float64 = 0.0)

    r1_seq = r1.seq
    r2_seq = r2.seq
    r1_length = length(r1_seq)
    r2_length = length(r2_seq_rc)

    # get the overlapped region
    if r2_length < insert_size
        r1_i = insert_size - r2_length + 1
        # check 8 bit alignment
        if r1_i % 2 == 1
            r2_seq_rc.len = r2_length
        else
            r1_i += 1
            deleteat!(r2_seq_rc, 1)
            BioBits.unsafe_extra_bits_to_zeros!(r2_seq_rc)  # deleteat is not bitsafe, so have to use it.
            # r2_seq_rc.part = 2:r2_length
        end
    else
        deleteat!(r2_seq_rc, 1:r2_length-insert_size)
        BioBits.unsafe_extra_bits_to_zeros!(r2_seq_rc)  # deleteat is not bitsafe, so have to use it.
        # r2_seq_rc.part = (r2_length-insert_size+1):r2_length
        r1_i = 1
    end
    length_overlap = min(length(r2_seq_rc), r1_length - r1_i + 1)
    length_overlap <= 0 && return false, 0.0

    # align r2_seq_rc.data and r1_seq.data
    if length(r2_seq_rc) != length_overlap  # when r1 length < insert size
        resize!(r2_seq_rc, length_overlap)
    end

    # Ptr{UInt32}: scan 8 bases each time
    p1 = get_pointer(0x0000000000000000, r1_seq)
    p2_rc = get_pointer(0x0000000000000000, r2_seq_rc)

    p1_offset = r1_i ÷ 2  # r1_i cannot be even, so r1_i -1 not necessary
    p2_rc_offset = 0
    offset_max = cld(length_overlap, 2)
    # the start of ncompatible, minus the score of extra tail match
    num_ones = length_overlap - cld(offset_max,8)*16

    max_num_ones = floor(Int, (1+min_ratio_mismatch) * length_overlap)

    while p2_rc_offset <= offset_max

        if num_ones > max_num_ones
            return false
        end

        # global ncompatible
        # global p1_offset
        # global p2_rc_offset
        # global num_ones
        bit1 = unsafe_load(p1 + p1_offset)
        num_ambiguous_bits = count_ones(bit1) - 16
        bit2 = N2gap(unsafe_load(p2_rc + p2_rc_offset))
        num_ones += count_ones(bit1|bit2) - num_ambiguous_bits

        p1_offset += 8
        p2_rc_offset += 8
    end

    # ratio_mismatch = (num_ones - length_overlap) / length_overlap
    # ratio_mismatch > min_ratio_mismatch && return false
    
    # equals to num_ones > (1-min_ratio_mismatch) * length_overlap && return false
    # see max_num_ones


    # start comsensus calling
    r1_end = min(r1_length, insert_size)
    r2_i = insert_size - r1_i + 1

    r1_qual = r1.qual
    r2_qual = r2.qual

    r1_prob = r1.prob
    r2_prob = r2.prob

    @inbounds while r1_i <= r1_end
        a = r1_seq[r1_i]
        b = r2_seq[r2_i]
        # if !((a | b) in (DNA_W, DNA_S)) # not complement
        if !iscomplement(a, b) # not complement
            a_prob = r1_prob[r1_i]
            b_prob = r2_prob[r2_i]
            if a_prob - b_prob > prob_diff
                # modify b to a
                r2_seq[r2_i] = complement(a)
                r2_qual[r2_i] = r1_qual[r1_i]
                r2_prob[r2_i] = a_prob
            elseif b_prob - a_prob > prob_diff
                r1_seq[r1_i] = complement(b)
                r1_qual[r1_i] = r2_qual[r2_i]
                r1_prob[r1_i] = b_prob
            end
        end
        r1_i += 1
        r2_i -= 1
    end
    return true
end

"""
    pe_consensus!(r1::FqRecord, r2::FqRecord, r1_seq_rc::LongDNA{4}, r2_seq_rc::LongDNA{4}; kmer_tolerance::Int64 = 2, overlap_score::Float64 = 0.0, min_ratio_mismatch::Float64 = 0.2, prob_diff::Float64 = 0.0)

Paired-end consensus calling for read pairs without adapters. Check whether the read pair has an overlapped region first. Return `is_consensused::Bool`.
"""
function pe_consensus!(r1::FqRecord, r2::FqRecord, r1_seq_rc::LongDNA{4}, r2_seq_rc::LongDNA{4}; kmer_tolerance::Int64 = 2, overlap_score::Float64 = 0.0, min_ratio_mismatch::Float64 = 0.2, prob_diff::Float64 = 0.0)

    r1_seq = r1.seq
    r2_seq = r2.seq
    r1_length = length(r1_seq)
    r2_length = length(r2_seq)

    # r1_overlap_from, r1_overlap_nmatch, ...
    r1_ms = bitwise_scan(r2_seq_rc, r1_seq, 1, kmer_tolerance)
    # r2_overlap_from, r2_overlap_nmatch, ...
    r2_ms = bitwise_scan(r1_seq_rc, r2_seq, 1, kmer_tolerance)

    # r1_overlap_from == 0 && return false, -1.0
    # r2_overlap_from == 0 && return false, -1.0

    r1_overlap_nbase = r1_length - r1_ms.idx + 1
    r2_overlap_nbase = r2_length - r2_ms.idx + 1
    r1_overlap_nbase != r2_overlap_nbase && return false, -1.0

    if overlap_score > 0
        r1_overlap_prob = probmean(r1, r1_ms.idx, r1_ms.idx + 15)
        r2_overlap_prob = probmean(r2, r2_ms.idx, r2_ms.idx + 15)

        max_overlap_score = max(r1_ms.ncompatible * r1_overlap_prob, r2_ms.ncompatible * r2_overlap_prob)
        max_overlap_score < overlap_score && return false, -1.0
    end

    insert_size = r1_length + r2_length - r1_overlap_nbase

    pe_consensus!(r1, r2, r2_seq_rc, insert_size; min_ratio_mismatch = min_ratio_mismatch, prob_diff = prob_diff)
end

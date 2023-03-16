
struct PrimerSet
    primer1_seqheadset::SeqHeadSet
    primer2_seqheadset::SeqHeadSet
    primer1_rc_seqheadset::SeqHeadSet
    primer2_rc_seqheadset::SeqHeadSet
    length1::Int
    length2::Int
    primer1::LongDNA{4}
    primer2::LongDNA{4}
end

function PrimerSet(primer1::LongDNA{4}, primer2::LongDNA{4})
    primer1_rc = reverse_complement(primer1)
    primer2_rc = reverse_complement(primer2)

    primer1_seqheadset = SeqHeadSet(primer1)
    primer2_seqheadset = SeqHeadSet(primer2)
    primer2_rc_seqheadset = SeqHeadSet(primer1_rc)
    primer2_rc_seqheadset = SeqHeadSet(primer2_rc)
    
    PrimerSet(
        primer1_seqheadset,
        primer2_seqheadset,
        primer2_rc_seqheadset,
        primer2_rc_seqheadset,
        length(primer1),
        length(primer2),
        primer1,
        primer2,
    )
end



@inline function adapter_match_and_trim_pe!(adapter1_seqheadsets::Vector{SeqHeadSet}, adapter2_seqheadsets::Vector{SeqHeadSet},
    adapter1s::Vector{LongDNA{4}}, adapter2s::Vector{LongDNA{4}},
    r1::FqRecord, r2::FqRecord, init_seq_rc::Bool,
    r1_seq_rc::LongDNA{4}, r2_seq_rc::LongDNA{4}, op::PEOptions)

    best_adapter_pe_res = AdapterPERes(false, 9223372036854775807, 9223372036854775807, -1.0, dna"", dna"")
    
    for (adapter1_seqheadset, adapter2_seqheadset, adapter1, adapter2) in zip(adapter1_seqheadsets, adapter2_seqheadsets, adapter1s, adapter2s)
        adapter_pe_res = adapter_match_pe(adapter1_seqheadset, adapter2_seqheadset, r1, r2, adapter1, adapter2, init_seq_rc, r1_seq_rc, r2_seq_rc, op)

        init_seq_rc = false

        best_adapter_pe_res = get_best(best_adapter_pe_res, adapter_pe_res)
        if best_adapter_pe_res.r12_score > 28.8  # 0.9*32
            break
        end
    end
    
    perform_consensus_and_trim!(r1, r2, r1_seq_rc, r2_seq_rc, best_adapter_pe_res, op)
end

@inline function primer_match_pe(ps::PrimerSet,
        r1::FqRecord, r2::FqRecord,
        init_seq_rc::Bool,
        r1_seq_rc::LongDNA{4}, r2_seq_rc::LongDNA{4}, op::PEOptions)
    # r1_seq_rc = r1_seq_rc_threads[thread_id]
    # r2_seq_rc = r2_seq_rc_threads[thread_id]

    # adapter1 = primer2_rc
    # adapter2 = primer1_rc

    # ms = MatchRes
    r1_adapter_ms = bitwise_scan(ps.primer2_rc_seqheadset, r1.seq, 1, op.kmer_tolerance)
    r2_adapter_ms = bitwise_scan(ps.primer1_rc_seqheadset, r2.seq, 1, op.kmer_tolerance)

    r2_seqheadset = SeqHeadSet(r2.seq)
    r1_seqheadset = SeqHeadSet(r1.seq)

    extra_tolerance = max(r1_adapter_ms.ncompatible, r2_adapter_ms.ncompatible) < op.kmer_n_match

    # rx_seq_rc is replaced with reverse complement of rx.seq if init_rc_dest
    r1_pe_ms = bitwise_scan_rc!(r1_seq_rc, r2_seqheadset, r1.seq, 1, op.kmer_tolerance + extra_tolerance, init_rc_dest = init_seq_rc)
    r2_pe_ms = bitwise_scan_rc!(r2_seq_rc, r1_seqheadset, r2.seq, 1, op.kmer_tolerance + extra_tolerance, init_rc_dest = init_seq_rc)

    # ms.idx means insert size, rather than adapter pos
    r1_adapter_ms.idx -= 1
    r2_adapter_ms.idx -= 1


    # if one hit, then check other matches with loosen kmer tomerance
    best_ms = max(r1_adapter_ms, r2_adapter_ms, r1_pe_ms, r2_pe_ms)
    if best_ms.ncompatible > op.kmer_n_match
        max_adapter_pos = best_ms.idx + 1

        if r1_adapter_ms.ncompatible < op.kmer_n_match
            r1_adapter_ms2 = bitwise_scan(adapter1_seqheadset, r1.seq, max_adapter_pos, op.kmer_tolerance+3, until=max_adapter_pos+1)
            if r1_adapter_ms2.ncompatible > op.kmer_n_match
                r1_adapter_ms2.idx -= 1  # insert size
                r1_adapter_ms = r1_adapter_ms2
            end
        end
        if r2_adapter_ms.ncompatible < op.kmer_n_match
            r2_adapter_ms2 = bitwise_scan(adapter2_seqheadset, r2.seq, max_adapter_pos, op.kmer_tolerance+3, until=max_adapter_pos+1)
            if r2_adapter_ms2.ncompatible > op.kmer_n_match
                r2_adapter_ms2.idx -= 1  # insert size
                r2_adapter_ms = r2_adapter_ms2
            end
        end
        if r1_pe_ms.ncompatible < op.kmer_n_match
            r1_rc_pos2 = length(r1.seq) - max_adapter_pos + 2
            r1_pe_ms2 = bitwise_scan(r2_seqheadset, r1_seq_rc, r1_rc_pos2, op.kmer_tolerance+3, until=r1_rc_pos2+1)
            if r1_pe_ms2.ncompatible > op.kmer_n_match
                r1_pe_ms.idx = r1_pe_ms2.idx == 0 ? 0 : length(r1.seq) - r1_pe_ms2.idx + 1
                r1_pe_ms.ncompatible = r1_pe_ms2.ncompatible
            end
        end
        if r2_pe_ms.ncompatible < op.kmer_n_match
            r2_rc_pos2 = length(r2.seq) - max_adapter_pos + 2
            r2_pe_ms2 = bitwise_scan(r1_seqheadset, r2_seq_rc, r2_rc_pos2, op.kmer_tolerance+3, until=r2_rc_pos2+1)
            if r2_pe_ms2.ncompatible > op.kmer_n_match
                r2_pe_ms.idx = r2_pe_ms2.idx == 0 ? 0 : length(r2.seq) - r2_pe_ms2.idx + 1
                r2_pe_ms.ncompatible = r2_pe_ms2.ncompatible
            end
        end
    end

    compute_prob_and_score!(r1_adapter_ms, r1, r1_adapter_ms.idx + 1, r1_adapter_ms.idx + 16)
    compute_prob_and_score!(r2_adapter_ms, r2, r2_adapter_ms.idx + 1, r2_adapter_ms.idx + 16)
    compute_prob_and_score!(r1_pe_ms, r1, r1_pe_ms.idx - 15, r1_pe_ms.idx, r2, 1, 16)
    compute_prob_and_score!(r2_pe_ms, r2, r2_pe_ms.idx - 15, r2_pe_ms.idx, r1, 1, 16)

    # choose possible insert size
    r1_insert_size_real, r1_score = insert_size_decision(r1_adapter_ms.idx, r1_adapter_ms.score, r1_pe_ms.idx, r1_pe_ms.score, insert_size_diff=op.pe_adapter_diff)
    r2_insert_size_real, r2_score = insert_size_decision(r2_adapter_ms.idx, r2_adapter_ms.score, r2_pe_ms.idx, r2_pe_ms.score, insert_size_diff=op.pe_adapter_diff)

    r1_insert_size_decision, r2_insert_size_decision, r12_score = insert_size_decision_separate(r1_insert_size_real, r1_score, r2_insert_size_real, r2_score, insert_size_diff=op.r1_r2_diff)

    # whether this is trimmed. for op.do_read_stats = true
    res_need_trim = false

    # v2.1.0: If a r1/2 adapter is found, but the region of r2/1 is missing or its quality too low (mean prob < 0.6), skip PE check and just trim like single-end

    #TODO: Check r2_insert_size_real or r1_insert_size_decision
    
    if abs(r1_insert_size_real - r2_insert_size_real) > op.r1_r2_diff
        if r1_score > op.trim_score_se
            r2_probs = probmean(r2, r1_insert_size_real + 1, r1_insert_size_real + 16)
            if r2_probs < 0.6
                r2_insert_size_real = r1_insert_size_real
                @goto trim
            end
        elseif r2_score > op.trim_score_se
            # TODO: old is r2_insert_size_pe
            # r1_probs = probmean(r1, r2_insert_size_real + 1, r2_insert_size_pe + 16)

            r1_probs = probmean(r1, r2_insert_size_real + 1, r2_insert_size_real + 16)

            if r1_probs < 0.6
                r1_insert_size_real = r2_insert_size_real
                @goto trim
            end
        end
    end

    if r12_score > op.trim_score
        # < 0: no adapter / pe matched
        is_true_positive = !is_false_positive(
            r1_adapter_ms.idx, r1_pe_ms.idx, length(r1.seq),
            r2_adapter_ms.idx, r2_pe_ms.idx, length(r2.seq),
            insert_size_diff=op.pe_adapter_diff, tail_length=op.tail_length)

        if is_true_positive
            @label trim
            res_need_trim = true
        end
    end
    return AdapterPERes(res_need_trim, r1_insert_size_decision, r2_insert_size_decision, r12_score, adapter1, adapter2)
end

@inline function perform_consensus_and_trim!(
    r1::FqRecord, r2::FqRecord, 
    r1_seq_rc::LongDNA{4}, r2_seq_rc::LongDNA{4},
    res::AdapterPERes, op::PEOptions)
    if res.need_trim
        # consensus
        if op.do_consensus_calling && res.r1_insert_size == res.r2_insert_size && res.r1_insert_size > 0
            is_concensused, ratio_mismatch = pe_consensus!(r1, r2, r2_seq_rc, res.r1_insert_size; min_ratio_mismatch=op.min_ratio_mismatch, prob_diff=op.prob_diff)
        else
            is_concensused = false
            ratio_mismatch = NaN
        end

        # v3.0.0-dev: If choose to trim adapter, check 1 bp offset of adapter sequences. It is because Atria might have 1 bp error in some cases.
        r1_insert_size = if res.r1_insert_size >= 1
            one_bp_check(r1.seq, res.adapter1, res.r1_insert_size, 4)
        else
            res.r1_insert_size
        end
        r2_insert_size = if res.r2_insert_size >= 1
            one_bp_check(r2.seq, res.adapter2, res.r2_insert_size, 4)
        else
            res.r2_insert_size
        end

        if r1_insert_size < 0
            r1_insert_size = 0
        end
        if r2_insert_size < 0
            r2_insert_size = 0
        end

        tail_trim!(r1, r1_insert_size)
        tail_trim!(r2, r2_insert_size)

    else
        if op.do_consensus_calling
            r1_insert_size = r2_insert_size = -1
            is_concensused, ratio_mismatch = pe_consensus!(r1, r2, r1_seq_rc, r2_seq_rc; kmer_tolerance=op.kmer_tolerance_consensus, overlap_score=op.overlap_score, min_ratio_mismatch=op.min_ratio_mismatch, prob_diff=op.prob_diff)
        else
            r1_insert_size = r2_insert_size = -1
            is_concensused = false
            ratio_mismatch = NaN
        end
    end
    return res.need_trim, r1_insert_size, r2_insert_size, is_concensused, ratio_mismatch
end
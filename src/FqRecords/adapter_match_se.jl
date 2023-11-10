


@inline function adapter_match_se(adapter1_seqheadset::SeqHeadSet,
                           r1::FqRecord,
                           kmer_tolerance::Int64,
                           trim_score::Float64)
    r1_adapter_match = bitwise_scan(adapter1_seqheadset, r1.seq, 1, kmer_tolerance)
    compute_prob_and_score!(r1_adapter_match, r1, r1_adapter_match.idx, r1_adapter_match.idx + 15)

    r1_insert_size = r1_adapter_match.idx - 1

    if r1_adapter_match.score > trim_score
        # r1_insert_size can be -1
        # trim
        r1_insert_size < 0 ? 0 : r1_insert_size
    else
        9223372036854775807  # typemax, no trim
    end
end

@inline function adapter_match_se(adapter1_seqheadsets::Vector{SeqHeadSet},
                           r1::FqRecord,
                           kmer_tolerance::Int64,
                           trim_score::Float64)

    nremain = 9223372036854775807  # typemax, no trim
    for adapter1_seqheadset in adapter1_seqheadsets
        nremain_new = adapter_match_se(adapter1_seqheadset, r1, kmer_tolerance, trim_score)
        if nremain_new < nremain
            nremain = nremain_new

            if nremain_new == 0
                break
            end
        end
    end
    nremain
end
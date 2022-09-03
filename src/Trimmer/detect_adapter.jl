
const adapter_candidates = [
    LongDNA{4}("AAGTCGGAGGCCAAGC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AAGTCGGATCGTAGCC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AATGATACGGCGACCA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("ACACTCTTTCCCTACA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AGATCGGAAGAGCACA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AGATCGGAAGAGCGGT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AGATCGGAAGAGCGTC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("AGATCGGAAGAGCTCG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CAAGCAGAAGACGGCA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCACTACGCCTCCGCT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCGACAGGTTCAGAGT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCGAGCCCACGAGACA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCGAGCCCACGAGACC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCGAGCCCACGAGACG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CCGAGCCCACGAGACT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CGACAGGTTCAGAGTT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CGGTCTCGGCATTCCT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CTAATACGACTCACTA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CTGAGCGGGCTGGCAA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CTGATGGCGCGAGGGA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CTGCCCCGGGTTCCTC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("CTGTCTCTTATACACA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGAAC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGAAG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGAAT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGACG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGACT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGAGC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGATA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GACGCTGCCGACGATC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GATCGGAAGAGCACAC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GATCGGAAGAGCGGTT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GATCGGAAGAGCGTCG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GATCGGAAGAGCTCGT") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GATCGTCGGACTGTAG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GTCTCGTGGGCTCGGA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("GTGACTGGAGTTCAGA") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("TACACTCTTTCCCTAC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("TCGGACTGTAGAACTC") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("TCGTCGGCAGCGTCAG") |> bitsafe! |> SeqHeadSet,
    LongDNA{4}("TGGAATTCTCGGGTGC") |> bitsafe! |> SeqHeadSet,
]

"""
    count_adapter(rs::Vector{FqRecord}, adapter::SeqHeadSet; kmer_tolerance::Int = 2)

Return `adapter_sequence::String, adapter_count::Int, average_identity::Float`
"""
function count_adapter(rs::Vector{FqRecord}, adapter::SeqHeadSet; kmer_tolerance::Int = 2)
    n_adapter = 0
    base_match = 0
    nmatch_threshold = 16 - kmer_tolerance
    for r in rs
        pos, nmatch = bitwise_scan(adapter, r.seq, 1, kmer_tolerance)
        if nmatch >= nmatch_threshold
            n_adapter += 1
            base_match += nmatch
        end
    end
    String(LongDNA{4}(adapter)), n_adapter, base_match / n_adapter / 16
end

function detect_adapter_threads!(rs::Vector{FqRecord}; kmer_tolerance::Int = 2, adapter_candidates::Vector{SeqHeadSet} = adapter_candidates)
    tasks = map(adapter_candidates) do adapter
        Threads.@spawn count_adapter(rs, adapter; kmer_tolerance=kmer_tolerance)
    end

    adapters = String[]
    counts = Int[]
    identities = Float64[]
    for task in tasks
        adapter, count, base_match = fetch(task)
        push!(adapters, adapter)
        push!(counts, count)
        push!(identities, base_match)
    end
    df = [adapters counts identities]
    perm = sortperm(counts, rev=true)
    df_sorted = @inbounds df[perm,:]
    top5 = @inbounds df_sorted[1:5,:]
    headers = ["Adapter", "Occurance", "Identity"]
    return top5, headers
end

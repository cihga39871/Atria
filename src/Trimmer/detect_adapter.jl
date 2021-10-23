
const adapter_candidates = [
    LongDNASeq("AAGTCGGAGGCCAAGC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AAGTCGGATCGTAGCC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AATGATACGGCGACCA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("ACACTCTTTCCCTACA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AGATCGGAAGAGCACA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AGATCGGAAGAGCGGT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AGATCGGAAGAGCGTC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("AGATCGGAAGAGCTCG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CAAGCAGAAGACGGCA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCACTACGCCTCCGCT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCGACAGGTTCAGAGT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCGAGCCCACGAGACA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCGAGCCCACGAGACC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCGAGCCCACGAGACG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CCGAGCCCACGAGACT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CGACAGGTTCAGAGTT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CGGTCTCGGCATTCCT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CTAATACGACTCACTA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CTGAGCGGGCTGGCAA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CTGATGGCGCGAGGGA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CTGCCCCGGGTTCCTC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("CTGTCTCTTATACACA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGAAC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGAAG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGAAT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGACG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGACT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGAGC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGATA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GACGCTGCCGACGATC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GATCGGAAGAGCACAC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GATCGGAAGAGCGGTT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GATCGGAAGAGCGTCG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GATCGGAAGAGCTCGT") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GATCGTCGGACTGTAG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GTCTCGTGGGCTCGGA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("GTGACTGGAGTTCAGA") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("TACACTCTTTCCCTAC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("TCGGACTGTAGAACTC") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("TCGTCGGCAGCGTCAG") |> bitsafe! |> SeqHeadSet,
    LongDNASeq("TGGAATTCTCGGGTGC") |> bitsafe! |> SeqHeadSet,
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
    String(LongDNASeq(adapter)), n_adapter, base_match / n_adapter / 16
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
    df_sorted = sortslices(df, dims=1, lt=(x,y)->!isless(x[2],y[2]))
    top3 = ["Adapter" "Occurance" "Identity"; df_sorted[1:3,:]]
end

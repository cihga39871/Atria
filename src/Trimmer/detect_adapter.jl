
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
        ms = bitwise_scan(adapter, r.seq, 1, kmer_tolerance)
        if ms.ncompatible >= nmatch_threshold
            n_adapter += 1
            base_match += ms.ncompatible
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

const empty_seq = dna""

"""
    check_pe_match(r1, r2)

Match pair-end to infer adapter sequences
"""
function check_pe_match(r1::FqRecord, r2::FqRecord; kmer_tolerance = 3, kmer_n_match = 9)
    r2_seqheadset = SeqHeadSet(r2.seq)
    r1_seqheadset = SeqHeadSet(r1.seq)

    ms1 = bitwise_scan_rc(r2_seqheadset, r1.seq, 1, kmer_tolerance)
    ms2 = bitwise_scan_rc(r1_seqheadset, r2.seq, 1, kmer_tolerance)

    r1_insert_size_pe = ms1.idx
    r1_pe_nmatch = ms1.ncompatible
    r2_insert_size_pe = ms2.idx
    r2_pe_nmatch = ms2.ncompatible
    # if r1_insert_size_pe != r2_insert_size_pe
    #     return empty_seq, 0.0, empty_seq, 0.0, 0
    # end
    
    adapter_start = r1_insert_size_pe + 1

    if r1_pe_nmatch < kmer_n_match
        adapter1 = empty_seq
        r1_adapter_prob = 0.0
    else
        adapter_end_r1 = min(adapter_start + 15, length(r1.seq))
        adapter1 = r1.seq[adapter_start:adapter_end_r1]
        r1_adapter_prob = probmean(r1, adapter_start, adapter_end_r1)
    end

    if r2_pe_nmatch < kmer_n_match
        adapter2 = empty_seq
        r2_adapter_prob = 0.0
    else
        adapter_end_r2 = min(adapter_start + 15, length(r2.seq))
        adapter2 = r2.seq[adapter_start:adapter_end_r2]
        r2_adapter_prob = probmean(r2, adapter_start, adapter_end_r2)
    end

    return adapter1, r1_adapter_prob, adapter2, r2_adapter_prob, adapter_start
end

function check_pe_match(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}; kmer_tolerance = 3, kmer_n_match = 9, occurance = 0.0004)
    adapter1s = LongDNA{4}[]
    adapter2s = LongDNA{4}[]
    r1_adapter_probs = Float64[]
    r2_adapter_probs = Float64[]
    r1_adapter_starts = Int[]
    r2_adapter_starts = Int[]
    for (r1,r2) in zip(r1s, r2s)
        adapter1, r1_adapter_prob, adapter2, r2_adapter_prob, adapter_start = check_pe_match(r1, r2; kmer_tolerance = kmer_tolerance, kmer_n_match = kmer_n_match)
        if adapter_start < 16
            continue
        end

        if length(adapter1) > 0 && r1_adapter_prob > 0.75
            push!(adapter1s, adapter1)
            push!(r1_adapter_probs, r1_adapter_prob)
            push!(r1_adapter_starts, adapter_start)
        end
        if length(adapter2) > 0 && r2_adapter_prob > 0.75
            push!(adapter2s, adapter2)
            push!(r2_adapter_probs, r2_adapter_prob)
            push!(r2_adapter_starts, adapter_start)
        end
    end

    nread = length(r1s)
    r1_df = DataFrame(:adapter => adapter1s, :adapter_prob => r1_adapter_probs, :adapter_start => r1_adapter_starts)
    r1_stats = read_adapter_stats(r1_df, nread; occurance = occurance)

    r2_df = DataFrame(:adapter => adapter2s, :adapter_prob => r2_adapter_probs, :adapter_start => r2_adapter_starts)
    r2_stats = read_adapter_stats(r2_df, nread; occurance = occurance)

    return r1_stats, r2_stats
end

function read_adapter_stats(r_df::DataFrame, nread::Int; occurance::Float64 = 0.0004)
    gdf = groupby(r_df, :adapter)
    stats = combine(gdf, :adapter => (x -> x[1]) => :adapter, nrow => :count, :adapter_prob => mean => :identity, :adapter_start => mean => :adapter_start_position)
    
    stats.occurance = stats.count ./ nread
    select!(stats, :adapter, :count, :occurance, :identity, :adapter_start_position)

    if nrow(stats) == 0
        return stats
    end

    merge_short_to_long_adapters!(stats, nread)

    filter!(:occurance => x -> x > occurance, stats)

    if nrow(stats) == 0
        return stats
    end

    sort!(stats, [:count, :identity], rev=true)
    count_cutoff = stats[1, :count] * 0.3

    filter!(:count => x -> x > count_cutoff, stats)
end

"""
    merge_short_to_long_adapters!(stats::DataFrame, nread::Int)
```
stats = 
 Row │ adapter           count  occurance   identity  adapter_start_position 
     │ LongSequ…         Int64  Float64     Float64   Float64                
─────┼───────────────────────────────────────────────────────────────────────
   1 │ AGATCGGAAGAGCGTC   3025  0.0144048   0.997003                 82.7977
   2 │ AGATCGGAAGAGCG      877  0.00417619  0.997014                 86.3683
   3 │ AGATCGGAAG          815  0.00388095  0.996967                 89.4417
   4 │ AGATCGGAAGAG        810  0.00385714  0.997015                 87.884
   5 │ AGATCG              748  0.0035619   0.996916                 91.1283
   6 │ AGATCGGA            709  0.00337619  0.997042                 90.2976
   7 │ AGAT                659  0.0031381   0.996979                 91.8285
   8 │ AG                  659  0.0031381   0.996918                 91.6464
   9 │ A                   480  0.00228571  0.997028                 90.0083
  10 │ AGA                 451  0.00214762  0.997028                 89.1242
  11 │ AGATC               398  0.00189524  0.997129                 88.3065
  12 │ AGATCGG             349  0.0016619   0.997093                 86.9427
```
"""
function merge_short_to_long_adapters!(stats::DataFrame, nread::Int)
    # stats_backup = deepcopy(stats)
    sort!(stats, :adapter, rev=true)
    for i in nrow(stats):-1:1
        seq = stats.adapter[i]
        if length(seq) == 16  # full length, cannot have parent adapters
            continue
        end
        seq_regex = BioRegex{DNA}("^" * String(stats.adapter[i]))
        parent_indices = map(x -> occursin(seq_regex, x), stats.adapter)
        parent_indices[i] = false
        n_parent = sum(parent_indices)
        if n_parent == 0  # skip
            continue
        end
        occurance_each = stats.occurance[i] / n_parent
        identity_each = stats.identity[i] * occurance_each
        adapter_start_position_each = stats.adapter_start_position[i] * occurance_each
        for (j, is_parent) in enumerate(parent_indices)
            is_parent || continue
            # weighted average for identity and adapter start position
            stats.identity[j] = (stats.identity[j] * stats.occurance[j] + identity_each) / (stats.occurance[j] + occurance_each)
            stats.adapter_start_position[j] = (stats.adapter_start_position[j] * stats.occurance[j] + adapter_start_position_each) / (stats.occurance[j] + occurance_each)
            
            stats.occurance[j] += occurance_each
        end
        deleteat!(stats, i)
    end

    # after all occurance calculation, compute count because of Int inaccuracy.
    stats.count = round.(Int, stats.occurance .* nread)
    stats
end

function show_paired_adapter_result(file::AbstractString, r_stats::DataFrame, n_reads::Int)
    print("$file: ")
    
    if nrow(r_stats) > 0
        println(r_stats)
    else
        println("no adapter found in the first $n_reads reads.")
    end
    println()
end

init_adapter_detection_summary() = DataFrame(
    :r1 => String[],
    :best_adapter1 => LongDNA{4}[],
    :occurance1 => Float64[],
    :identity1 => Float64[],
    :adapter_start_position1 => Float64[],
    :r2 => String[],
    :best_adapter2 => LongDNA{4}[],
    :occurance2 => Float64[],
    :identity2 => Float64[],
    :adapter_start_position2 => Float64[],
)

function push_adapter_detection_summary!(adapter_detection_summary::DataFrame, file1::String, r1_stats::DataFrame, file2::String, r2_stats::DataFrame)
    if nrow(r1_stats) == 0
        best_adapter1 = empty_seq
        occurance1 = NaN
        identity1 = NaN
        adapter_start_position1 = NaN
    else
        best_adapter1 = r1_stats.adapter[1]
        occurance1 = r1_stats.occurance[1]
        identity1 = r1_stats.identity[1]
        adapter_start_position1 = r1_stats.adapter_start_position[1]
    end
    if nrow(r2_stats) == 0
        best_adapter2 = empty_seq
        occurance2 = NaN
        identity2 = NaN
        adapter_start_position2 = NaN
    else
        best_adapter2 = r2_stats.adapter[1]
        occurance2 = r2_stats.occurance[1]
        identity2 = r2_stats.identity[1]
        adapter_start_position2 = r2_stats.adapter_start_position[1]
    end
    push!(adapter_detection_summary, [
        file1, best_adapter1, occurance1, identity1, adapter_start_position1,
        file2, best_adapter2, occurance2, identity2, adapter_start_position2
    ])
end
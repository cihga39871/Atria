

const qualpval_table = map(q -> 10 ^ (-q/10), 0:50)
const qualprob_table = 1.0 .- qualpval_table

"""
The quality offset used in FqRecords should be the real quality offset - 1, such as Illumina 1.8 => 33-1
"""
@inline function qualpval(Q, quality_offset)::Float64
    q = Q - quality_offset + 1
    q <= 0 && error("Input quality < 0 detected. Wrong --quality-format FORMAT or the input file in truncated.")
    @inbounds qualpval_table[q > 51 ? 51 : q]
end

"""
The quality offset used in FqRecords should be the real quality offset - 1, such as Illumina 1.8 => 33-1
"""
@inline function qualprob(Q, quality_offset)::Float64
    q = Q - quality_offset + 1
    q <= 0 && error("Input quality < 0 detected. Wrong --quality-format FORMAT or the input file in truncated.")
    @inbounds qualprob_table[q > 51 ? 51 : q]
end


@inline function update_prob_from_qual(r::FqRecord; quality_offset::Int64=33)::Nothing
    resize!(r.prob, length(r.qual))
    @inbounds for (i,Q) in enumerate(r.qual)
        r.prob[i] = qualprob(Q, quality_offset)
    end
    return
end



@inline function probsum(r::FqRecord, from::Int64, to::Int64)::Float64
    r_prob = r.prob
    nprob = length(r_prob)
    to   > nprob && (to   = nprob)
    from < 1     && (from = 1    )

    value = 0.0
    @inbounds while from <= to
        value += r_prob[from]
        from += 1
    end
    value
end

@inline function probmean(r::FqRecord, from::Int64, to::Int64)::Float64
    r_prob = r.prob
    nprob = length(r_prob)
    to   > nprob && (to   = nprob)
    from < 1     && (from = 1    )
    n = to - from + 1
    n <= 0 && return 0.0

    value = 0.0
    @inbounds while from <= to
        value += r_prob[from]
        from += 1
    end
    @fastmath value/n
end



@inline function compute_prob_and_score!(match_res::MatchRes, r::FqRecord, r_start::Int, r_end::Int; min_prob::Float64 = 0.75)
    match_res.prob = max(probmean(r, r_start, r_end), min_prob)
    match_res.score = @fastmath match_res.ncompatible * match_res.prob
end
@inline function compute_prob_and_score!(match_res::MatchRes, r1::FqRecord, r1_start::Int, r1_end::Int, r2::FqRecord, r2_start::Int, r2_end::Int; min_prob::Float64 = 0.75)
    prob1 = max(probmean(r1, r1_start, r1_end), min_prob)
    prob2 = max(probmean(r2, r2_start, r2_end), min_prob)
    match_res.prob = @fastmath prob1 * prob2
    match_res.score = @fastmath match_res.ncompatible * match_res.prob
end
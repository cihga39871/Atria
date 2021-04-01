
# @inline qualpval(Q, quality_offset) = 10 ^ ( -(Q - quality_offset) / 10) # ^ operation is slow!
# @inline qualprob(Q, quality_offset) = 1 - qualpval(Q, quality_offset)


const qualpval_table = map(q -> 10 ^ (-q/10), 0:50)
const qualprob_table = 1.0 .- qualpval_table
@inline function qualpval(Q, quality_offset)::Float64
    q = Q - quality_offset
    q < 0 && error("Input quality < 0 detected. Wrong --quality-format FORMAT or the input file in truncated.")
    @inbounds qualpval_table[q > 51 ? 51 : q]
end
@inline function qualprob(Q, quality_offset)::Float64
    q = Q - quality_offset
    q < 0 && error("Input quality < 0 detected. Wrong --quality-format FORMAT or the input file in truncated.")
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

# function qualcounttable(rs::Vector{FqRecord}, quality_offset::UInt8=0x21)
#     max_len = 0
#     for r in rs
#         len = length(r.qual)
#         if len > max_len
#             max_len = len
#         end
#     end
#
#     qualtable = zeros(Int, max_len, 30)
#     for r in rs
#         len = length(r.qual)
#         @inbounds for i = 1:len
#             q = min(r.qual[i] - quality_offset, 30)
#             qualtable[i, q] += 1
#         end
#     end
#     qualtable
# end

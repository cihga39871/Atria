
mutable struct DupCount
    @atomic count::Int
    id::String
end


function write_pcr_dedup_count(out_pcr_dedup_count::AbstractString, dup_dict::Dict{Tuple{LongDNA{4}, LongDNA{4}}, DupCount})
    dup_count = 0
    open(out_pcr_dedup_count, "w+") do io
        println(io, "count\tid\tr1_raw\tr2_raw")
        for (key, v) in dup_dict
            if v.count > 1
                @inbounds println(io, "$(v.count)\t$(v.id)\t$(key[1])\t$(key[2])")
                dup_count += v.count
            end
        end
    end
    dup_count
end

function write_pcr_dedup_count(out_pcr_dedup_count::AbstractString, dup_dict::Dict{LongDNA{4}, DupCount})
    dup_count = 0
    open(out_pcr_dedup_count, "w+") do io
        println(io, "count\tid\tread_raw")
        for (key, v) in dup_dict
            if v.count > 1
                println(io, "$(v.count)\t$(v.id)\t$(key)")
                dup_count += v.count
            end
        end
    end
    dup_count
end

function get_dup_count(dup_dict::Dict)
    dup_count = 0
    for v in values(dup_dict)
        if v.count > 1
            dup_count += v.count
        end
    end
    dup_count

end


# function pcr_dedup(dup_dict::Dict{Tuple{LongDNA{4}, LongDNA{4}}, DupCount}, r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
#     n_dup = 0

#     @inbounds for i in 1:n_reads
#         r1 = r1s[i]
#         r2 = r2s[i]
#         key = (r1.seq, r2.seq)
#         if haskey(dup_dict, key)
#             dup_dict[key].count += 1
#             isgoods[i] = false
#             n_dup += 1
#         else
#             new_key = (copy(r1.seq), copy(r2.seq))
#             dup_dict[new_key] = DupCount(1, String(copy(r1.id)))
#         end
#     end
#     n_dup
# end

# function pcr_dedup(dup_dict::Dict{LongDNA{4}, DupCount}, r1s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
#     n_dup = 0
#     @inbounds for i in 1:n_reads
#         r1 = r1s[i]
#         if haskey(dup_dict, r1.seq)
#             dup_dict[r1.seq].count += 1
#             isgoods[i] = false
#             n_dup += 1
#         else
#             new_key = copy(r1.seq)
#             dup_dict[new_key] = DupCount(1, String(copy(r1.id)))
#         end
#     end
#     n_dup
# end
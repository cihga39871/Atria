
function processing_reads!(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
    if length(isgoods) < n_reads
        resize!(isgoods, n_reads)
    end
    for i in 1:n_reads
        @inbounds isgoods[i] = read_processing!(r1s[i]::FqRecord, r2s[i]::FqRecord, 1)
    end
    nothing
end

# function processing_reads_range!(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, reads_range::UnitRange{Int64})
#     this_threadid = Threads.threadid()
#     for i in reads_range
#         @inbounds isgoods[i] = read_processing!(r1s[i], r2s[i], this_threadid)
#     end
#     nothing
# end

# function processing_reads_threads!(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
#     if length(isgoods) < n_reads
#         resize!(isgoods, n_reads)
#     end
#     # split reads to N reads per batch
#     @sync for reads_start in 1:256:n_reads
#         reads_end = min(reads_start + 255, n_reads)
#         reads_range = reads_start:reads_end

#         Threads.@spawn processing_reads_range!(r1s::Vector{FqRecord}, r2s::Vector{FqRecord}, isgoods::Vector{Bool}, reads_range)
#     end
#     nothing
# end

## single end
function processing_reads!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
    if length(isgoods) < n_reads
        resize!(isgoods, n_reads)
    end
    for i in 1:n_reads
        @inbounds isgoods[i] = read_processing!(r1s[i]::FqRecord, 1)
    end
    nothing
end

# function processing_reads_range!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, reads_range::UnitRange{Int64})
#     this_threadid = Threads.threadid()
#     for i in reads_range
#         @inbounds isgoods[i] = read_processing!(r1s[i], this_threadid)
#     end
#     nothing
# end

# function processing_reads_threads!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, n_reads::Int)
#     if length(isgoods) < n_reads
#         resize!(isgoods, n_reads)
#     end
#     # split reads to N reads per batch
#     @sync for reads_start in 1:512:n_reads
#         reads_end = min(reads_start + 511, n_reads)
#         reads_range = reads_start:reads_end

#         Threads.@spawn processing_reads_range!(r1s::Vector{FqRecord}, isgoods::Vector{Bool}, reads_range)
#     end
#     nothing
# end

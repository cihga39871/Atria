
function Base.copy(r::FqRecord)
    id = Vector{UInt8}(undef, length(r.id))
    des = Vector{UInt8}(undef, length(r.des))
    qual = Vector{UInt8}(undef, length(r.qual))
    prob = Vector{Float64}(undef, length(r.prob))

    copyto!(id, 1, r.id, 1, length(r.id))
    copyto!(des, 1, r.des, 1, length(r.des))
    copyto!(qual, 1, r.qual, 1, length(r.qual))
    copyto!(prob, 1, r.prob, 1, length(r.prob))

    seq = copy(r.seq)

    FqRecord(id, seq, des, qual, prob)
end

@inline function safe_copyto!(dest::Vector{UInt8}, src::T) where T <: AbstractArray
    resize!(dest, length(src))
    copyto!(dest, src)
end

@inline function safe_copyto!(dest::Vector{T}, src::Vector{T}, src_offset, N) where T <: Any
    resize!(dest, N)
    unsafe_copyto!(dest, 1, src, src_offset, N)
end

@inline function safe_copyto!(dest::LongDNASeq, src::Vector{UInt8}, src_offset, N)
    resize!(dest, N)
    # BioSequences.encode_chunks!(dest, 1, src, src_offset, N)
    copyto!(dest, 1, src, src_offset, N)
end
@inline function safe_copyto!(dest::LongDNASeq, src::Vector{UInt8})
    copy!(dest, src)
end

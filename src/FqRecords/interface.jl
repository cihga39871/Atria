
struct FqRecord
    id::Vector{UInt8}
    seq::LongDNASeq
    des::Vector{UInt8}
    qual::Vector{UInt8}
    prob::Vector{Float64}
    function FqRecord(id::Vector{UInt8}, seq::LongDNASeq, des::Vector{UInt8}, qual::Vector{UInt8}, prob::Vector{Float64})
        new(id::Vector{UInt8}, seq::LongDNASeq |> bitsafe!, des::Vector{UInt8}, qual::Vector{UInt8}, prob::Vector{Float64})
    end
end
FqRecord() = FqRecord(Vector{UInt8}(), LongDNASeq(), Vector{UInt8}(), Vector{UInt8}(), Vector{Float64}())

@inline function FqRecord(id::Vector{UInt8}, seq::LongDNASeq, des::Vector{UInt8}, qual::Vector{UInt8}; quality_offset=33)
    FqRecord(id, seq, des, qual, qualprob.(qual, quality_offset))
end

Base.:(==)(r1::FqRecord, r2::FqRecord) =
    r1.id == r2.id && r1.seq == r2.seq && r1.des == r2.des && r1.qual == r2.qual

function Base.isempty(r::FqRecord)::Bool
    isempty(r.id::Vector{UInt8}) && isempty(r.seq::LongDNASeq) && isempty(r.des::Vector{UInt8}) && isempty(r.qual::Vector{UInt8})
end

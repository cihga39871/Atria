
struct FqRecord
    id::Vector{UInt8}
    seq::LongDNA{4}
    des::Vector{UInt8}
    qual::Vector{UInt8}
    prob::Vector{Float64}
    function FqRecord(id::Vector{UInt8}, seq::LongDNA{4}, des::Vector{UInt8}, qual::Vector{UInt8}, prob::Vector{Float64})
        new(id::Vector{UInt8}, seq::LongDNA{4} |> bitsafe!, des::Vector{UInt8}, qual::Vector{UInt8}, prob::Vector{Float64})
    end
end
@inline FqRecord() = FqRecord(Vector{UInt8}(), LongDNA{4}(), Vector{UInt8}(), Vector{UInt8}(), Vector{Float64}())

@inline function FqRecord(id::Vector{UInt8}, seq::LongDNA{4}, des::Vector{UInt8}, qual::Vector{UInt8}; quality_offset=33)
    FqRecord(id, seq, des, qual, qualprob.(qual, quality_offset))
end

@inline Base.:(==)(r1::FqRecord, r2::FqRecord) =
    r1.id == r2.id && r1.seq == r2.seq && r1.des == r2.des && r1.qual == r2.qual

@inline function Base.isempty(r::FqRecord)::Bool
    isempty(r.id::Vector{UInt8}) && isempty(r.seq::LongDNA{4}) && isempty(r.des::Vector{UInt8}) && isempty(r.qual::Vector{UInt8})
end

mutable struct TrimStats
    @atomic polyG::Int
    @atomic polyT::Int
    @atomic polyA::Int
    @atomic polyC::Int
    @atomic complexity_filtered::Int
    @atomic hard_clip_after::Int
    @atomic tail_low_qual_trim::Int
    @atomic tail_N_trim::Int
    @atomic length_filtered::Int
    @atomic max_n_filtered::Int
    @atomic pcr_dedup_removed::Int
    @atomic quality_trim::Int
    @atomic adapter_trim::Int
end
TrimStats() = TrimStats(0,0,0,0,0,0,0,0,0,0,0,0,0)

function Base.empty!(t::TrimStats)
    @atomic t.polyG = 0
    @atomic t.polyT = 0
    @atomic t.polyA = 0
    @atomic t.polyC = 0
    @atomic t.complexity_filtered = 0
    @atomic t.hard_clip_after = 0
    @atomic t.tail_low_qual_trim = 0
    @atomic t.tail_N_trim = 0
    @atomic t.length_filtered = 0
    @atomic t.max_n_filtered = 0
    @atomic t.pcr_dedup_removed = 0
    @atomic t.quality_trim = 0
    @atomic t.adapter_trim = 0
end
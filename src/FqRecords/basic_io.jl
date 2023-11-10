
function Base.println(io::IO, r::FqRecord)
    println(String(copy(r.id)))
    println(r.seq)
    println(String(copy(r.des)))
    println(String(copy(r.qual)))
end

Base.print(io::IO, r::FqRecord) = println(io, r)
Base.display(r::FqRecord) = println(stdout, r)
Base.show(r::FqRecord) = println(stdout, r)


remove_blank(s::AbstractString) = replace(s, r"\n[ \t]*" => "\n")

"""
    fqreadrecord(s::IO    ; quality_offset=33)::FqRecord
    fqreadrecord(s::String; quality_offset=33)::FqRecord

It is very slow and not recommended. See also `load_fqs_threads!`.
"""
function fqreadrecord(s::IO; quality_offset=33)::FqRecord
    # 0x0a is \n
    # do not compatible with \r\n
    id = readuntil(s, 0x0a, keep=false)::Vector{UInt8}
    seq = LongDNA{4}(readuntil(s, 0x0a, keep=false))::LongDNA{4}
    des = readuntil(s, 0x0a, keep=false)::Vector{UInt8}
    qual = readuntil(s, 0x0a, keep=false)::Vector{UInt8}
    # nqual = length(qual::Vector{UInt8})::Int64
    FqRecord(id, seq, des, qual; quality_offset=quality_offset)::FqRecord
end
fqreadrecord(s::String; quality_offset=33)::FqRecord = fqreadrecord(IOBuffer(remove_blank(s)), quality_offset=quality_offset)

"""
    fqreadrecord!(r::FqRecord, s::IO)

It is very slow and not recommended. See also `load_fqs_threads!`.
"""
function fqreadrecord!(r::FqRecord, s::IO; quality_offset=33)
    safe_copyto!(r.id, readuntil(s, 0x0a, keep=false)::Vector{UInt8})
    safe_copyto!(r.seq, readuntil(s, 0x0a, keep=false)::Vector{UInt8})
    bitsafe!(r.seq)
    safe_copyto!(r.des, readuntil(s, 0x0a, keep=false)::Vector{UInt8})
    resize!(r.qual, length(r.seq))
    readfill!(s, r.qual)

    resize!(r.prob, length(r.seq))
    @inbounds for (i, q) in enumerate(r.qual)
        r.prob[i] = qualprob(q, quality_offset)
    end

    eof(s) && return
    read(s, UInt8) == 0xa || error("FASTQ is not valid: the lengths of sequence and quality are not the same for $(string(r.id)): $s")
    return
end


function fqwriterecord(io::IO, r::FqRecord)
    if isempty(r.seq::LongDNA{4})
        write(io, r.id::Vector{UInt8})
        write(io, '\n')
        write(io, 'N')
        write(io, '\n')
        write(io, r.des::Vector{UInt8})
        write(io, '\n')
        write(io, '!')
        write(io, '\n')
    else
        write(io, r.id::Vector{UInt8})
        write(io, '\n')
        print(io, r.seq::LongDNA{4}) # no write method for LongDNA{4}
        write(io, '\n')
        write(io, r.des::Vector{UInt8})
        write(io, '\n')
        write(io, r.qual::Vector{UInt8})
        write(io, '\n')
    end
end

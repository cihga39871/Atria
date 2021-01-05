
"""
    write_fqs_threads!(io1out::IO, io2out::IO,
        outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
        r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
        n_reads::Int, range_filter)

The interface to write paired FASTQ reads.
"""
function write_fqs_threads!(io1out::IO, io2out::IO,
    outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
    r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
    n_reads::Int, range_filter)

    FqRecord2StringVec!(outr1s::Vector{Vector{UInt8}}, r1s::Vector{FqRecord}, n_reads::Int)
    task_write = Threads.@spawn writebytes(io1out, outr1s, range_filter)

    FqRecord2StringVec!(outr2s::Vector{Vector{UInt8}}, r2s::Vector{FqRecord}, n_reads::Int)
    wait(task_write)
    writebytes(io2out, outr2s, range_filter)
end


"""
    FqRecord2StringVec!(out::Vector{UInt8}, r::FqRecord)

Empty `out`, and then convert `r` to it continuously. If empty sequence, write a N as sequence and a ! as quality.
"""
@inline function FqRecord2StringVec!(out::Vector{UInt8}, r::FqRecord)::Nothing
    # out = Base.StringVector(0)
    empty!(out)
    if isempty(r.seq::LongDNASeq)
        append!(out, r.id::Vector{UInt8})
        push!(out, 0x0a)  # \n
        push!(out, 0x4e)  # N
        push!(out, 0x0a)
        append!(out, r.des::Vector{UInt8})
        push!(out, 0x0a)
        push!(out, 0x21)  # !
        push!(out, 0x0a)
    else
        append!(out, r.id::Vector{UInt8})
        push!(out, 0x0a)
        length_out = length(out)
        r_seq = r.seq
        length_r = length(r_seq)
        resize!(out, length_out + length_r)
        @inbounds for (i, base) in enumerate(r_seq)
            out[length_out + i] = UInt8(convert(Char, base))
        end
        # slow:
        # for i in r.seq
        #     push!(out, UInt8(convert(Char, i)))
        # end
        push!(out, 0x0a)
        append!(out, r.des::Vector{UInt8})
        push!(out, 0x0a)
        append!(out, r.qual::Vector{UInt8})
        push!(out, 0x0a)
    end
    nothing
end

"""
    FqRecord2StringVec!(outrs::Vector{Vector{UInt8}}, rs::Vector{FqRecord}, stop::Int)
    FqRecord2StringVec!(outrs::Vector{Vector{UInt8}}, rs::Vector{FqRecord}, reads_range::UnitRange)

- `outrs`: the vector of string vectors to be modified in place.

- `rs`: the vector of reads to be converted.

- `stop::Int`: only convert `rs` in the range of `1:stop`.

- `reads_range::UnitRange`: only convert `rs` in the reads range.
"""
@inline function FqRecord2StringVec!(outrs::Vector{Vector{UInt8}}, rs::Vector{FqRecord}, stop::Int)::Nothing
    n_outrs = length(outrs)
    if n_outrs < stop
        # make outrs larger
        append!(outrs,
            Vector{UInt8}[Base.StringVector(0) for i = 1:(stop-n_outrs)])
    end
    @sync for reads_start in 1:4096:stop
        reads_end = min(reads_start + 4095, stop)
        reads_range = reads_start:reads_end
        Threads.@spawn FqRecord2StringVec!(outrs, rs, reads_range)
    end
    nothing
end

@inline function FqRecord2StringVec!(outrs::Vector{Vector{UInt8}}, rs::Vector{FqRecord}, reads_range::UnitRange)::Nothing
    @inbounds for i in reads_range
        FqRecord2StringVec!(outrs[i], rs[i])
    end
    nothing
end


@inline function writebytes(io::IOStream, outrs::Vector{Vector{UInt8}}, stop::Int)::Nothing
    @inbounds for i in 1:stop
        write(io, outrs[i])
    end
    nothing
end
@inline function writebytes(io::IOStream, outrs::Vector{Vector{UInt8}}, filters::SubArray{Bool,1,Array{Bool,1},Tuple{UnitRange{Int64}},true})::Nothing
    @inbounds for (i, val) in enumerate(filters)
        if val
            write(io, outrs[i])
        end
    end
    nothing
end

@inline function writebytes(io::IO, outrs::Vector{Vector{UInt8}}, stop::Int)::Nothing
    # for CodecZlib streams, call write once to increase speed (3.2X)
    # it is even faster than call pigz in shell.
    v_all = @inbounds outrs[1]
    @inbounds for i in 2:stop
        append!(v_all, outrs[i])
    end
    write(io, v_all)
    nothing
end

bytes_tmp = Vector{UInt8}(undef, 67108864) # 2^26

@inline function writebytes(io::IO, outrs::Vector{Vector{UInt8}}, filters::SubArray{Bool,1,Array{Bool,1},Tuple{UnitRange{Int64}},true}, bytes_tmp::Vector{UInt8} = bytes_tmp)::Nothing
    # for CodecZlib streams, call write once to increase speed

    # method 1: 112s
    # v_all = Base.StringVector(0)
    # @inbounds for (i, val) in enumerate(filters)
    #     if val
    #         append!(v_all, outrs[i])
    #     end
    # end
    # write(io, v_all)

    # method2: 213s
    # @inbounds for (i, val) in enumerate(filters)
    #     if val
    #         write(io, outrs[i])
    #     end
    # end

    # method3: 99s # the same speed as the natual pigz
    # bytes_tmp = Vector{UInt8}(undef, 67108864) # 2^26
    start = 1
    stop = length(bytes_tmp)
    @inbounds for (i, val) in enumerate(filters)
        if val
            outr = outrs[i]::Vector{UInt8}
            ncopy = length(outr)::Int
            new_stop = (start + ncopy - 1)::Int
            if new_stop > stop
                stop = max(stop + 2097152, new_stop)::Int
                resize!(bytes_tmp, stop)
            end
            unsafe_copyto!(bytes_tmp, start, outr, 1, ncopy)
            start += ncopy
        end
    end
    p_bytes_tmp = pointer(bytes_tmp)
    bytes_tmp_to_write = unsafe_wrap(Vector{UInt8}, p_bytes_tmp, start-1)
    write(io, bytes_tmp_to_write)
    nothing
end

@inline function fqwriterecord!(io::IO, outrs::Vector{Vector{UInt8}}, rs::Vector{FqRecord}, stop::Int)
    FqRecord2StringVec!(outrs, rs, stop)
    writebytes(io, outrs)
end

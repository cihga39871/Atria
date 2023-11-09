"""
    bytes_tmp1 = Vector{UInt8}(undef, 67108864) # 2^26

Used for writebytes(io1out::CodecZlibIO, outr1s, range_filter, bytes_tmp1)
"""
bytes_tmp1 = Vector{UInt8}(undef, 67108864) # 2^26
bytes_tmp2 = Vector{UInt8}(undef, 67108864) # 2^26

"""
    write_fqs_threads!(io1out::IO, io2out::IO,
        outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
        r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
        n_reads::Int, range_filter, task_write1, task_write2)

The interface to write paired FASTQ reads. 

- `r1s` and `r2s`: reads to write.
"""
function write_fqs_threads!(io1out::IOStream, io2out::IOStream,
    outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
    r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
    n_reads::Int, range_filter, task_write1, task_write2)

    task_r1s_unbox = Threads.@spawn begin
        wait(task_write1)  # last task
        # @info "write_fqs_threads! FqRecord2StringVec! - start - R1 : n_reads = $n_reads"
        FqRecord2StringVec!(outr1s::Vector{Vector{UInt8}}, r1s::Vector{FqRecord}, n_reads::Int)
        # @info "write_fqs_threads! FqRecord2StringVec! - done  - R1 : n_reads = $n_reads"
    end

    task_write1_new = Threads.@spawn begin
        wait(task_r1s_unbox)  # last task
        writebytes(io1out, outr1s, range_filter) # new task
    end

    task_r2s_unbox = Threads.@spawn begin
        wait(task_write2)
        # @info "write_fqs_threads! FqRecord2StringVec! - start - R2 : n_reads = $n_reads"
        FqRecord2StringVec!(outr2s::Vector{Vector{UInt8}}, r2s::Vector{FqRecord}, n_reads::Int)
        # @info "write_fqs_threads! FqRecord2StringVec! - done  - R2 : n_reads = $n_reads"
    end

    task_write2_new = Threads.@spawn begin
        wait(task_r2s_unbox)  # last task
        writebytes(io2out, outr2s, range_filter)
    end

    task_r1s_unbox, task_r2s_unbox, task_write1_new, task_write2_new
end
function write_fqs_threads!(io1out::IO, io2out::IO,
    outr1s::Vector{Vector{UInt8}}, outr2s::Vector{Vector{UInt8}},
    r1s::Vector{FqRecord}, r2s::Vector{FqRecord},
    n_reads::Int, range_filter, task_write1, task_write2)

    task_r1s_unbox = Threads.@spawn begin
        wait(task_write1)  # last task
        FqRecord2StringVec!(outr1s::Vector{Vector{UInt8}}, r1s::Vector{FqRecord}, n_reads::Int)
    end

    task_write1_new = Threads.@spawn begin
        wait(task_r1s_unbox)  # last task
        writebytes(io1out, outr1s, range_filter, bytes_tmp1) # new task
    end
    
    task_r2s_unbox = Threads.@spawn begin
        wait(task_write2)
        FqRecord2StringVec!(outr2s::Vector{Vector{UInt8}}, r2s::Vector{FqRecord}, n_reads::Int)
    end

    task_write2_new = Threads.@spawn begin
        wait(task_r2s_unbox)  # last task
        writebytes(io2out, outr2s, range_filter, bytes_tmp2)
    end

    task_r1s_unbox, task_r2s_unbox, task_write1_new, task_write2_new
end

function write_fqs_threads!(io1out::IOStream,
    outr1s::Vector{Vector{UInt8}},
    r1s::Vector{FqRecord},
    n_reads::Int, range_filter, task_write1)

    wait(task_write1)
    task_r1s_unbox = Threads.@spawn begin
        FqRecord2StringVec!(outr1s::Vector{Vector{UInt8}}, r1s::Vector{FqRecord}, n_reads::Int)
    end
    task_write1_new = Threads.@spawn begin
        wait(task_r1s_unbox)
        writebytes(io1out, outr1s, range_filter)
    end
    task_r1s_unbox, task_write1_new
end
function write_fqs_threads!(io1out::IO,
    outr1s::Vector{Vector{UInt8}},
    r1s::Vector{FqRecord},
    n_reads::Int, range_filter, task_write1)

    wait(task_write1)
    task_r1s_unbox = Threads.@spawn begin
        FqRecord2StringVec!(outr1s::Vector{Vector{UInt8}}, r1s::Vector{FqRecord}, n_reads::Int)
    end
    task_write1_new = Threads.@spawn begin
        wait(task_r1s_unbox)
        writebytes(io1out, outr1s, range_filter, bytes_tmp1)
    end
    task_r1s_unbox, task_write1_new
end


"""
    FqRecord2StringVec!(out::Vector{UInt8}, r::FqRecord)

Empty `out`, and then convert `r` to it continuously. If empty sequence, write a N as sequence and a ! as quality.
"""
@inline function FqRecord2StringVec!(out::Vector{UInt8}, r::FqRecord)::Nothing
    # out = Base.StringVector(0)
    empty!(out)
    if r.seq.len == 0x0000000000000000  # isempty(r.seq::LongDNA{4})
        append!(out, r.id::Vector{UInt8})
        append!(out, [0x0a, 0x4e, 0x0a])  # \nN\n
        # push!(out, 0x0a)  # \n
        # push!(out, 0x4e)  # N
        # push!(out, 0x0a)
        append!(out, r.des::Vector{UInt8})
        append!(out, [0x0a, 0x21, 0x0a])  # \n!\n
        # push!(out, 0x0a)
        # push!(out, 0x21)  # !
        # push!(out, 0x0a)
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
    if length(rs) < stop
        @error "length(rs) < stop" length(rs) stop
    end
    @sync for reads_start in 1:3072:stop
        reads_end = min(reads_start + 3071, stop)
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
        write_no_lock(io, outrs[i])
    end
    nothing
end
@inline function writebytes(io::IOStream, outrs::Vector{Vector{UInt8}}, filters::SubArray{Bool,1,Array{Bool,1},Tuple{UnitRange{Int64}},true})::Nothing
    @inbounds for (i, val) in enumerate(filters)
        if val
            write_no_lock(io, outrs[i])
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



@inline function writebytes(io::IO, outrs::Vector{Vector{UInt8}}, filters::SubArray{Bool,1,Array{Bool,1},Tuple{UnitRange{Int64}},true}, bytes_tmp::Vector{UInt8})::Nothing
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

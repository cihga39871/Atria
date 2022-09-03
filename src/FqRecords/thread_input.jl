
"""
    load_fqs_threads!(
        io1::IOStream                 , io2::IOStream                 ,
        in1bytes::Vector{UInt8}       , in2bytes::Vector{UInt8}       ,
        vr1s::Vector{Vector{FqRecord}}, vr2s::Vector{Vector{FqRecord}},
        r1s::Vector{FqRecord}         , r2s::Vector{FqRecord}         ;
        remove_first_n::Int64 = 0     , quality_offset::Int64=33      ,
        njobs = 2                     )

The interface to load FASTQ files. Design for `IOStream` only.

Caution: `r[12]s` might be new elements, or in-place modified ones, depending on the necessity to copy.

Return `n_r1, n_r2, r1s, r2s, ncopyed`.
"""
function load_fqs_threads!(
    io1::IOStream                 , io2::IOStream                 ,
    in1bytes::Vector{UInt8}       , in2bytes::Vector{UInt8}       ,
    vr1s::NTuple{N, Vector{FqRecord}}, vr2s::NTuple{N, Vector{FqRecord}},
    # vr1s::Vector{Vector{FqRecord}}, vr2s::Vector{Vector{FqRecord}},
    r1s::Vector{FqRecord}         , r2s::Vector{FqRecord}         ;
    remove_first_n::Int64 = 0     , quality_offset::Int64=33      ,
    njobs = 2                   ) where N

    #== IMPORTANT :: Shared elements conflicts in vrs and rs ::

    Currently, vrs and rs share elements (because of append! in append_loop! in load_fqs_threads!).
    So, in the next run, when modifying vrs, the elements in rs will change.
    It is not a problem if rs are all processed, but when some rs remain,
    it will cause a huge shared data conflict!
    Solution: replace the rs with a new one if rs is not all processed.
    =#
    n_r1 = length(r1s)
    n_r2 = length(r2s)
    ncopyed = 0
    @inbounds if remove_first_n != n_r1  # r1 has some remaining
        if eof(io1)
            # exception: no modification will be done in vr1s, so r1s will not change.
            n_remove_r1 = n_remove_r2 = remove_first_n
        else
            r1s = map(copy, view(r1s, remove_first_n + 1 : n_r1))
            ncopyed = length(r1s)
            # and do not remove in append_loop!
            n_remove_r1 = 0
            n_remove_r2 = remove_first_n
        end
    elseif remove_first_n != n_r2
        if eof(io2)
            n_remove_r1 = n_remove_r2 = remove_first_n
        else
            r2s = map(copy, view(r2s, remove_first_n + 1 : n_r2))
            ncopyed = length(r2s)
            n_remove_r1 = remove_first_n
            n_remove_r2 = 0
        end
    else
        n_remove_r1 = n_remove_r2 = remove_first_n
    end


    ## real thing starts
    njobs < 1 && (njobs = 1)

    idx_fq1_starts, idx_fq1_stops = read_chunks!(io1, in1bytes, njobs)
    task_read2 = Threads.@spawn read_chunks!(io2, in2bytes, njobs)

    tasks = Task[]
    @inbounds for i in eachindex(idx_fq1_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr1s[i], in1bytes, idx_fq1_starts[i], idx_fq1_stops[i]; remove_first_n=-1, quality_offset=quality_offset)
        push!(tasks, i_task)
    end

    idx_fq2_starts, idx_fq2_stops = fetch(task_read2)

    tasks2 = Task[]
    @inbounds for i in eachindex(idx_fq2_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr2s[i], in2bytes, idx_fq2_starts[i], idx_fq2_stops[i]; remove_first_n=-1, quality_offset=quality_offset)
        push!(tasks2, i_task)
    end
    # vr1s_stops = map(fetch, tasks)
    n_r1 = append_loop!(r1s, vr1s, tasks; remove_first_n=n_remove_r1)
    # vr2s_stops = map(fetch, tasks2)
    n_r2 = append_loop!(r2s, vr2s, tasks2; remove_first_n=n_remove_r2)

    (n_r1 == 0 || n_r2 == 0) && throw(error("Numbers of reads in R1 and R2 not the same!"))

    return n_r1, n_r2, r1s, r2s, ncopyed
end

"""
    load_fqs_threads!(
        io1::IOStream                    ,
        in1bytes::Vector{UInt8}          ,
        vr1s::NTuple{N, Vector{FqRecord}},
        r1s::Vector{FqRecord}            ;
        remove_first_n::Int64 = 0        ,
        quality_offset::Int64 = 33       ,
        njobs = 2                        ) where N

The interface to load FASTQ files. Design for `IOStream` only.

Caution: `r[12]s` might be new elements, or in-place modified ones, depending on the necessity to copy.

Return `n_r1, n_r2, r1s, r2s, ncopyed`.
"""
function load_fqs_threads!(
    io1::IOStream                    ,
    in1bytes::Vector{UInt8}          ,
    vr1s::NTuple{N, Vector{FqRecord}},
    r1s::Vector{FqRecord}            ;
    remove_first_n::Int64 = 0        ,
    quality_offset::Int64 = 33       ,
    njobs = 2                        ) where N

    #== IMPORTANT :: Shared elements conflicts in vrs and rs ::

    Currently, vrs and rs share elements (because of append! in append_loop! in load_fqs_threads!).
    So, in the next run, when modifying vrs, the elements in rs will change.
    It is not a problem if rs are all processed, but when some rs remain,
    it will cause a huge shared data conflict!
    Solution: replace the rs with a new one if rs is not all processed.
    =#
    n_r1 = length(r1s)
    ncopyed = 0
    @inbounds if remove_first_n != n_r1  # r1 has some remaining
        if eof(io1)
            # exception: no modification will be done in vr1s, so r1s will not change.
            n_remove_r1 = remove_first_n
        else
            r1s = map(copy, view(r1s, remove_first_n + 1 : n_r1))
            ncopyed = length(r1s)
            # and do not remove in append_loop!
            n_remove_r1 = 0
        end
    else
        n_remove_r1 = remove_first_n
    end


    ## real thing starts
    njobs < 1 && (njobs = 1)

    idx_fq1_starts, idx_fq1_stops = read_chunks!(io1, in1bytes, njobs)

    tasks = Task[]
    @inbounds for i in eachindex(idx_fq1_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr1s[i], in1bytes, idx_fq1_starts[i], idx_fq1_stops[i]; remove_first_n=-1, quality_offset=quality_offset)
        push!(tasks, i_task)
    end

    # vr1s_stops = map(fetch, tasks)
    n_r1 = append_loop!(r1s, vr1s, tasks; remove_first_n=n_remove_r1)

    (n_r1 == 0) && throw(error("Empty input from $io1"))

    return n_r1, r1s, ncopyed
end


"""
    load_fqs_threads!(
        io1::IO                       , io2::IO                       ,
        in1bytes::Vector{UInt8}       , in2bytes::Vector{UInt8}       ,
        in1bytes_nremain::Integer     , in2bytes_nremain::Integer     ,
        vr1s::NTuple{N, Vector{FqRecord}}, vr2s::NTuple{N, Vector{FqRecord}},
        r1s::Vector{FqRecord}         , r2s::Vector{FqRecord}         ;
        will_eof1::Bool = true        , will_eof2::Bool = true        ,
        in1bytes_resize_before_read::Integer = length(in1bytes)       ,
        in2bytes_resize_before_read::Integer = length(in2bytes)       ,
        remove_first_n::Int64 = 0     , njobs = 2                     ,
        quality_offset::Int64 = 33    ) where N

The interface to load FASTQ files. Design for general `IO` that does not support seeking.

- `will_eof1`: if not known, set to true.

Caution: `r[12]s` might be new elements, or in-place modified ones, depending on the necessity to copy.

Return `n_r1, n_r2, r1s, r2s, in1bytes_nremain, in2bytes_nremain, ncopyed`.
"""
function load_fqs_threads!(
    io1::IO                       , io2::IO                       ,
    in1bytes::Vector{UInt8}       , in2bytes::Vector{UInt8}       ,
    in1bytes_nremain::Integer     , in2bytes_nremain::Integer     ,
    vr1s::NTuple{N, Vector{FqRecord}}, vr2s::NTuple{N, Vector{FqRecord}},
    r1s::Vector{FqRecord}         , r2s::Vector{FqRecord}         ;
    will_eof1::Bool = true        , will_eof2::Bool = true        ,
    in1bytes_resize_before_read::Integer = length(in1bytes)       ,
    in2bytes_resize_before_read::Integer = length(in2bytes)       ,
    remove_first_n::Int64 = 0     , njobs = 2                     ,
    quality_offset::Int64 = 33    ) where N

    #== IMPORTANT :: Shared elements conflicts in vrs and rs ::

    Currently, vrs and rs share elements (because of append! in append_loop! in load_fqs_threads!).
    So, in the next run, when modifying vrs, the element in rs will change.
    It is not a problem if rs is all processed, but when some rs remains,
    it will cause a huge shared data conflict!
    Solution: replace the rs with a new one if rs is not all processed.
    =#
    n_r1 = length(r1s)
    n_r2 = length(r2s)
    ncopyed = 0
    @inbounds if remove_first_n != n_r1  # r1 has some remaining
        if eof(io1)
            # exception: no modification will be done in vr1s, so r1s will not change.
            n_remove_r1 = n_remove_r2 = remove_first_n
        else
            r1s = map(copy, view(r1s, remove_first_n + 1 : n_r1))
            ncopyed = length(r1s)
            # and do not remove in append_loop!
            n_remove_r1 = 0
            n_remove_r2 = remove_first_n
        end
    elseif remove_first_n != n_r2
        if eof(io2)
            n_remove_r1 = n_remove_r2 = remove_first_n
        else
            r2s = map(copy, view(r2s, remove_first_n + 1 : n_r2))
            ncopyed = length(r2s)
            n_remove_r1 = remove_first_n
            n_remove_r2 = 0
        end
    else
        n_remove_r1 = n_remove_r2 = remove_first_n
    end

    ## real thing starts
    njobs < 1 && (njobs = 1)

    if nthreads() < 2
        # read r1 and r2 one by one
        idx_fq1_starts, idx_fq1_stops, in1bytes, in1bytes_nremain = read_chunks!(io1, in1bytes, in1bytes_nremain, njobs; will_eof = will_eof1, resize_before_read = in1bytes_resize_before_read)
    else
        # read r1 and r2 at the same time
        task_read1 = Threads.@spawn read_chunks!(io1, in1bytes, in1bytes_nremain, njobs; will_eof = will_eof1, resize_before_read = in1bytes_resize_before_read)
    end

    task_read2 = Threads.@spawn read_chunks!(io2, in2bytes, in2bytes_nremain, njobs, will_eof = will_eof2, resize_before_read = in2bytes_resize_before_read)

    if nthreads() < 2
        # r1 has been done
        nothing
    else
        # wait for r1
        idx_fq1_starts, idx_fq1_stops, in1bytes, in1bytes_nremain = fetch(task_read1)
    end

    tasks = Task[]
    @inbounds for i in eachindex(idx_fq1_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr1s[i], in1bytes, idx_fq1_starts[i], idx_fq1_stops[i]; remove_first_n=-1, quality_offset = quality_offset)
        push!(tasks, i_task)
    end

    idx_fq2_starts, idx_fq2_stops, in2bytes, in2bytes_nremain = fetch(task_read2)

    tasks2 = Task[]
    @inbounds for i in eachindex(idx_fq2_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr2s[i], in2bytes, idx_fq2_starts[i], idx_fq2_stops[i]; remove_first_n=-1, quality_offset = quality_offset)
        push!(tasks2, i_task)
    end
    # vr1s_stops = map(fetch, tasks)
    n_r1 = append_loop!(r1s, vr1s, tasks; remove_first_n=n_remove_r1)
    # vr2s_stops = map(fetch, tasks2)
    n_r2 = append_loop!(r2s, vr2s, tasks2; remove_first_n=n_remove_r2)

    (n_r1 == 0 || n_r2 == 0) && throw(error("Numbers of reads in R1 and R2 not the same!"))

    return n_r1, n_r2, r1s, r2s, in1bytes_nremain, in2bytes_nremain, ncopyed
end


"""
    load_fqs_threads!(
        io1::IO                       ,
        in1bytes::Vector{UInt8}       ,
        in1bytes_nremain::Integer     ,
        vr1s::NTuple{N, Vector{FqRecord}},
        r1s::Vector{FqRecord}         ;
        will_eof1::Bool = true        ,
        remove_first_n::Int64 = 0     ,
        quality_offset::Int64 = 33    ,
        njobs = 2                     ) where N

The interface to load FASTQ files. Design for general `IO` that does not support seeking.

- `will_eof1`: if not known, set to true.

Caution: `r[12]s` might be new elements, or in-place modified ones, depending on the necessity to copy.

Return `n_r1, n_r2, r1s, r2s, in1bytes_nremain, in2bytes_nremain, ncopyed`.
"""
function load_fqs_threads!(
    io1::IO                       ,
    in1bytes::Vector{UInt8}       ,
    in1bytes_nremain::Integer     ,
    vr1s::NTuple{N, Vector{FqRecord}},
    r1s::Vector{FqRecord}         ;
    will_eof1::Bool = true        ,
    remove_first_n::Int64 = 0     ,
    quality_offset::Int64 = 33    ,
    njobs = 2                     ) where N

    #== IMPORTANT :: Shared elements conflicts in vrs and rs ::

    Currently, vrs and rs share elements (because of append! in append_loop! in load_fqs_threads!).
    So, in the next run, when modifying vrs, the element in rs will change.
    It is not a problem if rs is all processed, but when some rs remains,
    it will cause a huge shared data conflict!
    Solution: replace the rs with a new one if rs is not all processed.
    =#
    n_r1 = length(r1s)
    ncopyed = 0
    @inbounds if remove_first_n != n_r1  # r1 has some remaining
        if eof(io1)
            # exception: no modification will be done in vr1s, so r1s will not change.
            n_remove_r1 = remove_first_n
        else
            r1s = map(copy, view(r1s, remove_first_n + 1 : n_r1))
            ncopyed = length(r1s)
            # and do not remove in append_loop!
            n_remove_r1 = 0
        end
    else
        n_remove_r1 = remove_first_n
    end

    ## real thing starts
    njobs < 1 && (njobs = 1)

    idx_fq1_starts, idx_fq1_stops, in1bytes, in1bytes_nremain = read_chunks!(io1, in1bytes, in1bytes_nremain, njobs; will_eof = will_eof1)

    tasks = Task[]
    @inbounds for i in eachindex(idx_fq1_starts)
        i_task = Threads.@spawn StringChunk2FqRecord!(vr1s[i], in1bytes, idx_fq1_starts[i], idx_fq1_stops[i]; remove_first_n=-1, quality_offset = quality_offset)
        push!(tasks, i_task)
    end

    # vr1s_stops = map(fetch, tasks)
    n_r1 = append_loop!(r1s, vr1s, tasks; remove_first_n=n_remove_r1)

    (n_r1 == 0) && throw(error("Empty input from $io1"))

    return n_r1, r1s, in1bytes_nremain, ncopyed
end


"""
    check_filesize(file::String)

If `file` is gziped and can find decompressed file size, return `(decompressed file size, true)`. Otherwise, return `(original file size, false)`.
"""
function check_filesize end

@static if Sys.iswindows()
    @inline function check_filesize(file::String)::Tuple{Int, Bool}
        return (filesize(file), false)
    end
else
    @inline function check_filesize(file::String)::Tuple{Int, Bool}
        # check if file is a symlink
        while islink(file)
            file = joinpath(dirname(file), readlink(file))
        end
        description = try
            readchomp(`file $file`)
        catch
            ""
        end
        # regex to gzip compressed files
        size_match = match(r"original size modulo [^ ]* (\d+)", description)
        if size_match === nothing
            return (filesize(file), false)
        else
            try
                gziped_file_size = parse(Int, size_match.captures[1])
                if gziped_file_size == 0
                    # happens when use non standard compression program, such as fastp
                    return (filesize(file), false)
                else
                    return (gziped_file_size, true)
                end
            catch
                return (filesize(file), false)
            end
        end
    end
end


"""
    chunk_sizes(file1::String, file2::String, max_chunk_size::Int)

Compute how many bytes to read each time (chunk size) based on the file sizes.

Return chunk sizes of file 1 and 2, uncompressed sizes of file 1 and 2. If no compression is detected, or failed to find uncompressed sizes, return -1 for this file.
"""
@inline function chunk_sizes(file1::String, file2::String, max_chunk_size::Int)
    filesize1, isgz1 = check_filesize(file1)
    filesize2, isgz2 = check_filesize(file2)
    if isgz1 != isgz2
        # cannot determine gzip size using `file` for one sample.
        # just use size.
        filesize1 = filesize(file1)
        filesize2 = filesize(file2)
    end
    if filesize1 > filesize2
        chunk_size1 = max_chunk_size
        chunk_size2 = round(Int, max_chunk_size * (filesize2 / filesize1))
    else
        chunk_size1 = round(Int, max_chunk_size * (filesize1 / filesize2))
        chunk_size2 = max_chunk_size
    end
    uncompressed_size1 = isgz1 ? filesize1 : -1
    uncompressed_size2 = isgz2 ? filesize2 : -1
    chunk_size1, chunk_size2, uncompressed_size1, uncompressed_size2
end

"""
    adjust_inbyte_sizes(in1bytes::Vector{UInt8}, in2bytes::Vector{UInt8}, n_r1::Int, n_r2::Int, n_r1_before::Int, n_r2_before::Int, max_chunk_size::Int, default_chunk_size1::Int, default_chunk_size2::Int)

Return the new sizes of `in1bytes` and `in2bytes` according to number of reads in file 1 and file 2. The original sizes are not changed.

Return (chunk_size1, chunk_size2):Tuple{Int,Int}
"""
@inline function adjust_inbyte_sizes(in1bytes::Vector{UInt8}, in2bytes::Vector{UInt8}, n_r1::Int, n_r2::Int, n_r1_before::Int, n_r2_before::Int, max_chunk_size::Int, default_chunk_size1::Int, default_chunk_size2::Int)::Tuple{Int,Int}

    n_r1_read = n_r1 - n_r1_before
    n_r2_read = n_r2 - n_r2_before

    avg_length_r1 = length(in1bytes) / n_r1_read
    avg_length_r2 = length(in2bytes) / n_r2_read

    if isinf(avg_length_r1)
        avg_length_r1 = avg_length_r2 * default_chunk_size1 / default_chunk_size2
    elseif isinf(avg_length_r2)
        avg_length_r2 = avg_length_r1 * default_chunk_size2 / default_chunk_size1
    end

    if 0.975 < n_r1 / n_r2 < 1.025
        # do not care about n_r2 - n_r1
        if avg_length_r1 > avg_length_r2
            chunk_size2 = round(Int, (max_chunk_size/avg_length_r1) * avg_length_r2)
            chunk_size1 = max_chunk_size
        else
            chunk_size1 = round(Int, (max_chunk_size/avg_length_r2) * avg_length_r1)
            chunk_size2 = max_chunk_size
        end
        # resize!(in1bytes, chunk_size1)
        # resize!(in2bytes, chunk_size2)
        return chunk_size1, chunk_size2
    end

    # Try chunk_size2 = max_chunk_size
    chunk_size1 = (n_r2 - n_r1 + max_chunk_size/avg_length_r2) * avg_length_r1
    chunk_size1 = round(Int, chunk_size1)
    if chunk_size1 <= max_chunk_size
        # resize!(in1bytes, chunk_size1)
        # resize!(in2bytes, max_chunk_size)
        return chunk_size1, max_chunk_size
    end

    chunk_size2 = (n_r1 - n_r2 + max_chunk_size/avg_length_r1) * avg_length_r2
    chunk_size2 = round(Int, chunk_size2)
    if chunk_size2 > max_chunk_size
        @warn "Unexpected situation in adjust_inbyte_sizes!: chunk_size2 > max_chunk_size"
    end
    # resize!(in1bytes, max_chunk_size)
    # resize!(in2bytes, chunk_size2)
    return max_chunk_size, chunk_size2
end

"""
    adjust_inbyte_sizes!(in1bytes::Vector{UInt8}, in2bytes::Vector{UInt8}, n_r1::Int, n_r2::Int, n_r1_before::Int, n_r2_before::Int, max_chunk_size::Int, default_chunk_size1::Int, default_chunk_size2::Int)

Resize the sizes of `in1bytes` and `in2bytes` according to number of reads in file 1 and file 2.

Return (chunk_size1, chunk_size2):Tuple{Int,Int}
"""
@inline function adjust_inbyte_sizes!(in1bytes::Vector{UInt8}, in2bytes::Vector{UInt8}, n_r1::Int, n_r2::Int, n_r1_before::Int, n_r2_before::Int, max_chunk_size::Int, default_chunk_size1::Int, default_chunk_size2::Int)::Tuple{Int,Int}
    chunk_size1, chunk_size2 = adjust_inbyte_sizes(in1bytes, in2bytes, n_r1, n_r2, n_r1_before, n_r2_before, max_chunk_size, default_chunk_size1, default_chunk_size2)
    resize!(in1bytes, chunk_size1)
    resize!(in2bytes, chunk_size2)
    return (chunk_size1, chunk_size2)
end

function readfill!(s::IOStream, a::Vector{UInt8})::UInt64
    # read!(s, a) returns nothing or error, ignoring the number of bytes read.
    # the c call is lower wrapper of read!(::IOStream, a), and can return nbytes read.
    nbytes = GC.@preserve a ccall(:ios_readall, Csize_t, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), s, pointer(a), sizeof(a))
end

function Base.read!(s::IO, a::Vector{UInt8}, from::Integer)
    nbytes = length(a) - from + 1
    GC.@preserve a unsafe_read(s, pointer(a) + from - 1, nbytes)
    return a
end
function unsafe_fill!(a::Vector{UInt8}, x::UInt8, from::Integer, to::Integer = length(a))
    pb = pointer(a) + from - 1
    b_length = to - from + 1
    b = unsafe_wrap(Vector{UInt8}, pb, b_length)
    ccall(:memset, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t), b, x, b_length)
    return a
end
function readfill!(io::IO, a::Vector{UInt8}, from::Integer; will_eof::Bool=true)::UInt64
    # If io is not IOStream, things are complicated.
    # IO might not seekable, and read! will return error if EOF.
    # Solution: fill an abnormal character, and determin nbytes read use this.
    eof(io) && return 0x0000000000000000

    will_eof && unsafe_fill!(a, 0x00, from)
    try
        read!(io, a, from)
        nbytes = convert(UInt64, length(a))
    catch
        idx_first_0x00 = find_0x00(a, 1, length(a))
        a[idx_first_0x00] == 0x00 ||
            throw(error("Read to EOF: over-estimated uncompressed file size: $io"))
        nbytes = convert(UInt64, idx_first_0x00 - 1)
    end
end

function readfill!(io::IO, a::Vector{UInt8}; will_eof::Bool=true)::UInt64
    # If io is not IOStream, things are complicated.
    # IO might not seekable, and read! will return error if EOF.
    # Solution: fill an abnormal character, and determin nbytes read use this.
    eof(io) && return 0x0000000000000000

    will_eof && fill!(a, 0x00)
    try
        read!(io, a)
        nbytes = convert(UInt64, length(a))
    catch
        idx_first_0x00 = find_0x00(a, 1, length(a))
        a[idx_first_0x00] == 0x00 ||
            throw(error("Read to EOF: over-estimated uncompressed file size: $io"))
        nbytes = convert(UInt64, idx_first_0x00 - 1)
    end
end

"""
    find_0x00(a::Vector{UInt8}, start::Int, stop::Int)

Used in readfill!(io::IO, a::Vector{UInt8}).
"""
@inline function find_0x00(a::Vector{UInt8}, start::Int, stop::Int)::Int
    @inbounds if stop - start <= 1
        a[start] == 0x00 && return start
        a[stop ] == 0x00 && return stop
    end

    mid = ((stop - start) >> 1) + start
    @inbounds if a[mid] == 0x00
        find_0x00(a, start, mid - 1)
    else
        find_0x00(a, mid + 1, stop)
    end
end

"""
    read_chunks!(io::IOStream, inbytes::Vector{UInt8}, nthread::Int)

It is specialized for `::IOStream`. Because IOStream is seekable, it is possible to change the `io` position to the begining of truncated reads, rather than copying the truncated reads.
"""
function read_chunks!(io::IOStream, inbytes::Vector{UInt8}, nthread::Int)
    nbytes = readfill!(io, inbytes)

    if nbytes == 0x0000000000000000  # read nothing (already end of file before read)
        return UInt64[0x0000000000000001], UInt64[nbytes]
    end
    # Split inbytes into different sections for parallel computing:
    # We set up idx_search_step to get the rough length of section;
    # Then, we use seeklastfq() to retrospect the start of fastq;
    # unit ranges for sections were created.
    # The remaining bytes (idx_fq_starts[end] to the end of inbytes) are stored for the next analysis.
    idx_search_step = nbytes ÷ nthread
    idx_searches = 1:idx_search_step:(nbytes - idx_search_step + 1)

    idx_fq_starts = Vector{UInt}(undef, nthread + 1)
    @inbounds for (i,idx) in enumerate(idx_searches)
        idx_fq_starts[i] = seeklastfq(inbytes, idx)
    end

    @inbounds if nbytes == length(inbytes)
        # read till the end of string vector (inbytes)
        # so a truncated read can be found before nbytes
        idx_fq_starts[end] = seeklastfq(inbytes, nbytes)

        # deal with remaining bytes: change the position of IO
        io_offset = nbytes - idx_fq_starts[end]
        seek(io, position(io) - io_offset - 1)
    else
        # read till the end of file
        @inbounds if inbytes[nbytes] == 0x0a  # the last byte is \n
            # idx_fq_starts[end] should be the start index of the next read,
            # so pretend a read is following.
            idx_fq_starts[end] = nbytes + 1
        else
            idx_fq_starts[end] = nbytes + 2
        end
    end

    # create unit ranges for section
    @inbounds if idx_fq_starts[1] != 1
        throw(error("FASTQ file is invalid: in the begining of file chunk (byte length $(length(inbytes))), the first line does not start with @: $io"))
    end

    # indices to return
    unique!(idx_fq_starts)  # when nbytes is too small, it is useful

    idx_fq_ends = @inbounds idx_fq_starts[2:end] .- 2
    pop!(idx_fq_starts)  # delete the last item

    return (idx_fq_starts, idx_fq_ends)
end

"""
    read_chunks!(io::IO, inbytes::Vector{UInt8}, nremain::Integer, nthread::Int; will_eof::Bool = true)

Return vectors of chunk start indices and stop indices, inbytes (in-place), and number of bytes not processed in inbytes:
(idx_fq_starts, idx_fq_ends, inbytes, new_nremain)
"""
function read_chunks!(io::IO, inbytes::Vector{UInt8}, nremain::Integer, nthread::Int; will_eof::Bool = true, resize_before_read::Integer = length(inbytes))
    length_inbytes = length(inbytes)
    if nremain == length_inbytes
        # do not read IO, just process inbytes.
        nbytes = length_inbytes
    else
        # move the unprocessed part of inbytes (which is in the end) to the front.
        copyto!(inbytes, 1, inbytes, length_inbytes - nremain + 1, nremain)

        if resize_before_read != length_inbytes
            if resize_before_read < nremain
                # do not read anything
                resize!(inbytes, nremain)
                nbytes = length_inbytes
            else
                resize!(inbytes, resize_before_read)
                nbytes = readfill!(io, inbytes, nremain + 1, will_eof=will_eof)
            end
        else
            nbytes = readfill!(io, inbytes, nremain + 1, will_eof=will_eof)
        end
    end

    if nbytes == 0x0000000000000000  # read nothing (already end of file before read)
        return UInt64[0x0000000000000001], UInt64[0x0000000000000000], inbytes, 0
    end
    # Split inbytes into different sections for parallel computing:
    # We set up idx_search_step to get the rough length of section;
    # Then, we use seeklastfq() to retrospect the start of fastq;
    # unit ranges for sections were created.
    # The inbytes bytes (idx_fq_starts[end] to the end of inbytes) are stored for the next analysis.
    idx_search_step = nbytes ÷ nthread
    idx_searches = 1:idx_search_step:(nbytes - idx_search_step + 1)

    idx_fq_starts = Vector{UInt}(undef, nthread + 1)
    @inbounds for (i,idx) in enumerate(idx_searches)
        idx_fq_starts[i] = seeklastfq(inbytes, idx)
    end

    @inbounds if nbytes == length(inbytes)
        # read till the end of string vector (inbytes)
        # so a truncated read can be found before nbytes
        idx_fq_starts[end] = seeklastfq(inbytes, nbytes)

        # deal with unprocessed bytes: compute number of bytes not used
        # new_inbytes_remain = inbytes[idx_fq_starts[end]:end]
        new_nremain = length(inbytes) - idx_fq_starts[end] + 1
    else
        # read till the end of file
        @inbounds if inbytes[nbytes] == 0x0a  # the last byte is \n
            # idx_fq_starts[end] should be the start index of the next read,
            # so pretend a read is following.
            idx_fq_starts[end] = nbytes + 1
        else
            idx_fq_starts[end] = nbytes + 2
        end
        # new_inbytes_remain = Base.StringVector(0)
        new_nremain = 0
    end

    # create unit ranges for section
    if idx_fq_starts[1] != 1
        throw(error("FASTQ file is invalid: in the begining of file chunk (byte length $(length(inbytes))), the first line does not start with @: $io"))
    end

    # indices to return
    unique!(idx_fq_starts)  # when nbytes is too small, it is useful

    idx_fq_ends = @inbounds idx_fq_starts[2:end] .- 2
    pop!(idx_fq_starts)  # delete the last item

    return (idx_fq_starts, idx_fq_ends, inbytes, new_nremain)
end

"""
    read_chunks!(io::IOStream, inbytes::Vector{UInt8})

Designed for `::IOStream`. Not for parallel nor non-seekable `io`.

Return chunk start and stop indices.
"""
function read_chunks!(io::IOStream, inbytes::Vector{UInt8})
    nbytes = readfill!(io, inbytes)

    idx_fq_start = 0x0000000000000001

    if nbytes == 0x0000000000000000  # read nothing (already end of file before read)
        return idx_fq_start, nbytes
    elseif nbytes == length(inbytes)
        # read till the end of string vector (inbytes)
        # so a truncated read can be found before nbytes
        idx_fq_start_next = seeklastfq(inbytes, nbytes)

        idx_fq_stop = idx_fq_start_next - 0x0000000000000002

        # deal with remaining bytes: change the position of IO
        io_offset = nbytes - idx_fq_start_next
        seek(io, position(io) - io_offset - 0x0000000000000001)
    else
        # read till the end of file
        @inbounds if inbytes[nbytes] == 0x0a  # the last byte is \n
            idx_fq_stop = nbytes - 0x0000000000000001
        else
            idx_fq_stop = nbytes
        end
    end
    return (idx_fq_start, idx_fq_stop)
end

"""
    first_byte_of_this_line(inbytes::Vector{UInt8}, i)

Find the first character (UInt8 format) of the the line with the char indexed `i`.

Caution: last index of the previous line can be 0.

Return (last index of the previous line, the value of the first char of this line)
"""
@inline function first_byte_of_this_line(inbytes::Vector{UInt8}, i::T)::Tuple{T,UInt8} where T
    j = i - 1
    while j > 0
        @inbounds if inbytes[j] == 0x0a
            # return (j+1, inbytes[j+1])
            break
        else
            j -= 1
        end
    end
    # j == 0
    return (j, @inbounds inbytes[j+1])
end

@inline function index_next_line_and_copy_this_line(inbytes::Vector{UInt8}, start::UInt64, chunk_end::UInt64)::Tuple{UInt64,Vector{UInt8}}
    stop = start
    while stop <= chunk_end
        @inbounds if inbytes[stop] == 0x0a
            break
        end
        stop += 1
    end
    return (stop + 1, @inbounds inbytes[start:(stop - 1)])
end
@inline function index_next_line_and_n_this_line(inbytes::Vector{UInt8}, start::UInt64, chunk_end::UInt64)::Tuple{UInt64,UInt64}
    stop = start
    while stop <= chunk_end
        @inbounds if inbytes[stop] == 0x0a
            break
        end
        stop += 1
    end
    return (stop + 1, (stop - start))
end
@inline function index_next_line_and_view_this_line(inbytes::Vector{UInt8}, start::UInt64, chunk_end::UInt64)::Tuple{UInt64,SubArray{UInt8,1,Array{UInt8,1},Tuple{UnitRange{Int64}},true}}
    stop = start
    while stop <= chunk_end
        @inbounds if inbytes[stop] == 0x0a
            break
        end
        stop += 1
    end
    return (stop + 1, @inbounds view(inbytes, start:(stop - 1)))
end


"""
    seeklastfq(inbytes::Vector{UInt8}, nbytes::UInt)

Seek the index of the last fastq read.

Caution: This function does not check whether the last fastq read is truncated or complete!
"""
@inline function seeklastfq(inbytes::Vector{UInt8}, nbytes::UInt)::UInt
    # idx_last_fq = nbytes
    i = nbytes
    while i > 0
        i_last_but_1, val_last = first_byte_of_this_line(inbytes, i)

        if val_last == 0x40  # @
            if i_last_but_1 == 0
                # it is the first line start with @. so return the first index
                return 0x0000000000000001
            end
            # val_last_but_1 can be anything
            i_last_but_2, val_last_but_1 = first_byte_of_this_line(inbytes, i_last_but_1)

            #! if bounds error occurs here (i_last_but_2 == 0),
            #! the input file is mot a valid FASTQ!
            #1 but because --check-bounds=no is open,
            #1 we need to throw error, too.
            if i_last_but_2 == 0
                throw(error("FASTQ file is invalid: in the begining of file chunk (byte length $(length(inbytes))), the first two lines are start with (\\n,@) or (+,@) or (@,@)."))
            end
            i_last_but_3, val_last_but_2 = first_byte_of_this_line(inbytes, i_last_but_2)
            if val_last_but_2 == 0x2b  # +
                # Bingo! found (i_last_but_1 + 1) is the the index of the last fastq!
                return i_last_but_1 + 0x0000000000000001
            else
                # last_but_2 shoud be in SEQ line
                # so last_but_3 is ID line!
                i = i_last_but_3
            end
        else
            i = i_last_but_1
        end
    end
    i
end

"""
    append_loop!(rs::Vector{FqRecord}, vrs::Vector{Vector{FqRecord}}, vrs_stops_tasks::Vector{Task}; remove_first_n::Int64=0)

- `vrs_stops_tasks`: it can have less elements than `vrs`. In this way, remaining elements in `vrs` are ignored.

- `remove_first_n`: remove the first NUM elements of `rs` before append.
"""
@inline function append_loop!(rs::Vector{FqRecord}, vrs, vrs_stops_tasks::Vector{Task}; remove_first_n::Int64=0)::Int64

    if remove_first_n != 0
        deleteat!(rs, 1:remove_first_n)
    end

    for (i, task) in enumerate(vrs_stops_tasks)  # cannot enumerate vrs. because vrs_stop might have less elements (especially read till end of file)
        stop = fetch(task)
        append!(rs, @inbounds view(vrs[i], 1:stop))
        # append!(rs, vec[1:vrs_stops[i]])
    end
    length(rs)
end

"""
    StringChunk2FqRecord!(rs::Vector{FqRecord}, inbytes::Vector{UInt8}, idx_fq_start::UInt64, idx_fq_end::UInt64; remove_first_n::Int64=0, quality_offset::Int64=33)

- `rs`: it will be replaced by new FqRecords in place.

- `inbytes`: the chunk of raw bytes. Only analyze `inbytes[idx_fq_start:idx_fq_end]`.

- `remove_first_n`: remove the first N elements of `rs`. `0` = no remove; `-1` = remove all (inline replace); `>0` = deleteat! the first N elements.

Return number of reads in `rs` are valid.
"""
@inline function StringChunk2FqRecord!(rs::Vector{FqRecord}, inbytes::Vector{UInt8}, idx_fq_start::UInt64, idx_fq_end::UInt64; remove_first_n::Int64=0, quality_offset::Int64=33)::Int64

    if remove_first_n == 0
        i_read = n_read = length(rs)
    elseif remove_first_n == -1
        i_read = 0
        n_read = length(rs)
    else
        deleteat!(rs, 1:remove_first_n)
        i_read = n_read = length(rs)
    end

    while idx_fq_start <= idx_fq_end
        i_read += 1
        # @info "test" i_read idx_fq_start
        idx_id_start = idx_fq_start
        idx_seq_start, n_id = index_next_line_and_n_this_line(inbytes, idx_id_start, idx_fq_end)
        idx_des_start, n_seq = index_next_line_and_n_this_line(inbytes, idx_seq_start, idx_fq_end)
        idx_qual_start, n_des = index_next_line_and_n_this_line(inbytes, idx_des_start, idx_fq_end)
        idx_fq_start, n_qual = index_next_line_and_n_this_line(inbytes, idx_qual_start, idx_fq_end)
        if i_read <= n_read
            # TODO: rs_i_read = rs[i_read]
            r = @inbounds rs[i_read]
            safe_copyto!(r.id, inbytes, idx_id_start, n_id)
            safe_copyto!(r.seq, inbytes, idx_seq_start, n_seq)
            bitsafe!(r.seq)
            safe_copyto!(r.des, inbytes, idx_des_start, n_des)
            safe_copyto!(r.qual, inbytes, idx_qual_start, n_qual)
            update_prob_from_qual(r, quality_offset=quality_offset)
        else
            id = @inbounds inbytes[idx_id_start:(idx_id_start + n_id - 1)]
            des = @inbounds inbytes[idx_des_start:(idx_des_start + n_des - 1)]
            qual = @inbounds inbytes[idx_qual_start:(idx_qual_start + n_qual - 1)]

            seq = LongDNA{4}(undef, n_seq)
            safe_copyto!(seq, inbytes, idx_seq_start, n_seq)

            i_fq = FqRecord(id, seq, des, qual, quality_offset=quality_offset)
            push!(rs, i_fq)
        end
    end
    i_read
end

# """
#     StringChunk2FqRecord!(rs::Vector{FqRecord}, inbytes::Vector{UInt8}, idx_fq_start::UInt64, idx_fq_end::UInt64; remove_first_n::Int64=0, quality_offset::Int64=33)
#
# - `rs`: it will be replaced by new FqRecords in place.
#
# - `inbytes`: the chunk of raw bytes. Only analyze `inbytes[idx_fq_start:idx_fq_end]`.
#
# - `remove_first_n`: remove the first N elements of `rs`. `0` = no remove; `-1` = remove all (inline replace); `>0` = deleteat! the first N elements. In reality, those are not removed, but copied from index from `remove_first_n+1`
#
# Return number of reads in `rs` are valid.
# """
# @inline function StringChunk2FqRecord!(rs::Vector{FqRecord}, inbytes::Vector{UInt8}, idx_fq_start::UInt64, idx_fq_end::UInt64; remove_first_n::Int64=0, quality_offset::Int64=33)::Int64
#
#     if remove_first_n == 0
#         i_read = n_read = length(rs)
#     elseif remove_first_n == -1
#         i_read = 0
#         n_read = length(rs)
#     else
#         deleteat!(rs, 1:remove_first_n)
#         i_read = n_read = length(rs)
#
#         n_read = length(rs)
#         i_copy_from = remove_first_n
#         i_copy_to = 0
#         @inbounds while i_copy_from < n_read
#             i_copy_from += 1
#             i_copy_to += 1
#             safe_copyto!(rs[i_copy_to], rs[i_copy_from])
#         end
#         i_read = i_copy_to
#     end
#
#     while idx_fq_start <= idx_fq_end
#         i_read += 1
#         # @info "test" i_read idx_fq_start
#         idx_id_start = idx_fq_start
#         idx_seq_start, n_id = index_next_line_and_n_this_line(inbytes, idx_id_start, idx_fq_end)
#         idx_des_start, n_seq = index_next_line_and_n_this_line(inbytes, idx_seq_start, idx_fq_end)
#         idx_qual_start, n_des = index_next_line_and_n_this_line(inbytes, idx_des_start, idx_fq_end)
#         idx_fq_start, n_qual = index_next_line_and_n_this_line(inbytes, idx_qual_start, idx_fq_end)
#         if i_read <= n_read
#             # TODO: rs_i_read = rs[i_read]
#             r = @inbounds rs[i_read]
#             safe_copyto!(r.id, inbytes, idx_id_start, n_id)
#             safe_copyto!(r.seq, inbytes, idx_seq_start, n_seq)
#             bitsafe!(r.seq)
#             safe_copyto!(r.des, inbytes, idx_des_start, n_des)
#             safe_copyto!(r.qual, inbytes, idx_qual_start, n_qual)
#             update_prob_from_qual(r, quality_offset=quality_offset)
#         else
#             id = @inbounds inbytes[idx_id_start:(idx_id_start + n_id - 1)]
#             des = @inbounds inbytes[idx_des_start:(idx_des_start + n_des - 1)]
#             qual = @inbounds inbytes[idx_qual_start:(idx_qual_start + n_qual - 1)]
#
#             seq = LongDNA{4}(n_seq)
#             safe_copyto!(seq, inbytes, idx_seq_start, n_seq)
#
#             i_fq = FqRecord(id, seq, des, qual, quality_offset=quality_offset)
#             push!(rs, i_fq)
#         end
#     end
#     i_read
# end

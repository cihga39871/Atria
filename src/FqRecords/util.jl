
#=
Some functions, such as BioSequences.throw_encode_error
were modified from BioSequences package developped by BioJulia.
Those functions have their own license:

MIT License

Copyright (c) 2018: BioJulia.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

@inline function check_identifier(r1_id::Vector{UInt8}, r2_id::Vector{UInt8})::Bool
    stop1 = findfirst(x -> x == 0x20 || x == 0x2f, r1_id)  # ' ' or '/'
    if isnothing(stop1)
        return r1_id == r2_id
    end
    stop1 -= 1  # do not count ' ' or '/'
    if length(r2_id) < stop1
        return false
    end
    @inbounds for i in 1:stop1
        if r1_id[i] != r2_id[i]
            return false
        end
    end
    true
end

@inline function check_identifier(r1::FqRecord, r2::FqRecord)::Bool
    check_identifier(r1.id, r2.id)
end

@noinline function throw_identifier_error(r1::FqRecord, r2::FqRecord)
    error("Identifiers of r1 and r2 are not the same!\n   R1: $(String(copy(r1.id)))\n   R2: $(String(copy(r2.id)))")
end

# only modify the error message.
@noinline function BioSequences.throw_encode_error(A::BioSequences.Alphabet, src::AbstractArray{UInt8}, soff::Integer)
    for i in 1:div(64, BioSequences.bits_per_symbol(A))
        sym = src[soff+i-1]
        if BioSequences.stringbyte(A, sym) & 0x80 == 0x80
            # find the context around the error: one previous line and the current line
            nsrc = length(src)
            context_start = soff + i - 2
            context_previous_line = true
            context_end = soff + i
            while context_start > 1
                char = src[context_start]
                if char == 0x0a # \n
                    if context_previous_line
                        context_previous_line = false
                    else
                        context_start += 1
                        break
                    end
                elseif soff - context_start > 300 + 300 * !context_previous_line
                    break
                end
                context_start -= 1
            end
            while context_end < nsrc
                char = src[context_end]
                if char == 0x0a # \n or \r
                    context_end -= 1
                    break
                elseif context_end - soff > 100
                    break
                end
                context_end += 1
            end
            context = String(copy(src[context_start:context_end]))
        
            error("Cannot encode $sym ('$(Char(sym))') to $A. Is the input file valid? Does the disk have bad sections? The error is found in the following context:\n\n$context\n")
        end
    end
end

const complement_table = map(0x01:0xff) do bits
    reinterpret(DNA,
    (bits & 0x01) << 3 | (bits & 0x08) >> 3 |
    (bits & 0x02) << 1 | (bits & 0x04) >> 1)
end

@inline function BioSequences.complement(nt::DNA)  # FASTER than native one
    bits = BioSequences.compatbits(nt)
    @inbounds complement_table[bits]
end

@inline function iscomplement(a::DNA, b::DNA)
    complement(a) === b
end


# codes modified from Julia Base

function write_no_lock(s::IOStream, b::UInt8)
    Int(ccall(:ios_putc, Cint, (Cint, Ptr{Cvoid}), b, s.ios))
end
function write_no_lock(s::IOStream, a::Vector{UInt8})
    GC.@preserve a unsafe_write_no_lock(s, pointer(a), UInt64(sizeof(a)))
end
# """
#     unsafe_write_no_lock(io::IO, ref, nbytes::UInt)
#
# Copy `nbytes` from `ref` (converted to a pointer) into the `IO` object.
#
# It is recommended that subtypes `T<:IO` override the following method signature
# to provide more efficient implementations:
# `unsafe_write_no_lock(s::T, p::Ptr{UInt8}, n::UInt)`
# """
# function unsafe_write_no_lock(s::IO, p::Ptr{UInt8}, n::UInt)
#     written::Int = 0
#     for i = 1:n
#         written += write(s, unsafe_load(p, i))
#     end
#     return written
# end
function unsafe_write_no_lock(s::IOStream, p::Ptr{UInt8}, nb::UInt)
    Int(ccall(:ios_write, Csize_t, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), s.ios, p, nb))
end

# write(io::AbstractPipe, byte::UInt8) = write(Base.pipe_writer(io)::IO, byte)

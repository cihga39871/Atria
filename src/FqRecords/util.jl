
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
        BioSequences.stringbyte(A, sym) & 0x80 == 0x80 && error("Cannot encode $sym ('$(Char(sym))') in $(String(copy(src))) to $A. Is the input file valid? Does the disk have bad sections?")
    end
end

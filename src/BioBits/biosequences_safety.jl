
#=
Some functions, such as BioSequences._orphan!, Base.resize!, and
BioSequences.reverse_complement! were modified from BioSequences package
developped by BioJulia. Those functions have their own license:

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

@inline function isbitsafe(seq::LongDNA{4})
    unsafe_isbitsafe(seq) && 
    seq.data[end] == 0x0000000000000000 &&
    if length(seq.data) > 1 && seq.len % 16 != 0x0000000000000000
        (seq.data[end-1] >> (seq.len % 16 * 4) == 0x0000000000000000)
    else
        true
    end
end

@inline function unsafe_isbitsafe(seq::LongDNA{4})
    length(seq.data) == cld(seq.len, 16) + 1
end

"""
    unsafe_extra_bits_to_zeros!(seq::LongDNA{4})

Caution: use only in bitsafe seq!
"""
@inline function unsafe_extra_bits_to_zeros!(seq::LongDNA{4})
    if !isempty(seq)
        remain = (seq.len % 16)
        @inbounds if remain != 0
            seq.data[end-1] &= ~(0xffffffffffffffff << (remain * 4))
        end
    end
    @inbounds seq.data[end] = 0x0000000000000000
    return seq
end

"""
    bitsafe!(seq::LongDNA{4})

Resize `seq.data` to allow loading a pointer `Ptr{UInt64}` safely at the end of `seq`.

Caution: bitsafe LongDNA{4} may not be compatible on all BioSequences functions, especially those do in-place replacement.
"""
@inline function bitsafe!(seq::LongDNA{4})
    if !unsafe_isbitsafe(seq)
        resize!(seq.data, cld(seq.len, 16) + 1)
    end
    unsafe_extra_bits_to_zeros!(seq)
end

"""
    resize!(seq::LongDNA{4}, size::Int[, force::Bool=false])

It overrides `resize!` in BioSequences. Resize a biological sequence `seq`, to a given `size`. The underlying data is bitsafe.
"""
@inline function Base.resize!(seq::LongSequence{A}, size::Int, force::Bool=false) where {A}
    if size < 0
        throw(ArgumentError("size must be non-negative"))
    else
        if force | (BioSequences.seq_data_len(A, size) > BioSequences.seq_data_len(A, length(seq)))
            resize!(seq.data, BioSequences.seq_data_len(A, size))
        end
        seq.len = size
        bitsafe!(seq)
    end
end

function BioSequences.reverse_complement!(seq::LongSequence{<:NucleicAcidAlphabet})
    pred = x -> BioSequences.complement_bitpar(x, Alphabet(seq))
    BioSequences.reverse_data!(pred, seq.data, BioSequences.seq_data_len(seq) % UInt, BioSequences.BitsPerSymbol(seq))
    BioSequences.zero_offset!(seq)
    bitsafe!(seq)
end

function BioSequences.reverse_complement(seq::LongSequence{<:NucleicAcidAlphabet})
    cp = typeof(seq)(undef, unsigned(length(seq)))
    pred = x -> BioSequences.complement_bitpar(x, Alphabet(seq))
    BioSequences.reverse_data_copy!(pred, cp.data, seq.data, BioSequences.seq_data_len(seq) % UInt, BioSequences.BitsPerSymbol(seq))
    BioSequences.zero_offset!(cp)
    bitsafe!(cp)
end

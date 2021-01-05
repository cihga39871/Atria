
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

"""
Re-write BioSequences._orphan! because original function allocates a vector each time but that is not necessary in Atria.

Caution: it is only tested in Atria, and may not be compatible in all BioSequences functions!

Caution: after using it, r1, r2 and r1_seq_rc, r2_seq_rc can only run `pe_consensus!` once! If run twice, unexpected things can happen.
"""
function BioSequences._orphan!(seq::BioSequences.LongSequence{A},
		 size::Integer = length(seq)) where {A}

    j, r = BioSequences.bitindex(seq, 1)
	new_seq_data_len = BioSequences.seq_data_len(A, size)

    @inbounds if !isempty(seq) && new_seq_data_len > 0
        x = seq.data[j] >> r
        m = BioSequences.index(BioSequences.bitindex(seq, lastindex(seq))) - j + 1
        l = min(new_seq_data_len, m)
        @simd for i in 1:l-1
            y = seq.data[j + i]
            seq.data[i] = x | y << (64 - r)
            x = y >> (r & 63)
        end
        if m <= l
            seq.data[l] = x
        else
            y = seq.data[j + l]
            seq.data[l] = x | y << (64 - r)
        end
    end
	resize!(seq.data, new_seq_data_len)
    seq.part = 1:size
    seq.shared = false
    return seq
end


"""
    bitsafe!(seq::LongDNASeq)

Resize `seq.data` to allow loading a pointer `Ptr{UInt64}` safely at the end of `seq`.

Caution: bitsafe LongDNASeq may not be compatible on all BioSequences functions, especially those do in-place replacement.
"""
@inline function bitsafe!(seq::LongDNASeq)::LongDNASeq
    BioSequences.orphan!(seq)
    seq_data = seq.data
    idx_stop = seq.part.stop
    @inbounds if idx_stop == 0
        resize!(seq_data, 1)
        seq_data[1] = 0x0000000000000000
        return seq
    end
    idx_data_stop = cld(idx_stop, 16)
    length_data = length(seq_data)

    # mask overflow bits to 0
    @inbounds if idx_stop % 16 != 0
        seq_data[idx_data_stop] &= ~(0xffffffffffffffff << (idx_stop % 16 * 4))
    end

    # make `unsafe_bitseq` safe by adding one extra element to `seq.data`
    # FEATURE: `unsafe_bitseq` does not check the boundary,
    # and it might `unsafe_load` the pointer to last base of seq intentionally.
    idx_data_stop_safe = idx_data_stop + 1
    @inbounds if length_data == idx_data_stop_safe
        seq_data[idx_data_stop_safe] = 0x0000000000000000
        return seq
    else
        resize!(seq_data, idx_data_stop_safe)
        seq_data[idx_data_stop_safe] = 0x0000000000000000
        return seq
    end
end


function isbitsafe(seq::LongDNASeq)::Bool
    seq_data = seq.data
    idx_stop = seq.part.stop
    idx_data_stop_safe = cld(idx_stop, 16) + 1
    length_seq_data = length(seq_data)
    length_seq_data >= idx_data_stop_safe || return false
    @inbounds for i = idx_data_stop_safe:length_seq_data
        seq_data[i] == 0x0000000000000000 || return false
    end
    true
end

"""
    resize!(seq::LongDNASeq, size::Int[, force::Bool=false])

It overrides `resize!` in BioSequences. Resize a biological sequence `seq`, to a given `size`. The underlying data is bitsafe.
"""
@inline function Base.resize!(seq::LongDNASeq, size::Int, force::Bool=false)::LongDNASeq
    size < 0 && throw(ArgumentError("resize sequence to a length < 0"))
    size_before = length(seq)
    BioSequences.orphan!(seq, size, force || seq.part.start != 1)
    seq_data = seq.data
    seq.part = 1:size
    @inbounds if size == 0
        resize!(seq_data, 1)
        seq_data[1] = 0x0000000000000000
        return seq
    end
    if size_before >= size
        bitsafe!(seq)
    else
        idx_data_stop_before = cld(size_before, 16)

        # mask original overflow bits to 0
        @inbounds if idx_data_stop_before % 16 != 0
            seq_data[idx_data_stop_before] &= ~(0xffffffffffffffff << (size_before % 16 * 4))
        end

        idx_data_stop_safe = cld(size, 16) + 1
        resize!(seq_data, idx_data_stop_safe)
        @inbounds for i in idx_data_stop_before+1:idx_data_stop_safe
            seq_data[i] = 0x0000000000000000
        end
    end
    seq
end

function BioSequences.reverse_complement!(seq::LongDNASeq)
    # is_seq_bitsafe = isbitsafe(seq)

    bitsafe!(seq)
    seq_data = seq.data
    pred = x -> BioSequences.complement_bitpar(x, DNAAlphabet{4}())
    BioSequences.reverse_data!(pred, seq_data, BioSequences.BitsPerSymbol(seq))
    BioSequences.zero_offset!(seq)

    # if is_seq_bitsafe
        Base.copyto!(seq_data, 1, seq_data, 2, length(seq_data)-1)
        @inbounds seq_data[end] = 0x0000000000000000
    # end
    return seq
end


module BioBits

export bitsafe!,
SeqHead,
SeqHeadSet,
get_pointer,
get_unsafe_index_of_last_bitseq,
unsafe_bitseq,
N2gap,
bin,
bitwise_scan,
_bitwise_scan_fullseq,
bitwise_scan_rc!,
bitwise_scan_rc,
isbitsafe,
insert_size_decision,
insert_size_decision_separate,
is_false_positive
#TruncSeq

# NOTE:TruncSeq has some unknown accuracy problems.

using BioSymbols
using BioSequences

include("biosequences_safety.jl")
include("get_seq.jl")
include("bit_match.jl")
include("insert_size_decision.jl")

end

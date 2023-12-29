
module BioBits

using Reexport

@reexport using BioSymbols
@reexport using BioSequences

# include("biosequences_safety.jl")
# export bitsafe!, 
# isbitsafe

include("get_seq.jl")
export N2gap,
SeqHead,
SeqHeadSet,
get_pointer,
get_unsafe_index_of_last_bitseq,
unsafe_bitseq,
bin

include("bit_match.jl")
export MatchRes, 
bitwise_scan,
_bitwise_scan_fullseq,
bitwise_scan_rc!,
bitwise_scan_rc

include("insert_size_decision.jl")
export insert_size_decision,
insert_size_decision_separate,
is_false_positive,
one_bp_check

end

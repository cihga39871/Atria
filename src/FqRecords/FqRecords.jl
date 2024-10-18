
module FqRecords

export FqRecord,
qualpval,
qualprob,
update_prob_from_qual,
probsum,
probmean,
copyto!,
safe_copyto!,
fqreadrecord,
fqreadrecord!,
fqwriterecord,
check_identifier,
throw_identifier_error,
iscomplement,
load_fqs_threads!,
read_chunks!,
StringChunk2FqRecord!,
chunk_sizes,
get_ideal_inbyte_sizes,
get_ideal_inbyte_sizes!,
write_fqs_threads!,
isinreadlength!,
count_N,
isnotmuchN!,
front_trim!,
tail_trim!,
tail_N_trim!,
tail_low_qual_trim!,
qualitymatch,
seq_complexity,
polyX_tail_scan,
pe_consensus!

using Reexport

@reexport using Base.Threads
using ..BioBits
using ..BioBits.BioSymbols
using ..BioBits.BioSequences

include("interface.jl")
include("quality.jl")
export compute_prob_and_score!

include("copy.jl")
include("basic_io.jl")
include("util.jl")
include("consensus.jl")

include("thread_input.jl")
include("thread_output.jl")

include("check_and_trim.jl")

include("adapter_match_se.jl")
export adapter_match_se

include("adapter_match_pe.jl")
export adapter_match_and_trim_pe!,
adapter_match_pe, 
PEOptions, AdapterPERes

include("pcr_dedup.jl")
export DupCount, pcr_dedup, write_pcr_dedup_count

end

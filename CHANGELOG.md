# Atria Change Log

## v2.1.1

- Fixing wrapper_single_end.jl: cannot trim true adapter position at index of -1.

## v2.1.0

- If a r1/2 adapter is found, but the region of r2/1 is missing or its quality too low (mean prob < 0.6), skip PE check and just trim like single-end. With this, trim_score do not need to compensate for the situation, so rise the default trim-score-pe (10->19).

## v2.0.0

- Supporting low-complexity filtration.

- Supporting polyX tail trimming.

- Supporting single-end fastq.

- Supporting bzip2 compression/decompression.

- Supporting non standardized gzip compression files.

- Optimizing default parameters. (r1-r2-diff 0->0, trim-score-pe 8->10, score-diff removed, kmer-n-match 8->9)

- Robustness optimization: the lower bound of match probability is set to 0.75 because match probability lower than 0.75 is outlier and affect trim score strongly.

## v1.1.1

- Performance optimization: adapter and PE trimming: following v1.1.0-1, if the loosen match's nmatch > trim_score, replace the old one.

## v1.1.0

- Performance optimization: adapter and PE trimming: if no adapters were matched, the number of errors of PE match is loosen.

- Performance optimization: consensus calling: new arg `--kmer-tolerance-consensus 2->10`; optimized arg `--min-ratio-mismatch 0.2->0.28`.

- Speed optimization: check `overlap_score > 0` before computing score (`pe_consensus!`).

## v1.0.3

- More detailed error output when encoding a non-nucleotide character (`throw_encode_error(...)`).

- Following symbolic link before checking file size for non-Windows platforms (`check_filesize(::String)`).

- When run in multi-file parallel mode, write stdout and stderr to a 'stdlog' file (`julia_wrapper_atria(...)`).

- Add option `--check-identifier` to check whether the identifiers of r1 and r2 are the same.

## v1.0.2

- First mature version of Atria.

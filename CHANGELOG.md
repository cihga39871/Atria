# Atria Change Log

## TODO

- Feature: multiple primer trimming.
- Feature: UMI trimming.

## v4.1.2

- Fix: do not throw error if input paired end files are empty when doing `--detect-adapter`.

## v4.1.1

- Fix: `-z NUM -Z NUM` error when length to trim < 0.

## v4.1.0

- Change: --length-range default change from 50:500 to 30:999999.
- Feature: HardClipEnd: new process to hard remove the last N bases.
- Change: names in processing order (--order -O) changed.
- Feature: PCRDedup: remove PCR duplicates from fastq files. The entire paired sequence is compared and hashed. This method require large memory because it stores hashes of reads. To enable, use `--pcr-dedup`.
- Feature: processing stats are recorded in the json file.
- Fix: `polyX_tail_scan` algorithm now is more precise, and tailing Ns also count.



## v4.0.3

- Fix: `--order` or `-O` option should accept multiple arguments.

## v4.0.2

- Fix: `--detect-adapter` for paired reads: refer to index 1 of empty vector when no adapter is found.

## v4.0.1

- Fix: dep cihga39871/BioSequences.jl: detailed error message if input files' line break is '\r\n'.

## v4.0.0

- Optimize: algorithm: now the non-overtrim rate for reads without adapters are higher.
- Feature: re-write trimming to allow trim multiple adapters at the same time. This change is adjusted for metabarcoding data.
- Feature: hard-clip: now hard-clip arguments do differently for r1 and r2. This change is adjusted for metabarcoding data. Remove `-C --clip-after -c --clip5`; add `-b --clip-after-r1 -B --clip-after-r2 -e --clip5-r1 -E --clip5-r2`.
- Optimize: --detect-adapter for paired-end reads now guess adapters from pair information, rather than the existing adapter pool.
- Feature: users can customize order of processing: `-O | --order`.

## v3.2.2-1

- Fix: undef error of is_concensused when enabling --stat (thanks to kalavattam, #3)

## v3.2.2

- Optimize: speed up for threads <= 2.
- Fix: `atria test` should not depend on source files.

## v3.2.1

- Feature: automatically skip completed analyses. Use --force or -f to disable the feature.

## v3.2.0

- Remove multi-proc mode since it is unstable.

## v3.1.4

- Logging: new logging for versions and sample completion.
- Fix v3.1.3: multi-proc mode: Julia v1.8.1 does not allow assign new ARGS, and add `-t nthread` in `julia_args`.
- Fix v3.1.3: pe-consensus: error when `insert_size = -1`; fix trimming when `insert_size = -1`.
- Benchmark `iscomplement` in Atria v3.1.2 and that in BioSequences, and found it is good to stick to BioSequences.

## v3.1.3

- Compatible: Julia v1.8 and BioSequences v3.1.0.
- Fix: quality offset not changed in some places when providing a different --quality-format.
- Fix: use `Base.invokelatest` to bypass world age for functions evaluated at run time.
- Docs: update.

## v3.1.2

- Fix: optimize output file names if ending with .bz2.

## v3.1.1

- Fix: when reporting an encode error, report the previous and current lines instead of the whole chunk of data.

## v3.1.0

- New feature: `--detect-adapter` for adapter determination.

- Fix: when input is an empty compressed fastq, atria exits with error because `read_chunks!(::IO, ...)` should return 4 elements, but returned 2.

## v3.0.3

- Fix v3.0.2: `will_eof` should be true when unknown.

- Do not resize chunk sizes before cycle 1 when inputs are compressed and cannot determine uncompressed sizes. Just assume data are not trimmed before.

## v3.0.2

- Fix uncompressed_size1 not defined on gzipped single-end input (#2).

## v3.0.1

- Avoid to lock `IOStream` when write fastq in thread_output.jl: replace `write(::IOStream, ...)` with `write_no_lock(::IOStream, ...)`. It is slightly faster.

- Speed optimization for consensus calling: overwrite `BioSequences.complement(::DNA)` (1.40X), and define `iscomplement(::DNA, ::DNA)` (1.79X).

- Other minor parallel implementations.

## v3.0.0

- If users choose to trim adapter, check 1 bp offset of adapter sequences. It is because Atria might have 1 bp error in some cases.

## v2.1.2

- Parameter optimization using `atria simulate`: --trim-score-pe 19->10, --tail-length 8->12.

- Development of Atria simulation methods.

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

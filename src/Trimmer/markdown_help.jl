
const atria_markdown_help_text = md"""
# Atria $atria_version

An ultra-fast and accurate adapter and quality trimming software designed for paired-end sequencing data.

If you use Atria, please cite
> Jiacheng Chuan, Aiguo Zhou, Lawrence Richard Hale, Miao He, Xiang Li, Atria: an ultra-fast and accurate trimmer for adapter and quality trimming, Gigabyte, 1, 2021 https://doi.org/10.46471/gigabyte.31

Github: https://github.com/cihga39871/Atria

## Usage

Try `atria -h` or `atria --help` for more information.

### Input and Output

The input files should be paired-end FastQ(.gz|.bz2) files (in the same order), or single-end fastqs:

1. Read 1 files: `-r XXXX_R1.fastq YYYY_R1.fastq.gz ...`

2. Read 2 files (optional): `-R XXXX_R2.fastq YYYY_R2.fastq.gz ...`

Output all files to a directory: `-o PATH` or `--output-dir PATH`. Default is the current directory.

Atria skips completed analysis by default. Use `-f` or `--force` to disable the feature.

### Order of processing

Order of trimming and filtration processing methods. Unlisted process will not be done. See default for process names.

- `--order PROCESS...` or `-O PROCESS...` (default: `POLY-X ADAPTER-CONSENSUS UMI PRIMER CLIP3 CLIP5 QUALITY TAIL-N FILTER-N FILTER-LENGTH FILTER-COMPLEXITY`)

### Trimming methods

#### Poly X Tail Trimming

Remove poly-X tails. Suggest to enable `--polyG` for Illumina NextSeq/NovaSeq data.

- Enable: `--polyG`, `--polyT`, `--polyA`, and/or `--polyC` (default: disabled)

- Trim poly X tail if length > INT: `--poly-length 10`

#### Adapter Trimming

- Read 1 adapter(s) (after DNA insert): `-a SEQ...` or `--adapter1 SEQ...` (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)

- Read 2 adapter(s) (after DNA insert): `-A SEQ...` or `--adapter2 SEQ...` (default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT) (if paired-end)

- Disable: `--no-adapter-trim`

- `--detect-adapter` if you do not know adapter sequences.
   >Atria does not trim detected adapters automatically, please check results first.

#### Paired-end Consensus Calling

The overlapped regions of read pairs are checked and corrected. **It is available only when input files are paired-end and Adapter Trimming is on.**

- Disable: `--no-consensus`

#### UMI (Unique Molecular Identifier)

Trim and extract UMI to the first part of read names, so they can be presented in BAM records after mapping.

- Enable and specify UMI location(s): `--umi LOC...`, and LOC can be:
   + `INDEX1`: the R1 index is UMI.
   + `INDEX2`: the R2 index is UMI.
   + `READ1`: the head of read1 is UMI.
   + `READ2`: the head of read2 is UMI.
   (default: disabled)

- If UMI locations contain `READ1` and/or `READ2`:
   + UMI length argument `--umi-len INT` is required. 
   + Skip several bases after UMI: `--umi-skip INT` (default: 0) 

#### Primer Trimming

Trim primers from 5' and 3' ends (default: no primer trimming)

- Directly provide primer sequence(s):
   + `-m SEQ...` or `--primer1 SEQ...`: primers(s) at 5' end of read 1, and their reverse complement appended to 3' end of read 2.

   + `-M SEQ...` or `--primer1 SEQ...`: primers(s) at 5' end of read 1, and their reverse complement appended to 3' end of read 2.

- Or provide a primer table: `-P FILE` or `--primers FILE`. Format of primer table:
   + Each line is a primer set.
   + Columns are primer1, primer2, primer name.
   + Deliminator is TAB (`\t`).
   + No header line; Lines starts with `#` are ignored.

#### Hard Clip 3' End

Resize reads to a fixed length by discarding extra bases in 3' end (tail).

- Number of bases to keep in read 1: `-b INT` or `--clip-after-r1 INT` (default: disabled)

- Number of bases to keep in read 2: `-B INT` or `--clip-after-r2 INT` (default: disabled)

#### Hard Clip 5' End

Remove the first INT bases from 5' end (front).

- Number of bases to remove in read 1: `-e INT` or `--clip5-r1 INT` (default: disabled)

- Number of bases to remove in read 2: `-E INT` or `--clip5-r2 INT` (default: disabled)

#### Quality Trimming

Trim low-quality tails. Trimming read tails when the average quality of bases in a sliding window is low.

- Average quality threshold: `-q 20` or `--quality-score 20` (default: 20)

- Sliding window length: `--quality-kmer 5` (default: 5)

- FastQ quality format: `--quality-format Illumina1.8`, or `--quality-format 33` (default: 33, ie. Illumina1.8)

- Disable: `--no-quality-trim`

#### Tail N Trimming

Trim N tails.
   - Disable: `--no-tail-n-trim`

#### N Filtration

Discard a read pair if the number of N in one read is greater than a certain amount. N tails are ignored if Tail N Trimming is on.

- Number of N allowed in each read: `-n 15` or `--max-n 15` (default: 15)

- Disable: `-n -1` or `--max-n -1`

#### Read Length Filtration

Filter read pair length in a range.
- Read length range: `--length-range 50:500` (default: 50:500)

- Disable: `--no-length-filtration`

#### Read Complexity Filtration

Discard reads with low complexity. Complexity is the percentage of base that is different from its next base.

    - Enable: `--enable-complexity-filtration` (default: disabled)

    - Complexity threshold: `--min-complexity 0.3` (default: 0.3)

### Parallel (multi-threading) computing

- Use INT threads: `-t 8` or `--threads 8` (default: 8)

- If memory is not sufficient, use `--log2-chunk-size INT` where INT is from 23 to 25. Memory usage reduces exponentially as it decreases.
"""

function atria_markdown_help()
   println(stderr)
   show(stderr, "text/plain", atria_markdown_help_text)
   println(stderr)
end

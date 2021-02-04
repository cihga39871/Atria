
const atria_markdown_help_text = md"""
# Atria

An ultra-fast and accurate adapter and quality trimming software designed for paired-end sequencing data.

Author: Jiacheng (Eric) Chuan

## Usage

Try `atria -h` or `atria --help` for more information.

### Input and Output

The input files should be paired-end FastQ(.gz) files (in the same order), or single-end fastqs:

1. Read 1 files: `-r X_R1.FQ Y_R1.FQ.GZ ...`

2. Read 2 files (optional): `-R X_R2.FQ Y_R2.FQ.GZ ...`

Output all files to a directory: `-o PATH` or `--output-dir PATH`. Default is the current directory.

### Trimming methods

Atria integrated several trimming and read filtration methods. It does the following sequentially.

1. **Poly X Tail Trimming**: remove read tail with poly X.

   suggest to enable `--polyG` for Illumina NextSeq/NovaSeq data.

   - enable: `--polyG`, `--polyT`, `--polyA`, and/or `--polyC` (default: disabled)

   - trim poly X tail if length > INT: `--poly-length 10`

2. **Adapter Trimming**

   - specify read 1 adapter: `-a SEQ` or ` --adapter1 SEQ` (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA)

   - specify read 2 adapter: `-A SEQ` or ` --adapter2 SEQ` (default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT) (if paired-end)

   - disable: `--no-adapter-trim`

3. **Paired-end Consensus Calling**: the overlapped regions of read pairs are checked and corrected. *It is available only when input files are paired-end and Adapter Trimming is on.*
   - disable: `--no-consensus`

4. **Hard Clip 3' end**: resize reads to a fixed length by discarding extra bases in 3' end.
   - specify the number of bases to keep: `-C INT` or `--clip-after INT` (default: disabled)

5. **Hard Clip 5' end**: remove the first INT bases from 5' end.
   - specify the number of bases to remove: `-c INT` or `--clip5 INT` (default: disabled)

6. **Quality Trimming**: trim low-quality tails. (Trimming read tails when the average quality of bases in a sliding window is low.)

   - specify average quality threshold: `-q 20` or `--quality-score 20` (default: 20)

   - specify sliding window length: `--quality-kmer 5` (default: 5)

   - specify FastQ quality format: `--quality-format Illumina1.8`, or `--quality-format 33` (default: 33, ie. Illumina1.8)

   - disable: `--no-quality-trim`

7. **Tail N Trimming**: trim N tails.
   - disable: `--no-tail-n-trim`

8. **N Filtration**: discard a read pair if the number of N in one read is greater than a certain amount. N tails are ignored if Tail N Trimming is on.
   - specify # N allowed in each read: `-n 15` or `--max-n 15` (default: 15)

   - disable: `-n -1` or `--max-n -1`

9. **Read Length Filtration**: filter read pair length in a range.
   - specify read length range: `--length-range 50:500` (default: 50:500)

   - disable: `--no-length-filtration`

10. **Read Complexity Filtration**: filter reads with low complexity.

   Complexity is the percentage of base that is different from its next base.

   - enable: `--enable-complexity-filtration` (default: disabled)

   - specify complexity threshold: `--min-complexity 0.3` (default: 0.3)

### Parallel computing

Atria has two parallel methods: **multi-threading** and **multi-core**. You can enable them both.

1. **Multi-threading** uses multiple threads to process one sample. The speed marginal gain drops as the threads increase: `-t 8` or `--threads 8`. (Default: 8)
2. **Multi-core** processes multiple samples at the same time. No marginal gain loses if IO and memory are sufficient: `-p 4` or `--procs 4`. (Default: disabled)
3. If memory is not sufficient, use `--log2-chunk-size INT` where INT is from 23 to 25. Memory usage reduces exponentially as it decreases.
"""

function atria_markdown_help()
   println(stderr)
   show(stderr, "text/plain", atria_markdown_help_text)
   println(stderr)
end

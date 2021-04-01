#!bash

atria=/home/jiacheng/projects/atria

. $atria/benchmark/trimming-functions.bash

run_all_trimmer() {
    rm -f stderr-simple-time.log
    run_atria 8 2>> stderr-simple-time.log
    run_atria_consensus 8 2>> stderr-simple-time.log
    run_adapterremoval 8 2>> stderr-simple-time.log
    run_skewer 8 2>> stderr-simple-time.log
    run_trim_galore 8 2>> stderr-simple-time.log
    run_trimmomatic 8 2>> stderr-simple-time.log
    run_ktrim 8 2>> stderr-simple-time.log
    pigz -f Ktrim/ktrim.read1.fq Ktrim/ktrim.read2.fq
    run_fastp 8 2>> stderr-simple-time.log
    run_seqpurge 8 2>> stderr-simple-time.log
    run_atropos  8 2>> stderr-simple-time.log
    pasteTimeOutput stderr-simple-time.log > time_benchmark-simple_time.txt
}


####### human data

working_dir=~/analysis/atria-benchmark/ERR4695159
cd $working_dir

r1=ERR4695159_1.fastq.gz
r2=ERR4695159_2.fastq.gz
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
bwa_ref=`pwd`/genomes/hg38.fasta.gz

run_all_trimmer


######## SRR330569: RNA-seq D. simulans

working_dir=~/analysis/atria-benchmark/SRR330569
cd $working_dir

r1=SRR330569.3_1.fastq.gz
r2=SRR330569.3_2.fastq.gz
a1=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
bwa_ref=`pwd`/genomes/dsim-all-chromosome-r2.02.fasta

run_all_trimmer


##### Ensifer spp associated with Medicago whole genome sequencing
working_dir=~/analysis/atria-benchmark/SRR7243169
cd $working_dir

r1=SRR7243169_1.fastq.gz
r2=SRR7243169_2.fastq.gz
a1=CTGTCTCTTATACACATCT
a2=CTGTCTCTTATACACATCT
bwa_ref=`pwd`/genomes/Pseudomonas.sp.Z003-0.4C.fasta

run_all_trimmer

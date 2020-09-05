#!bash
working_dir=~/analysis/atria-benchmark/simulate

cd $working_dir

atria=/home/jiacheng/projects/atria
profile1=`pwd`/ART/art_profiler_illumina/profile_SRR7243169R1.txt
profile2=`pwd`/ART/art_profiler_illumina/profile_SRR7243169R2.txt

replicate=replicate_`date +%s`

mkdir $replicate
cd $replicate


../ART/art_illumina_dir/src/art_illumina --paired --in ../genomes/TAIR10.1.fasta --out TAIR10.sim_ --len 200 --fcov 20 --mflen 250 --sdev 80 -1 $profile1 -2 $profile2

r1=TAIR10.sim_1.fq
r2=TAIR10.sim_2.fq
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# fix a bug of art_illumina: somethimes generating 0x00 in sequences.
# tr -s '\000' 'N' < $r1 > $r1.without_null
# tr -s '\000' 'N' < $r2 > $r2.without_null
#
# mv $r1.without_null $r1
# mv $r2.without_null $r2

# AdapterRemoval only support reading ATCGN, no other ambiguous reads
awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim_1.fq > TAIR10.sim_1.fq.a
awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim_2.fq > TAIR10.sim_2.fq.a

mv $r1.a $r1
mv $r2.a $r2

# get length stats
$atria/benchmark/aln2len.pl TAIR10.sim_1.aln > TAIR10.sim_1.length.tab
$atria/benchmark/aln2len.pl TAIR10.sim_2.aln > TAIR10.sim_2.length.tab

rm TAIR10.sim_1.aln TAIR10.sim_2.aln

. $atria/benchmark/trimming-functions.bash


#### benchmark


rm -f stderr.log


run_atria 1 2>> stderr.log
run_atria 2 2>> stderr.log
run_atria 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab Atria/TAIR10.sim_1.atria.fq Atria/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.Atria
pigz Atria/*



run_atria_consensus 1 2>> stderr.log
run_atria_consensus 2 2>> stderr.log
run_atria_consensus 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab Atria-consensus/TAIR10.sim_1.atria.fq Atria-consensus/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.Atria-consensus
pigz Atria-consensus/*



run_adapterremoval 1 2>> stderr.log
run_adapterremoval 2 2>> stderr.log
run_adapterremoval 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab AdapterRemoval-3/adapterremoval.pair1.truncated AdapterRemoval-3/adapterremoval.pair2.truncated TAIR10.sim_2.length.tab > summary.AdapterRemoval
pigz AdapterRemoval-3/*



run_skewer 1 2>> stderr.log
run_skewer 2 2>> stderr.log
run_skewer 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab Skewer/Skewer-trimmed-pair1.fastq Skewer/Skewer-trimmed-pair2.fastq TAIR10.sim_2.length.tab > summary.Skewer
pigz Skewer/*



run_trim_galore 1 2>> stderr.log
run_trim_galore 2 2>> stderr.log
run_trim_galore 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab TrimGalore/TAIR10.sim_1_val_1.fq TrimGalore/TAIR10.sim_2_val_2.fq TAIR10.sim_2.length.tab > summary.TrimGalore
pigz TrimGalore/*



run_trimmomatic 1 2>> stderr.log
run_trimmomatic 2 2>> stderr.log
run_trimmomatic 4 2>> stderr.log

$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab Trimmomatic/out-pair1.paired.fq Trimmomatic/out-pair2.paired.fq TAIR10.sim_2.length.tab > summary.Trimmomatic
pigz Trimmomatic/*



pasteTimeOutput(){
    paste \
        <(grep -E "Command being timed" $1 | sed 's/.*Command being timed://' ) \
        <(grep -E "^\sUser time" $1) \
        <(grep -E "^\sSystem time" $1) \
        <(grep -E "^\sPercent of CPU this job got" $1) \
        <(grep -E "^\sElapsed" $1) \
        <(grep -E "^\sMaximum resident set size" $1)
}
pasteTimeOutput stderr.log > time_benchmark.txt

#!bash
working_dir=~/analysis/atria-benchmark/art_simulate

cd $working_dir

atria=/home/jiacheng/projects/atria
profile1=`pwd`/ART/art_profiler_illumina/profile_SRR330569R1.txt
profile2=`pwd`/ART/art_profiler_illumina/profile_SRR330569R2.txt

replicate=replicate_`date +%s`

mkdir $replicate
cd $replicate

# # not do it parallel
# parallel -j 15 ../ART/art_illumina_dir/src/art_illumina --paired --in ../genomes/TAIR10.1.fasta --len 100 --fcov 5 --mflen 150 --sdev 80 -1 $profile1 -2 $profile2 --out ::: TAIR10.sim{1..15}_
#
# aln2len(){
#     # aln2len.pl is too slow. this function replaces it.
#     nheader_lines=`grep -n "##Header End" $1 | sed 's/:.*//'`
#     awk "NR>$nheader_lines" $1 | paste - - -  | awk '{sub(/\/.*/, "", $2); print $2"\t"length($5)}' > $1.length.tab
# }
# export -f aln2len
# parallel -j 30 aln2len ::: *aln
#
#
# # get length stats
# cat TAIR10.sim[1-9]*_1.aln.length.tab > TAIR10.sim_1.aln.length.tab
# cat TAIR10.sim[1-9]*_2.aln.length.tab > TAIR10.sim_2.aln.length.tab
#
# # AdapterRemoval only support reading ATCGN, no other ambiguous reads
# awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim[1-9]*_1.fq > TAIR10.sim_1.fq
# awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim[1-9]*_2.fq > TAIR10.sim_2.fq
#
# rm TAIR10.sim[1-9]*

##
../ART/art_illumina_dir/src/art_illumina --paired --in ../genomes/TAIR10.1.fasta --len 100 --fcov 75 --mflen 150 --sdev 80 -1 $profile1 -2 $profile2 --out TAIR10.sim_
aln2len(){
    # aln2len.pl is too slow. this function replaces it.
    nheader_lines=`grep -n "##Header End" $1 | sed 's/:.*//'`
    awk "NR>$nheader_lines" $1 | paste - - -  | awk '{sub(/\/.*/, "", $2); print $2"\t"length($5)}' > $1.length.tab
}
export -f aln2len
parallel -j 2 aln2len ::: TAIR10.sim_1.aln TAIR10.sim_2.aln

# AdapterRemoval only support reading ATCGN, no other ambiguous reads
awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim_1.fq > tmp1
awk '{if (NR%4==2) { gsub("[^ACTGN]", "N") ; print} else {print}}' TAIR10.sim_2.fq > tmp2
mv tmp1 TAIR10.sim_1.fq
mv tmp2 TAIR10.sim_2.fq

# benchmark starts

r1=TAIR10.sim_1.fq
r2=TAIR10.sim_2.fq
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

. $atria/benchmark/trimming-functions.bash


#### benchmark
rm -f stderr.log

run_atria 1 2> stderr.log
run_atria 2 2>> stderr.log
run_atria 4 2>> stderr.log
run_atria 8 2>> stderr.log
run_atria 16 2>> stderr.log
run_atria 32 2>> stderr.log

run_atria_consensus 1 2>> stderr.log
run_atria_consensus 2 2>> stderr.log
run_atria_consensus 4 2>> stderr.log
run_atria_consensus 8 2>> stderr.log
run_atria_consensus 16 2>> stderr.log
run_atria_consensus 32 2>> stderr.log

run_adapterremoval 1 2>> stderr.log
run_adapterremoval 2 2>> stderr.log
run_adapterremoval 4 2>> stderr.log
run_adapterremoval 8 2>> stderr.log
run_adapterremoval 16 2>> stderr.log
run_adapterremoval 32 2>> stderr.log

run_skewer 1 2>> stderr.log
run_skewer 2 2>> stderr.log
run_skewer 4 2>> stderr.log
run_skewer 8 2>> stderr.log
run_skewer 16 2>> stderr.log
run_skewer 32 2>> stderr.log

run_trim_galore 1 2>> stderr.log
run_trim_galore 2 2>> stderr.log
run_trim_galore 4 2>> stderr.log
run_trim_galore 8 2>> stderr.log
run_trim_galore 16 2>> stderr.log
run_trim_galore 32 2>> stderr.log

run_trimmomatic 1 2>> stderr.log
run_trimmomatic 2 2>> stderr.log
run_trimmomatic 4 2>> stderr.log
run_trimmomatic 8 2>> stderr.log
run_trimmomatic 16 2>> stderr.log
run_trimmomatic 32 2>> stderr.log

run_ktrim 1 2>> stderr.log
run_ktrim 2 2>> stderr.log
run_ktrim 4 2>> stderr.log
run_ktrim 8 2>> stderr.log
run_ktrim 16 2>> stderr.log
run_ktrim 32 2>> stderr.log

run_fastp 2 2>> stderr.log
run_fastp 4 2>> stderr.log
run_fastp 8 2>> stderr.log
run_fastp 16 2>> stderr.log
run_fastp 1 2>> stderr.log # multi-thread fastp disorganizes read order, leading to mulfunction of evalTrimming.pl

run_atropos 1 2>> stderr.log
run_atropos 2 2>> stderr.log
run_atropos 4 2>> stderr.log
run_atropos 8 2>> stderr.log
run_atropos 16 2>> stderr.log
run_atropos 32 2>> stderr.log

# SeqPurge cannot output without gzip, so leave it in gzip compression test
# run_seqpurge 1 2>> stderr.log
# run_seqpurge 2 2>> stderr.log
# run_seqpurge 4 2>> stderr.log
# run_seqpurge 8 2>> stderr.log
# run_seqpurge 16 2>> stderr.log
# run_seqpurge 32 2>> stderr.log

pids=

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Atria/TAIR10.sim_1.atria.fq Atria/TAIR10.sim_2.atria.fq TAIR10.sim_2.aln.length.tab > summary.Atria &
pids[1]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Atria-consensus/TAIR10.sim_1.atria.fq Atria-consensus/TAIR10.sim_2.atria.fq TAIR10.sim_2.aln.length.tab > summary.Atria-consensus &
pids[2]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab AdapterRemoval-3/adapterremoval.pair1.truncated AdapterRemoval-3/adapterremoval.pair2.truncated TAIR10.sim_2.aln.length.tab > summary.AdapterRemoval &
pids[3]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Skewer/Skewer-trimmed-pair1.fastq Skewer/Skewer-trimmed-pair2.fastq TAIR10.sim_2.aln.length.tab > summary.Skewer &
pids[4]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab TrimGalore/TAIR10.sim_1_val_1.fq TrimGalore/TAIR10.sim_2_val_2.fq TAIR10.sim_2.aln.length.tab > summary.TrimGalore &
pids[5]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Trimmomatic/out-pair1.paired.fq Trimmomatic/out-pair2.paired.fq TAIR10.sim_2.aln.length.tab > summary.Trimmomatic &
pids[6]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab fastp/out.fastp.r1.fq fastp/out.fastp.r2.fq TAIR10.sim_2.aln.length.tab > summary.Fastp &
pids[7]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Ktrim/ktrim.read1.fq Ktrim/ktrim.read2.fq TAIR10.sim_2.aln.length.tab > summary.Ktrim &
pids[8]=$!

$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Atropos/TAIR10.sim_1.fq.atropos.fq Atropos/TAIR10.sim_2.fq.atropos.fq TAIR10.sim_2.aln.length.tab > summary.Atropos &
pids[9]=$!

# $atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab Atria-src/TAIR10.sim_1.atria.fq Atria-src/TAIR10.sim_2.atria.fq TAIR10.sim_2.aln.length.tab > summary.Atria-src


# $atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab FastqPuri/out.fastqpuri1_good.fq FastqPuri/out.fastqpuri2_good.fq TAIR10.sim_2.aln.length.tab > summary.FastqPuri &
# TP	FP_ft	FP_ot	FN_fr	FN_ut	TN	PPV	Sen.	Spec.	mCC
# 0	185805	666885	11167455	0	32505135	0	0	0.99431631516255	-0.0808480134645243
# (FPR, TPR) = (0.00568368483745041, 0)


for pid in ${pids[*]}
do
    wait $pid
done

##### time of compressed.
pigz $r1 $r2
rm Trimmomatic/*fq
rm TrimGalore/*fq
rm Skewer/*fastq
rm AdapterRemoval-3/*truncated
rm Ktrim/*fq
rm fastp/*fq
rm Atropos/*fq


rm -f stderr.pigz.log

/usr/bin/time pigz -p 8 -c Atria-consensus/*atria.fq 1>/dev/null 2>> stderr.pigz.log
/usr/bin/time pigz -p 8 -c Atria/*atria.fq 1>/dev/null 2>> stderr.pigz.log

rm Atria-consensus/*atria.fq Atria/*atria.fq

r1=$r1.gz
r2=$r2.gz

/usr/bin/time pigz -cd $r1 $r2 > /dev/null 2>> stderr.pigz.log


rm -f stderr.gz.log

run_atria 1 2> stderr.gz.log
run_atria 2 2>> stderr.gz.log
run_atria 4 2>> stderr.gz.log
run_atria 8 2>> stderr.gz.log
run_atria 16 2>> stderr.gz.log
run_atria 32 2>> stderr.gz.log

run_atria_consensus 1 2>> stderr.gz.log
run_atria_consensus 2 2>> stderr.gz.log
run_atria_consensus 4 2>> stderr.gz.log
run_atria_consensus 8 2>> stderr.gz.log
run_atria_consensus 16 2>> stderr.gz.log
run_atria_consensus 32 2>> stderr.gz.log

run_adapterremoval 1 2>> stderr.gz.log
run_adapterremoval 2 2>> stderr.gz.log
run_adapterremoval 4 2>> stderr.gz.log
# run_adapterremoval 8 2>> stderr.gz.log  # too slow
# run_adapterremoval 16 2>> stderr.gz.log

run_skewer 1 2>> stderr.gz.log
run_skewer 2 2>> stderr.gz.log
run_skewer 4 2>> stderr.gz.log
# run_skewer 8 2>> stderr.gz.log  # too slow
# run_skewer 16 2>> stderr.gz.log

run_trim_galore 1 2>> stderr.gz.log
run_trim_galore 2 2>> stderr.gz.log
run_trim_galore 4 2>> stderr.gz.log
run_trim_galore 8 2>> stderr.gz.log
run_trim_galore 16 2>> stderr.gz.log
run_trim_galore 32 2>> stderr.gz.log

run_trimmomatic 1 2>> stderr.gz.log
run_trimmomatic 2 2>> stderr.gz.log
run_trimmomatic 4 2>> stderr.gz.log
# run_trimmomatic 8 2>> stderr.gz.log  # too slow
# run_trimmomatic 16 2>> stderr.gz.log

run_seqpurge 1 2>> stderr.gz.log

pigz -d SeqPurge/TAIR10.sim_1.fq.gz.seqpurge.fq.gz
pigz -d SeqPurge/TAIR10.sim_2.fq.gz.seqpurge.fq.gz
$atria/benchmark/evalTrimming.pl 100 TAIR10.sim_1.aln.length.tab SeqPurge/TAIR10.sim_1.fq.gz.seqpurge.fq SeqPurge/TAIR10.sim_2.fq.gz.seqpurge.fq TAIR10.sim_2.aln.length.tab > summary.SeqPurge

run_seqpurge 2 2>> stderr.gz.log
run_seqpurge 4 2>> stderr.gz.log
run_seqpurge 8 2>> stderr.gz.log
run_seqpurge 16 2>> stderr.gz.log
run_seqpurge 32 2>> stderr.gz.log

run_fastp 4 2>> stderr.gz.log
run_fastp 2 2>> stderr.gz.log
run_fastp 1 2>> stderr.gz.log

# atropos too slow and errors
# run_atropos 4 2>> stderr.gz.log
# run_atropos 2 2>> stderr.gz.log
# run_atropos 1 2>> stderr.gz.log

##### stat
pasteTimeOutput stderr.log > time_benchmark.txt
pasteTimeOutput stderr.gz.log > time_benchmark_gz.txt

NUM_BASES=`grep "read pairs processed" Skewer/Skewer-trimmed.log | awk '{print $1*200}'`
$atria/benchmark/time_stats.jl time_benchmark.txt $NUM_BASES > time_benchmark.df.txt
$atria/benchmark/time_stats.jl time_benchmark_gz.txt $NUM_BASES stderr.pigz.log > time_benchmark_gz.df.txt

$atria/benchmark/time_stats_plot.R -i time_benchmark.df.txt time_benchmark_gz.df.txt


cd ~/analysis/atria-benchmark/julia1.8

a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

r1="reads_diff_indel.R1.fastq"
r2="reads_diff_indel.R2.fastq"

r1="reads_diff_indel.R1.fastq.gz"
r2="reads_diff_indel.R2.fastq.gz"

atria_old=/home/jc/projects/atria/app-3.2.1/bin/atria
atria_new=atria

. $atria/benchmark/trimming-functions.bash

atria simulate --prefix reads_diff_indel --adapter1 $a1 --adapter2 $a2 --repeat 30000 --subsitution-rate 0.001 0.002 0.003 0.004 0.005 --insertion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 --deletion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 -s 100 -i `seq 66 2 120`

NUM_BASES=`echo "4200000 * 200" | bc`


run_atria(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v $atria_old --no-consensus \
        -r $r1 $r1 -R $r2 $r2 \
        -o Atria-old \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads -f
}

run_atria_consensus(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v $atria_old \
        -r $r1 $r1 -R $r2 $r2 \
        -o Atria-consensus-old \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads -f
}

run_atria_new(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v $atria_new --no-consensus \
        -r $r1 $r1 -R $r2 $r2 \
        -o Atria-new \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads -f
}

run_atria_consensus_new(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v $atria_new \
        -r $r1 $r1 -R $r2 $r2 \
        -o Atria-consensus-new \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads -f
}

echo "" 2> stderr.log
run_atria 1 2>> stderr.log
run_atria 2 2>> stderr.log
run_atria 4 2>> stderr.log
run_atria 8 2>> stderr.log
run_atria 16 2>> stderr.log

run_atria_new 1 2>> stderr.log
run_atria_new 2 2>> stderr.log
run_atria_new 4 2>> stderr.log
run_atria_new 8 2>> stderr.log
run_atria_new 16 2>> stderr.log

run_atria_consensus 1 2>> stderr.log
run_atria_consensus 2 2>> stderr.log
run_atria_consensus 4 2>> stderr.log
run_atria_consensus 8 2>> stderr.log
run_atria_consensus 16 2>> stderr.log

run_atria_consensus_new 1 2>> stderr.log
run_atria_consensus_new 2 2>> stderr.log
run_atria_consensus_new 4 2>> stderr.log
run_atria_consensus_new 8 2>> stderr.log
run_atria_consensus_new 16 2>> stderr.log

# run_atria 16 2>> stderr.log
# run_atria_new 16 2>> stderr.log
# run_atria_consensus 16 2>> stderr.log
# run_atria_consensus_new 16 2>> stderr.log

ll */*fastq

for i in *
do
    if [[ -d $i ]]
    then
        atria readstat $i/*.f*q &
    fi
done
wait

atria statplot -i auto -l DIR -F

pasteTimeOutput stderr.log > time_benchmark.txt
$atria/benchmark/time_stats.jl time_benchmark.txt $NUM_BASES > time_benchmark.df.txt

cd ~/analysis/atria-benchmark/julia1.8.5-atria4.0.0

a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

r1="reads_diff_indel.R1.fastq"
r2="reads_diff_indel.R2.fastq"

# r1="reads_diff_indel.R1.fastq.gz"
# r2="reads_diff_indel.R2.fastq.gz"

atria_old=/export/home/CFIA-ACIA/chuanj/software/Atria/app-3.2.1/bin/atria
atria_new=atria

atria_old=/home/jc/projects/Atria/app-3.2.2-1/bin/atria
atria_new=/home/jc/projects/Atria/app-4.0.0-devFixWriteSpeed/bin/atria


. $atria/benchmark/trimming-functions.bash

# atria simulate --prefix reads_diff_indel --adapter1 $a1 --adapter2 $a2 --repeat 300000 --subsitution-rate 0.001 0.002 0.003 0.004 0.005 --insertion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 --deletion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 -s 100 -i `seq 66 2 120`
atria simulate --prefix reads_diff_indel --adapter1 $a1 --adapter2 $a2 --repeat 300000 --subsitution-rate 0.001 0.003 0.005 --insertion-rate 1.0e-5 3.0e-5 5.0e-5 --deletion-rate 1.0e-5 3.0e-5 5.0e-5 -s 100 -i `seq 78 2 108`

NUM_READS=`wc -l reads_diff_indel.R1.fastq | cut -f1 -d" "`
NUM_BASES=`echo "$NUM_READS / 4 * 200" | bc`
echo NUM_BASES=$NUM_BASES


run_atria_src(){
    local num_threads=1
    local outdir=Atria-src
    if [[ $1 ]]; then
        num_threads=$1
    fi
    if [[ $2 ]]; then
        outdir=$2
    fi
	export JULIA_NUM_THREADS=$num_threads
	@time $atria/src/atria --no-consensus -t $num_threads -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 --force
}

run_atria(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    $time -v $atria_old --no-consensus -t $num_threads \
        -r $r1 -R $r2 \
        -o Atria-old \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --force
}

run_atria_consensus(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    $time -v $atria_old -t $num_threads \
        -r $r1 -R $r2 \
        -o Atria-consensus-old \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --force
}

run_atria_new(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    $time -v $atria_new --no-consensus -t $num_threads \
        -r $r1 -R $r2 \
        -o Atria-new-4.0.0-devFixWriteSpeed \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --force
}

run_atria_consensus_new(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    $time -v $atria_new \
        -r $r1 -R $r2 \
        -o Atria-consensus-new-4.0.0-devFixWriteSpeed \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 -t $num_threads --force
}


echo "" 2> stderr.base.log
run_atria 1 2>> stderr.base.log
run_atria 2 2>> stderr.base.log
run_atria 4 2>> stderr.base.log
run_atria 8 2>> stderr.base.log
run_atria 16 2>> stderr.base.log

run_atria_consensus 1 2>> stderr.base.log
run_atria_consensus 2 2>> stderr.base.log
run_atria_consensus 4 2>> stderr.base.log
run_atria_consensus 8 2>> stderr.base.log
run_atria_consensus 16 2>> stderr.base.log

echo "" 2> stderr.dev.log


run_atria_new 1 2>> stderr.dev.log
run_atria_new 2 2>> stderr.dev.log
run_atria_new 4 2>> stderr.dev.log
run_atria_new 8 2>> stderr.dev.log
run_atria_new 16 2>> stderr.dev.log

run_atria_consensus_new 1 2>> stderr.dev.log
run_atria_consensus_new 2 2>> stderr.dev.log
run_atria_consensus_new 4 2>> stderr.dev.log
run_atria_consensus_new 8 2>> stderr.dev.log
run_atria_consensus_new 16 2>> stderr.dev.log

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

cat stderr.base.log stderr.dev.log > std_all.log
pasteTimeOutput std_all.log > time_benchmark.txt
$atria/benchmark/time_stats.jl time_benchmark.txt $NUM_BASES > time_benchmark.df.txt
wps time_benchmark.df.txt & 
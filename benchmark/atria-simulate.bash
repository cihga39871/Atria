#! bash

adapter_length=33

if [[ $1 -lt 33 ]]
then
    adapter_length=$1
fi

working_dir=~/analysis/atria-benchmark/atria_simulate

mkdir -p $working_dir
cd $working_dir

mkdir adapter_length_$adapter_length
cd adapter_length_$adapter_length

# select first adapter_length bp of adapters
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

a1=${a1:0:$adapter_length}
a2=${a2:0:$adapter_length}

#### simulate data with different indels
# docs: https://github.com/cihga39871/Atria/blob/master/docs/3.Benchmark_toolkit.md#data-simulation
atria simulate --prefix reads_diff_indel --adapter1 $a1 --adapter2 $a2 --repeat 30000 --subsitution-rate 0.001 0.002 0.003 0.004 0.005 --insertion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 --deletion-rate 1.0e-5 2.0e-5 3.0e-5 4.0e-5 5.0e-5 -s 100 -i `seq 66 2 120`

r1="reads_diff_indel.R1.fastq"
r2="reads_diff_indel.R2.fastq"

# load trimming functions
. $atria/benchmark/trimming-functions.bash

rm -f stderr.log
run_atria 16 2>> stderr.log
run_adapterremoval 16 2>> stderr.log
run_skewer 16 2>> stderr.log
run_trim_galore 16 2>> stderr.log
run_trimmomatic 16 2>> stderr.log
# run_ktrim 16 2>> stderr.log # ktrim fails to output validate fastq
run_fastp 1 2>> stderr.log
run_atropos 16 2>> stderr.log
run_seqpurge 1 2>> stderr.log
run_cutadapt 16 2>> stderr.log

pigz -d SeqPurge/*gz

mv AdapterRemoval-3/adapterremoval.pair1.truncated AdapterRemoval-3/adapterremoval.pair1.fq
mv AdapterRemoval-3/adapterremoval.pair2.truncated AdapterRemoval-3/adapterremoval.pair2.fq

# cat Trimmomatic/out-pair1.unpaired.fq >> Trimmomatic/out-pair1.paired.fq
# cat Trimmomatic/out-pair2.unpaired.fq >> Trimmomatic/out-pair2.paired.fq
# rm Trimmomatic/out-pair1.unpaired.fq Trimmomatic/out-pair2.unpaired.fq

ll */*fastq */*fq

for i in *
do
    if [[ -d $i ]]
    then
        julia -L $atria/src/Atria.jl -e "Atria.Benchmark.julia_wrapper_readstat(ARGS)" $i/*.f*q &
    fi
done

# atria readstat Cutadapt/out.cutadapt.R*.fq

ps -x | grep -c "Atria.Benchmark.julia_wrapper_readstat"


### Adapter length 16, 20, 24, 28, 33
# atria statplot -i */*r12.stat.tsv

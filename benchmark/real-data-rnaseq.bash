
##### SRR330569: RNA-seq D. simulans
working_dir=~/analysis/atria-benchmark/SRR330569
cd $working_dir

r1=SRR330569.3_1.fastq.gz
r2=SRR330569.3_2.fastq.gz
a1=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
bwa_ref=genomes/Drosophila.simulans.fasta

# download reference
mkdir genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/754/195/GCF_000754195.2_ASM75419v2/GCF_000754195.2_ASM75419v2_genomic.fna.gz -O $bwa_ref.gz
gzip -d $bwa_ref.gz

# build reference
hisat2-build $bwa_ref $bwa_ref-hisat2


##### Pipelines

. $atria/benchmark/trimming-functions.bash

rm -f stderr.log
run_atria 4 2>> stderr.log

run_atria_consensus 4 2>> stderr.log

run_adapterremoval 4 2>> stderr.log

run_skewer 4 2>> stderr.log

run_trim_galore 4 2>> stderr.log

run_trimmomatic 4 2>> stderr.log


# quality trimming
QSCORE=15
qualtrim Atria/${r1/.fastq.gz}.atria.fastq.gz Atria/${r2/.fastq.gz}.atria.fastq.gz
qualtrim Atria-consensus/${r1/.fastq.gz}.atria.fastq.gz Atria-consensus/${r2/.fastq.gz}.atria.fastq.gz
qualtrim AdapterRemoval-3/adapterremoval.pair1.truncated.gz AdapterRemoval-3/adapterremoval.pair2.truncated.gz
qualtrim Skewer/Skewer-trimmed-pair1.fastq.gz Skewer/Skewer-trimmed-pair2.fastq.gz
qualtrim TrimGalore/${r1/.fastq.gz}_val_1.fq.gz TrimGalore/${r2/.fastq.gz}_val_2.fq.gz
qualtrim Trimmomatic/out-pair1.paired.fq.gz Trimmomatic/out-pair2.paired.fq.gz

for d in */qualtrim
do
	echo $d
	if [[ `ls $d/*.qual$QSCORE.truncated.gz | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.truncated.gz
		mapping_hisat2 $d/*.qual$QSCORE.truncated.gz
	elif [[ `ls $d/*.qual$QSCORE.fq.gz | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fq.gz
		mapping_hisat2 $d/*.qual$QSCORE.fq.gz
	elif [[ `ls $d/*.qual$QSCORE.fastq.gz | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fastq.gz
		mapping_hisat2 $d/*.qual$QSCORE.fastq.gz
	fi
done



for i in */qualtrim/*[sb]am
do
    echo `date` - $i
    samtools stats $i > $i.samtools-stats
done

# pasteSamtoolsStats */qualtrim/*bam.samtools-stats
pasteSamtoolsStats */qualtrim/*.samtools-stats


##### Ensifer spp associated with Medicago whole genome sequencing
# Miseq 301bp paired reads
working_dir=~/analysis/atria-benchmark/SRR7243169
cd $working_dir

r1=SRR7243169_1.fastq
r2=SRR7243169_2.fastq
a1=CTGTCTCTTATACACATCT
a2=CTGTCTCTTATACACATCT
bwa_ref=genomes/Pseudomonas.sp.Z003-0.4C.fasta

# download reference
mkdir genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/104/975/GCF_900104975.1_IMG-taxon_2675902959_annotated_assembly/GCF_900104975.1_IMG-taxon_2675902959_annotated_assembly_genomic.fna.gz -O $bwa_ref.gz
gzip -d $bwa_ref.gz

# build reference
bwa index $bwa_ref
bowtie2-build -f $bwa_ref $bwa_ref-bowtie2


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
qualtrim Atria/${r1/.fastq}.atria.fastq Atria/${r2/.fastq}.atria.fastq
qualtrim Atria-consensus/${r1/.fastq}.atria.fastq Atria-consensus/${r2/.fastq}.atria.fastq
qualtrim AdapterRemoval-3/adapterremoval.pair1.truncated AdapterRemoval-3/adapterremoval.pair2.truncated
qualtrim Skewer/Skewer-trimmed-pair1.fastq Skewer/Skewer-trimmed-pair2.fastq
qualtrim TrimGalore/${r1/.fastq}_val_1.fq TrimGalore/${r2/.fastq}_val_2.fq
qualtrim Trimmomatic/out-pair1.paired.fq Trimmomatic/out-pair2.paired.fq

for d in */qualtrim
do
	echo $d
	if [[ `ls $d/*.qual$QSCORE.truncated | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.truncated
		mapping_bowtie2 $d/*.qual$QSCORE.truncated
	elif [[ `ls $d/*.qual$QSCORE.fq | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fq
		mapping_bowtie2 $d/*.qual$QSCORE.fq
	elif [[ `ls $d/*.qual$QSCORE.fastq | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fastq
		mapping_bowtie2 $d/*.qual$QSCORE.fastq
	fi
done

bowtie2stat 15



for d in */qualtrim
do
	echo $d
	if [[ `ls $d/*.qual$QSCORE.truncated | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.truncated
		mapping $d/*.qual$QSCORE.truncated
	elif [[ `ls $d/*.qual$QSCORE.fq | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fq
		mapping $d/*.qual$QSCORE.fq
	elif [[ `ls $d/*.qual$QSCORE.fastq | wc -l` == 2 ]]
	then
		ll $d/*.qual$QSCORE.fastq
		mapping $d/*.qual$QSCORE.fastq
	fi
done

for i in */qualtrim/*bam
do
    echo `date` - $i
    samtools stats $i > $i.samtools-stats
done

# pasteSamtoolsStats */qualtrim/*bam.samtools-stats
pasteSamtoolsStats */qualtrim/*{d,q}.bam.samtools-stats */qualtrim/*bowtie2.bam.samtools-stats


samalign(){
	samtools view $1 | awk '{print $1"\t"$3"\t"$4"\t"$6}'
}

sampaste(){
	paste <(samtools view $1 | awk '{print $1"\t"$3"\t"$4"\t"$6}') <(samtools view $2 | awk '{print $1"\t"$3"\t"$4"\t"$6}') | less -SN
}

akadiff(){
	diff <(samalign Atria-src/qualtrim/SRR7243169_1.removeN.atria.qual15.fastq.bowtie2.sort.bam) <(samalign Skewer/qualtrim/Skewer-trimmed-pair1.qual15.fastq.bowtie2.sort.bam) -y | less
}

# samsort AdapterRemoval-3/qualtrim/adapterremoval.pair1.qual15.truncated.bowtie2.bam
# samsort Atria-src/qualtrim/SRR7243169_1.removeN.atria.qual15.fastq.bowtie2.bam
# samsort Skewer/qualtrim/Skewer-trimmed-pair1.qual15.fastq.bowtie2.bam

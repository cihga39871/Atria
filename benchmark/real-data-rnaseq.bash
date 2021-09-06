
##### SRR330569: RNA-seq D. simulans
working_dir=~/analysis/atria-benchmark/SRR330569
cd $working_dir

r1=SRR330569.3_1.fastq.gz
r2=SRR330569.3_2.fastq.gz
a1=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
bwa_ref=`pwd`/genomes/dsim-all-chromosome-r2.02.fasta

# download reference
mkdir genomes
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/754/195/GCF_000754195.2_ASM75419v2/GCF_000754195.2_ASM75419v2_genomic.fna.gz -O $bwa_ref_genbank.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.02_FB2020_03/fasta/dsim-all-chromosome-r2.02.fasta.gz -O $bwa_ref.gz

gzip -d $bwa_ref.gz

# build reference
hisat2-build $bwa_ref $bwa_ref-hisat2


##### Pipelines

. $atria/benchmark/trimming-functions.bash

rm -f stderr.log
run_atria 8 2>> stderr.log

run_atria_consensus 8 2>> stderr.log

run_adapterremoval 8 2>> stderr.log

run_skewer 8 2>> stderr.log

run_trim_galore 8 2>> stderr.log

run_trimmomatic 8 2>> stderr.log

run_ktrim 8 2>> stderr.log
pigz Ktrim/ktrim.read1.fq Ktrim/ktrim.read2.fq

run_fastp 8 2>> stderr.log

run_seqpurge 8 2>> stderr.log

run_atropos  8 2>> stderr.log

# mapping without qualtrim
mkdir -p trimmed
ln -s ../AdapterRemoval-3/adapterremoval.pair1.truncated.gz trimmed/adapterremoval.R1.fastq.gz
ln -s ../AdapterRemoval-3/adapterremoval.pair2.truncated.gz trimmed/adapterremoval.R2.fastq.gz

ln -s ../Atria/${r1/.fastq*/}.atria.fastq.gz trimmed/atria.R1.fastq.gz
ln -s ../Atria/${r2/.fastq*/}.atria.fastq.gz trimmed/atria.R2.fastq.gz

ln -s ../Atria-consensus/${r1/.fastq*/}.atria.fastq.gz trimmed/atria-consensus.R1.fastq.gz
ln -s ../Atria-consensus/${r2/.fastq*/}.atria.fastq.gz trimmed/atria-consensus.R2.fastq.gz

ln -s ../Skewer/Skewer-trimmed-pair1.fastq.gz trimmed/Skewer.R1.fastq.gz
ln -s ../Skewer/Skewer-trimmed-pair2.fastq.gz trimmed/Skewer.R2.fastq.gz

ln -s ../TrimGalore/${r1/.fastq*/}_val_1.fq.gz trimmed/trimgalore.R1.fastq.gz
ln -s ../TrimGalore/${r2/.fastq*/}_val_2.fq.gz trimmed/trimgalore.R2.fastq.gz

ln -s ../Trimmomatic/out-pair1.paired.fq.gz trimmed/trimmomatic.R1.fastq.gz
ln -s ../Trimmomatic/out-pair2.paired.fq.gz trimmed/trimmomatic.R2.fastq.gz

ln -s ../Ktrim/ktrim.read1.fq.gz trimmed/ktrim.R1.fastq.gz
ln -s ../Ktrim/ktrim.read2.fq.gz trimmed/ktrim.R2.fastq.gz

ln -s ../fastp/out.fastp.r1.fq.gz trimmed/fastp.R1.fastq.gz
ln -s ../fastp/out.fastp.r2.fq.gz trimmed/fastp.R2.fastq.gz

ln -s ../SeqPurge/SRR330569.3_1.fastq.gz.seqpurge.fq.gz trimmed/seqpurge.R1.fastq.gz
ln -s ../SeqPurge/SRR330569.3_2.fastq.gz.seqpurge.fq.gz trimmed/seqpurge.R2.fastq.gz

ln -s ../Atropos/SRR330569.3_1.fastq.gz.atropos.fq.gz trimmed/atropos.R1.fastq.gz
ln -s ../Atropos/SRR330569.3_2.fastq.gz.atropos.fq.gz trimmed/atropos.R2.fastq.gz


# mapping after qualtrim
QSCORE=15
time atria -r trimmed/*.R1.fastq.gz -R trimmed/*.R2.fastq.gz -t 5 -p 6 -o trimmed-qualtrim --no-tail-n-trim --max-n=-1 --no-adapter-trim --no-length-filtration --quality-score $QSCORE --check-identifier
rename --force "s/atria.fastq/qual$QSCORE.fastq/" trimmed-qualtrim/*fastq*

# time atria -r trimmed/ktrim.R1.fastq.gz trimmed/fastp.R1.fastq.gz trimmed/seqpurge.R1.fastq.gz trimmed/atropos.R1.fastq.gz -R trimmed/ktrim.R2.fastq.gz trimmed/fastp.R2.fastq.gz trimmed/seqpurge.R2.fastq.gz trimmed/atropos.R2.fastq.gz -t 7 -p 4 -o trimmed-qualtrim --no-tail-n-trim --max-n=-1 --no-adapter-trim --no-length-filtration --quality-score $QSCORE --check-identifier
# rename --force "s/atria.fastq/qual$QSCORE.fastq/" trimmed-qualtrim/*fastq*

# trimmed*/{atropos,fastp,ktrim,seqpurge}.R1*fastq.gz
for i in trimmed*/*.R1*fastq.gz
do
	echo $i
	mapping_hisat2 $i ${i/.R1/.R2}

	samtools stats $i.hisat2.sam > $i.hisat2.sam.samtools-stats &
	pigz $i.hisat2.sam
done 2>&1 | tee mapping.log


cd trimmed
pasteSamtoolsStats *samtools-stats
cd ..


cd trimmed-qualtrim
pasteSamtoolsStats *samtools-stats
cd ..

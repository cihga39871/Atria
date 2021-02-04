
##### SRR330569: RNA-seq D. simulans
working_dir=~/analysis/atria-benchmark/ERR4695159
cd $working_dir

r1=ERR4695159_1.fastq.gz
r2=ERR4695159_2.fastq.gz
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
bwa_ref=`pwd`/genomes/hg38.fasta.gz

# download reference
mkdir genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -O $bwa_ref

gzip -d $bwa_ref.gz

# build reference
bowtie2-build $bwa_ref $bwa_ref-bowtie2


##### Pipelines

. $atria/benchmark/trimming-functions.bash

rm -f stderr.log
run_atria 16 2>> stderr.log

run_atria_consensus 16 2>> stderr.log

run_adapterremoval 16 2>> stderr.log

run_skewer 16 2>> stderr.log

run_trim_galore 16 2>> stderr.log

run_trimmomatic 16 2>> stderr.log

run_ktrim 16 2>> stderr.log
pigz Ktrim/ktrim.read1.fq Ktrim/ktrim.read2.fq

run_fastp 16 2>> stderr.log

run_seqpurge 16 2>> stderr.log

run_atropos  16 2>> stderr.log

# mapping without qualtrim
mkdir -p trimmed
ln -s ../AdapterRemoval-3/adapterremoval.pair1.truncated.gz trimmed/adapterremoval.R1.fastq.gz
ln -s ../AdapterRemoval-3/adapterremoval.pair2.truncated.gz trimmed/adapterremoval.R2.fastq.gz

ln -s ../Atria/${r1/.fastq*/}.atria.fastq.gz trimmed/atria.R1.fastq.gz
ln -s ../Atria/${r2/.fastq*/}.atria.fastq.gz trimmed/atria.R2.fastq.gz

ln -s ../Atria-consensus/${r1/.fastq*/}.atria.fastq.gz trimmed/atria-consensus.R1.fastq.gz
ln -s ../Atria-consensus/${r2/.fastq*/}.atria.fastq.gz trimmed/atria-consensus.R2.fastq.gz

# ln -s ../Atria-src/${r1/.fastq*/}.atria.fastq.gz trimmed/atria-src.R1.fastq.gz
# ln -s ../Atria-src/${r2/.fastq*/}.atria.fastq.gz trimmed/atria-src.R2.fastq.gz
# ln -s ../Atria-consensus-src/${r1/.fastq*/}.atria.fastq.gz trimmed/atria-consensus-src.R1.fastq.gz
# ln -s ../Atria-consensus-src/${r2/.fastq*/}.atria.fastq.gz trimmed/atria-consensus-src.R2.fastq.gz

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

ln -s ../SeqPurge/ERR4695159_1.fastq.gz.seqpurge.fq.gz trimmed/seqpurge.R1.fastq.gz
ln -s ../SeqPurge/ERR4695159_2.fastq.gz.seqpurge.fq.gz trimmed/seqpurge.R2.fastq.gz

ln -s ../Atropos/ERR4695159_1.fastq.gz.atropos.fq.gz trimmed/atropos.R1.fastq.gz
ln -s ../Atropos/ERR4695159_2.fastq.gz.atropos.fq.gz trimmed/atropos.R2.fastq.gz


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
	mapping_bowtie2 $i ${i/.R1/.R2}

	samtools stats $i.bowtie2.bam > $i.bowtie2.bam.samtools-stats &
done 2>&1 | tee mapping.log


cd trimmed
pasteSamtoolsStats *samtools-stats
cd ..


cd trimmed-qualtrim
pasteSamtoolsStats *samtools-stats
cd ..

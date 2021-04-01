
##### Ensifer spp associated with Medicago whole genome sequencing
# Miseq 301bp paired reads
working_dir=~/analysis/atria-benchmark/SRR7243169
cd $working_dir

r1=SRR7243169_1.fastq
r2=SRR7243169_2.fastq
a1=CTGTCTCTTATACACATCT
a2=CTGTCTCTTATACACATCT
bwa_ref=`pwd`/genomes/Pseudomonas.sp.Z003-0.4C.fasta

# download reference
mkdir genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/104/975/GCF_900104975.1_IMG-taxon_2675902959_annotated_assembly/GCF_900104975.1_IMG-taxon_2675902959_annotated_assembly_genomic.fna -O $bwa_ref
gzip -d $bwa_ref

# build reference
bwa index $bwa_ref
# bowtie2-build -f $bwa_ref $bwa_ref-bowtie2


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

run_fastp 8 2>> stderr.log

run_seqpurge 8 2>> stderr.log
pigz -d SeqPurge/*gz

run_atropos 8 2>> stderr.log


# mapping without qualtrim
mkdir -p trimmed
ln -s ../AdapterRemoval-3/adapterremoval.pair1.truncated trimmed/adapterremoval.R1.fastq
ln -s ../AdapterRemoval-3/adapterremoval.pair2.truncated trimmed/adapterremoval.R2.fastq

ln -s ../Atria/${r1/.fastq*/}.atria.fastq trimmed/atria.R1.fastq
ln -s ../Atria/${r2/.fastq*/}.atria.fastq trimmed/atria.R2.fastq

ln -s ../Atria-consensus/${r1/.fastq*/}.atria.fastq trimmed/atria-consensus.R1.fastq
ln -s ../Atria-consensus/${r2/.fastq*/}.atria.fastq trimmed/atria-consensus.R2.fastq

ln -s ../Skewer/Skewer-trimmed-pair1.fastq trimmed/Skewer.R1.fastq
ln -s ../Skewer/Skewer-trimmed-pair2.fastq trimmed/Skewer.R2.fastq

ln -s ../TrimGalore/${r1/.fastq*/}_val_1.fq trimmed/trimgalore.R1.fastq
ln -s ../TrimGalore/${r2/.fastq*/}_val_2.fq trimmed/trimgalore.R2.fastq

ln -s ../Trimmomatic/out-pair1.paired.fq trimmed/trimmomatic.R1.fastq
ln -s ../Trimmomatic/out-pair2.paired.fq trimmed/trimmomatic.R2.fastq

ln -s ../Ktrim/ktrim.read1.fq trimmed/ktrim.R1.fastq
ln -s ../Ktrim/ktrim.read2.fq trimmed/ktrim.R2.fastq

ln -s ../fastp/out.fastp.r1.fq trimmed/fastp.R1.fastq
ln -s ../fastp/out.fastp.r2.fq trimmed/fastp.R2.fastq

ln -sf ../SeqPurge/SRR7243169_1.fastq.seqpurge.fq trimmed/seqpurge.R1.fastq
ln -sf ../SeqPurge/SRR7243169_2.fastq.seqpurge.fq trimmed/seqpurge.R2.fastq

ln -sf ../Atropos/SRR7243169_1.fastq.atropos.fq trimmed/atropos.R1.fastq
ln -sf ../Atropos/SRR7243169_2.fastq.atropos.fq trimmed/atropos.R2.fastq

# trimmed/{atropos,fastp,ktrim,seqpurge}.R1*fastq
for i in trimmed/*.R1.fastq
do
	echo $i
	mapping $i ${i/.R1./.R2.}
	samtools stats $i.bam > $i.bam.samtools-stats
done 2>&1 | tee mapping.log

cd trimmed
pasteSamtoolsStats *samtools-stats
cd ..

# mapping after qualtrim
QSCORE=15
time atria -r trimmed/*.R1.fastq -R trimmed/*.R2.fastq -t 5 -p 6 -o trimmed-qualtrim --no-tail-n-trim --max-n=-1 --no-adapter-trim --no-length-filtration --quality-score $QSCORE
rename --force "s/atria.fastq/qual$QSCORE.fastq/" trimmed-qualtrim/*fastq

# time atria -r trimmed/ktrim.R1.fastq trimmed/fastp.R1.fastq trimmed/seqpurge.R1.fastq trimmed/atropos.R1.fastq -R trimmed/ktrim.R2.fastq trimmed/fastp.R2.fastq trimmed/seqpurge.R2.fastq trimmed/atropos.R2.fastq -t 7 -p 4 -o trimmed-qualtrim --no-tail-n-trim --max-n=-1 --no-adapter-trim --no-length-filtration --quality-score $QSCORE --check-identifier
# rename --force "s/atria.fastq/qual$QSCORE.fastq/" trimmed-qualtrim/*fastq*

# trimmed-qualtrim/{atropos,fastp,ktrim,seqpurge}.R1*fastq
for i in trimmed-qualtrim/*.R1.qual$QSCORE.fastq
do
	echo $i
	mapping $i ${i/.R1./.R2.}
	samtools stats $i.bam > $i.bam.samtools-stats
done 2>&1 | tee mapping.log

cd trimmed-qualtrim
pasteSamtoolsStats *samtools-stats
cd ..




working_dir=~/analysis/atria-benchmark/SRR7243169


cd $working_dir

# Ensifer spp associated with Medicago whole genome sequencing
# Miseq 301bp paired reads
r1=SRR7243169_1.fastq
r2=SRR7243169_2.fastq
a1=CTGTCTCTTATACACATCT
a2=CTGTCTCTTATACACATCT
bwa_ref=genomes/Pseudomonas.sp.Z003-0.4C.fasta


. $atria/benchmark-aka3/trimming-functions.bash


atria-src

QSCORE=15
qualtrim Atria-src/${r1/.fastq}.atria.fastq Atria-src/${r2/.fastq}.atria.fastq

mapping_bowtie2 Atria-src/qualtrim/*.qual$QSCORE.fastq
bowtie2stat 15

mapping Atria-src/qualtrim/*.qual$QSCORE.fastq
for i in Atria-src/qualtrim/*.qual$QSCORE.fastq.*bam
do
    samtools stats $i > $i.samtools-stats
done

pasteSamtoolsStats */qualtrim/*{d,q}.bam.samtools-stats */qualtrim/*bowtie2.bam.samtools-stats

wps samtools-stats.collection

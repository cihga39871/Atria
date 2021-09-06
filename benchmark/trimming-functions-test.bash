
atria-trim-arg(){
	kmer_tolerance=$1
	diff1=$2  # pe-adapter-diff
	diff2=$3  # r1-r2-diff
	kmer_n_match=$4  # r1-r2-score-diff
	trim_score=$5
	tail_length=$6 # 8
	outdir="Atria-$(julia $atria/src/atria --version)-$1-$2-$3-$4-$5-$6"

	local r1=TAIR10.sim_1.fq
	local r2=TAIR10.sim_2.fq
	local a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
	local a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
	julia $atria/src/atria --no-consensus -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 -T $kmer_tolerance -d $diff1 -D $diff2 --kmer-n-match $kmer_n_match --trim-score-pe $trim_score --stats
	$atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab $outdir/TAIR10.sim_1.atria.fq $outdir/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.$outdir
	cat summary.$outdir
}
atria-trim-arg-330569(){
	kmer_tolerance=$1
	diff1=$2  # pe-adapter-diff
	diff2=$3  # r1-r2-diff
	kmer_n_match=$4  # r1-r2-score-diff
	trim_score=$5
	tail_length=$6 # 8
	outdir="Atria-$(julia $atria/src/atria --version)-$1-$2-$3-$4-$5-$6"

	local r1=SRR330569.3_1.fastq.gz
	local r2=SRR330569.3_2.fastq.gz
	local a1=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG
	local a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT
	local bwa_ref=`pwd`/genomes/dsim-all-chromosome-r2.02.fasta

	julia $atria/src/atria --no-consensus -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 -T $kmer_tolerance -d $diff1 -D $diff2 --kmer-n-match $kmer_n_match --trim-score-pe $trim_score --tail-length $tail_length
	i=$outdir/SRR330569.3_1.atria.fastq.gz
	mapping_hisat2 $i $outdir/SRR330569.3_2.atria.fastq.gz
	samtools stats $i.hisat2.sam > $i.hisat2.sam.samtools-stats
	rm $i.hisat2.sam
}
atria-trim-arg-7243169(){
	kmer_tolerance=$1
	diff1=$2  # pe-adapter-diff
	diff2=$3  # r1-r2-diff
	kmer_n_match=$4  # r1-r2-score-diff
	trim_score=$5
	tail_length=$6 # 8
	outdir="Atria-$(julia $atria/src/atria --version)-$1-$2-$3-$4-$5-$6"

	local r1=SRR7243169_1.fastq
	local r2=SRR7243169_2.fastq
	local a1=CTGTCTCTTATACACATCT
	local a2=CTGTCTCTTATACACATCT
	local bwa_ref=`pwd`/genomes/Pseudomonas.sp.Z003-0.4C.fasta

	julia $atria/src/atria --no-consensus -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 -T $kmer_tolerance -d $diff1 -D $diff2 --kmer-n-match $kmer_n_match --trim-score-pe $trim_score --tail-length $tail_length
	# $atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab $outdir/TAIR10.sim_1.atria.fq $outdir/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.$outdir
	# cat summary.$outdir
	local i=$outdir/SRR7243169_1.atria.fastq
	mapping $i $outdir/SRR7243169_2.atria.fastq
	samtools stats $i.bam > $i.bam.samtools-stats
	rm $i.bam
}
atria-trim-arg-human(){
	kmer_tolerance=$1
	diff1=$2  # pe-adapter-diff
	diff2=$3  # r1-r2-diff
	kmer_n_match=$4  # r1-r2-score-diff
	trim_score=$5
	tail_length=$6 # 8
	outdir="Atria-$(julia $atria/src/atria --version)-$1-$2-$3-$4-$5-$6"

	r1=ERR4695159_1.fastq.gz
	r2=ERR4695159_2.fastq.gz
	a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
	a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
	bwa_ref=`pwd`/genomes/hg38.fasta.gz

	julia $atria/src/atria --no-consensus -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 -T $kmer_tolerance -d $diff1 -D $diff2 --kmer-n-match $kmer_n_match --trim-score-pe $trim_score --stats --compress no --tail-length $tail_length
	# $atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab $outdir/TAIR10.sim_1.atria.fq $outdir/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.$outdir
	# cat summary.$outdir
	local i=$outdir/ERR4695159_1.atria.fastq
	mapping_bowtie2 $i $outdir/ERR4695159_2.atria.fastq
	samtools stats $i.bowtie2.bam > $i.bowtie2.bam.samtools-stats
	rm $i.bowtie2.bam
}

atria-trim-arg-simulate(){
	kmer_tolerance=$1 # 2
	diff1=$2  # pe-adapter-diff 0
	diff2=$3  # r1-r2-diff 0
	kmer_n_match=$4  # 9
	trim_score=$5 # 10
	tail_length=$6 # 12
	outdir="Atria-$(julia $atria/src/atria --version)-$1-$2-$3-$4-$5-$6"

	julia $atria/src/atria --no-consensus -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 -T $kmer_tolerance -d $diff1 -D $diff2 --kmer-n-match $kmer_n_match --trim-score-pe $trim_score --tail-length $tail_length
	# $atria/benchmark/evalTrimming.pl 200 TAIR10.sim_1.length.tab $outdir/TAIR10.sim_1.atria.fq $outdir/TAIR10.sim_2.atria.fq TAIR10.sim_2.length.tab > summary.$outdir
	# cat summary.$outdir
	parallel -j 2 atria readstat ::: $outdir/*fastq
}

atria-consensus-arg(){
	outdir="Atria-consensus-$1-$2-$3"
	$atria/src/atria -r $r1 -R $r2 -o $outdir --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 --overlap-score $1 --min-ratio-mismatch $2 --prob-diff $3
}
atria-src(){
	export JULIA_NUM_THREADS=15
	$atria/src/atria --no-consensus -r $r1 -R $r2 -o Atria-src --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 --stats
}
atria-consensus-src(){
		export JULIA_NUM_THREADS=15
    $atria/src/atria \
        -r $r1 -R $r2 \
        -o Atria-consensus-src \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration --adapter1 $a1 --adapter2 $a2 --stats
}

mapping_noendclip() {
    bwa mem -v 1 -t 25 -L 50,50 $bwa_ref $1 $2 |\
        samtools view -@ 20 -b -o $1.noendclip.bam
}


qualtrim_trimmomatic(){
	if [[ $1 = *gz ]]
	then
		local num_threads=20
		local gzext=.gz
	else
		local num_threads=8
		local gzext=
	fi
	if [[ $3 ]]
	then
		local num_threads=$3
	fi

	local OUTDIR=`dirname "$1"`/../trimmed-qualtrim_3rd_party
	mkdir -p $OUTDIR

	if [[ $r1 = *gz ]]; then
		local isgz=.gz
	else
		local isgz=
	fi

	local output_r1=$OUTDIR/`basename "${1/.fastq*/}"`.qual$QSCORE.fq$gzext
	local output_r1_up=$OUTDIR/`basename "${1/.fastq*/}"`.qual$QSCORE.unpair.fq$gzext
	local output_r2=$OUTDIR/`basename "${2/.fastq*/}"`.qual$QSCORE.fq$gzext
	local output_r2_up=$OUTDIR/`basename "${2/.fastq*/}"`.qual$QSCORE.unpair.fq$gzext

	/usr/bin/time -v java -jar /usr/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $num_threads -phred33 $1 $2 $output_r1 $output_r1_up $output_r2 $output_r2_up SLIDINGWINDOW:5:$QSCORE MINLEN:0
}


samalign(){
	samtools view $1 | awk '{print $1"\t"$3"\t"$4"\t"$6}'
}

sampaste(){
	paste <(samtools view $1 | awk '{print $1"\t"$3"\t"$4"\t"$6}') <(samtools view $2 | awk '{print $1"\t"$3"\t"$4"\t"$6}') | less -SN
}

akadiff(){
	diff <(samalign $1) <(samalign $2) -y | less
}

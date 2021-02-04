#!bash
export JULIA_NUM_THREADS=16

run_atria(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v atria --no-consensus \
        -r $r1 -R $r2 \
        -o Atria \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads
}

run_atria_consensus(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    /usr/bin/time -v atria \
        -r $r1 -R $r2 \
        -o Atria-consensus \
        --no-tail-n-trim --max-n=-1 --no-quality-trim --no-length-filtration \
        --adapter1 $a1 --adapter2 $a2 --threads $num_threads
}

run_adapterremoval() {
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    local err=3
    local folder="AdapterRemoval-$err"
	mkdir -p "$folder"
	if [[ $r1 = *gz ]]; then
		/usr/bin/time -v AdapterRemoval --file1 $r1 --file2 $r2 \
	        --basename "$folder"/adapterremoval \
	        --adapter1 $a1 --adapter2 $a2 \
	        --mm $err --minlength 0 --threads $num_threads --gzip
	else
		/usr/bin/time -v AdapterRemoval --file1 $r1 --file2 $r2 \
			--basename "$folder"/adapterremoval \
			--adapter1 $a1 --adapter2 $a2 \
			--mm $err --minlength 0 --threads $num_threads
	fi
}

run_skewer(){
    local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    local OUTDIR="Skewer"
    mkdir -p $OUTDIR
	if [[ $r1 = *gz ]]; then
	    /usr/bin/time -v skewer --quiet \
	        -x $a1 -y $a2 -m pe  \
	        -l 0 -o $OUTDIR/$OUTDIR $r1 $r2 --threads $num_threads --compress
	else
		/usr/bin/time -v skewer --quiet \
	        -x $a1 -y $a2 -m pe  \
	        -l 0 -o $OUTDIR/$OUTDIR $r1 $r2 --threads $num_threads
	fi
}

run_trim_galore(){
	local num_threads=1
    if [[ $1 ]]; then
        num_threads=$1
    fi
    local OUTDIR="TrimGalore"
	mkdir -p $OUTDIR
	/usr/bin/time -v trim_galore --cores $num_threads \
	    --quality 0 \
	    -o $OUTDIR \
	    --adapter $a1 \
	    --adapter2 $a2 \
	    -e 0.1 --stringency 1 \
	    --max_n 100 --length 0 \
	    --paired $r1 $r2
}

run_trimmomatic(){
	local num_threads=1
	if [[ $1 ]]; then
		num_threads=$1
	fi
	local OUTDIR="Trimmomatic"
	mkdir -p $OUTDIR
	rm -f adapters.fa
	echo '>TruSeq3/1' >> adapters.fa
	echo $a1 >> adapters.fa
	echo '>TruSeq3/2' >> adapters.fa
	echo $a2 >> adapters.fa
	output=$OUTDIR/out
	if [[ $r1 = *gz ]]; then
		local isgz=.gz
	else
		local isgz=
	fi

	/usr/bin/time -v java -jar /usr/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $num_threads -phred33 $r1 $r2 $output-pair1.paired.fq$isgz $output-pair1.unpaired.fq$isgz $output-pair2.paired.fq$isgz $output-pair2.unpaired.fq$isgz ILLUMINACLIP:adapters.fa:2:30:10:1:TRUE:keepBothReads MINLEN:1
}

run_ktrim(){
	local num_threads=1
	if [[ $1 ]]; then
		num_threads=$1
	fi
	local OUTDIR="Ktrim"
	mkdir -p $OUTDIR
	/usr/bin/time -v ktrim -1 $r1 -2 $r2 -t $num_threads -p 33 -q 1 -s 10 -a $a1 -b $a2 -o Ktrim/ktrim
}

run_fastp(){
	local num_threads=1
	if [[ $1 ]]; then
		num_threads=$1
	fi
	local OUTDIR="fastp"
	output=$OUTDIR/out.fastp
	if [[ $r1 = *gz ]]; then
		local isgz=.gz
	else
		local isgz=
	fi
	mkdir -p $OUTDIR
	/usr/bin/time -v fastp --in1 $r1 --in2 $r2 --out1 $output.r1.fq$isgz --out2 $output.r2.fq$isgz \
		-z 6 --adapter_sequence $a1 --adapter_sequence_r2 $a2 --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --thread $num_threads
}

run_seqpurge() {
	local num_threads=1
	if [[ $1 ]]; then
		num_threads=$1
	fi
    local folder=SeqPurge
    mkdir -p "$folder"
	# output always gziped
    /usr/bin/time -v SeqPurge -in1 $r1 -in2 $r2 -out1 "$folder"/$r1.seqpurge.fq.gz -out2 "$folder"/$r2.seqpurge.fq.gz \
        -a1 $a1 -a2 $a2 -mep 0.1 \
        -qcut 0 -min_len 0 -summary "$folder"/seqpurge.summary -threads $num_threads
}

run_atropos() {
	# Atropos 1.1.29 with Python 3.8.5
	local num_threads=1
	if [[ $1 ]]; then
		num_threads=$1
	fi
	if [[ $r1 = *gz ]]; then
		local isgz=.gz
	else
		local isgz=
	fi
    local folder="Atropos"
    mkdir -p "$folder"
	if [[ $num_threads == 1 ]]
	then
	    /usr/bin/time -v atropos trim -a $a1 -A $a2 \
	        -o "$folder"/$r1.atropos.fq$isgz -p "$folder"/$r2.atropos.fq$isgz -pe1 $r1 -pe2 $r2 \
	        --aligner insert -e 0.1
	else
		/usr/bin/time -v atropos trim -a $a1 -A $a2 \
	        -o "$folder"/$r1.atropos.fq$isgz -p "$folder"/$r2.atropos.fq$isgz -pe1 $r1 -pe2 $r2 \
	        --aligner insert -e 0.1 --threads $num_threads --preserve-order
	fi
}


mapping() {
    bwa mem -v 1 -t 25 $bwa_ref $1 $2 |\
	samtools view -@ 10 -b -o $1.bam
}

mapping_bowtie2(){
	bowtie2 --maxins 800 --threads 25 -x $bwa_ref-bowtie2 -1 $1 -2 $2 2> $1.bowtie2.stat |\
	samtools view -@ 10 -b -o $1.bowtie2.bam
}
mapping_hisat2(){
	hisat2 --threads 25 -x $bwa_ref-hisat2 -1 $1 -2 $2 -S $1.hisat2.sam 2> $1.hisat2.stat
}
qualtrim(){
		local DIR=`dirname "$1"`/../trimmed-qualtrim
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
		time atria -r "$1" -R "$2" -t $num_threads --check-identifier \
		-o "$DIR" \
		--no-tail-n-trim --max-n=-1 --no-adapter-trim --no-length-filtration \
		--quality-score $QSCORE
		rename --force "s/atria.fastq/qual$QSCORE.fastq/" "$DIR"/*fastq$gzext
		rename --force "s/atria.fq/qual$QSCORE.fq/" "$DIR"/*fq$gzext
		rename --force "s/atria.truncated/qual$QSCORE.truncated/" "$DIR"/*truncated$gzext
		rename --force "s/atria.log/qual$QSCORE.log/" "$DIR"/*log*
}

bowtie2stat(){
	if [[ $1 ]]
	then
		local QSCORE=$1
	else
		local QSCORE=
	fi
	grep -v Warning */*qualtrim/*qual$QSCORE*bowtie2.stat | sed 's#/[^:]*#\t#' | grep "aligned concordantly exactly 1 time" | column -ts$'\t'
	echo
	grep -v Warning */*qualtrim/*qual$QSCORE*bowtie2.stat | sed 's#/[^:]*#\t#' | grep "aligned 0 times concordantly or discordantly" | column -ts$'\t'
}

pasteSamtoolsStats(){
    grep ^SN $1| cut -f 2,4 | sed 's/\t# \(.*\)/ [\1]/' | sed 's/://' | awk 'BEGIN{print "sample"};{print}' > samtools-stats.collection.txt
    for i in "$@"
    do
        paste samtools-stats.collection.txt <(grep ^SN $i| cut -f 3 | awk -v var=${i/.fastq*/} 'BEGIN{print var};{print}') > samtools-stats.collection.tmp
        mv samtools-stats.collection.tmp samtools-stats.collection.txt
    done
    echo Output: samtools-stats.collection.txt
}

pasteTimeOutput(){
    paste \
        <(grep -E "Command being timed" $1 | sed 's/.*Command being timed://') \
        <(grep -E "^\sUser time" $1) \
        <(grep -E "^\sSystem time" $1) \
        <(grep -E "^\sPercent of CPU this job got" $1) \
        <(grep -E "^\sElapsed" $1) \
        <(grep -E "^\sMaximum resident set size" $1)
}


sam2bam(){
	for i in "$@"
	do
		echo `date` - $i
		samtools view -b $i > ${i:0:-3}bam
		if [[ $? == 0 ]]
		then
			rm $i
		else
			rm ${i:0:-3}bam
			echo SamToBam failed: $i
		fi
	done
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

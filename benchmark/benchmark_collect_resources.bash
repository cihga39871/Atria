#!bash

pasteTimeOutput(){
    paste \
        <(grep -E "Command being timed" $1 | sed 's/.*Command being timed://') \
        <(grep -E "^\sUser time" $1) \
        <(grep -E "^\sSystem time" $1) \
        <(grep -E "^\sPercent of CPU this job got" $1) \
        <(grep -E "^\sElapsed" $1) \
        <(grep -E "^\sMaximum resident set size" $1)
}

pasteSamtoolsStats(){
    grep ^SN $1| cut -f 2,4 | sed 's/\t# \(.*\)/ [\1]/' | sed 's/://' | awk 'BEGIN{print "sample"};{print}' > samtools-stats.collection
    for i in "$@"
    do
        paste samtools-stats.collection <(grep ^SN $i| cut -f 3 | awk -v var=${i/.fastq*/} 'BEGIN{print var};{print}') > samtools-stats.collection.tmp
        mv samtools-stats.collection.tmp samtools-stats.collection
    done
    echo Output: samtools-stats.collection
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

samsort(){
	samtools sort -@ 12 $1 -Obam -n -o ${1/bam/sort.bam}
}

#!/bin/bash

## Script that takes a representative read, extracts the inserted sequence with 2000 bp flanking it both sides, maps supporting reads to it
## and utilizes racon to polish the sequence.


line=$1		#line of vcf containing the insertion
outputdir=$2	#directory where directory redone is
reference=$3	#reference genome file
index=$4	#name for directory
racon=$5 #path for racon build

if [[ $line == \#* ]]; then
	continue
fi
id=$(echo $line | cut -f3 -d' ')
echo $id
echo $index
mkdir -p $outputdir/redone
mkdir -p $outputdir/redone/$index
rm $outputdir/redone/$index/$id.fq
for f in $outputdir/reads_in_insertions/*/$id.txt;	#extract supporting reads and makes a fq file
	do
	sample=$(echo ${f} | rev | cut -f1,2 -d '/' | rev | cut -f1 -d'/')
	echo $f
	echo  $outputdir/$sample".fastq"
	grep -A3 -f $f <(samtools fastq $outputdir/$sample*"_rm.bam") | sed '/^--$/d' >> $outputdir/redone/$index/$id.fq
	done;
echo $outputdir/reads_in_insertions/*/$id".txt_largest.txt"
for best in $outputdir/reads_in_insertions/*/$id".txt_largest.txt";	#get the median length read as representative read
	do
	grep median $best | grep -w -v "." >> $outputdir/redone/$index/$id"_best_reads.txt"
	done;
filesize=$(wc -l $outputdir/redone/$index/$id"_best_reads.txt" | cut -f1 -d' ')
if [[ $filesize == 0 ]]
	then
        continue
	fi
half=$((filesize / 2))
if [[ $half == 0 ]]
	then
        half=1
        fi
sequence_name=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 2)
echo $outputdir"/redone/fastq/"$id"_best_reads.txt"
echo $sequence_name
sequence=$(echo $line | cut -f5 -d' ')
chr=$(echo $line | cut -f1 -d' ')
start=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 4)
insertion=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 3)
printf ">" > $outputdir/redone/$index/$id"_1.fasta"
echo $id >> $outputdir/redone/$index/$id"_1.fasta"
printf ">" > $outputdir/redone/$index/$id"_read.fasta"
echo $sequence_name >> $outputdir/redone/$index/$id"_read.fasta"
echo  $sequence_name $outputdir"/redone/"$index/$id".fq"
grep -A1 $sequence_name $outputdir/redone/$index/$id.fq | tail -n1 >> $outputdir/redone/$index/$id"_read.fasta"


#extract sequence form representative read

start_before="$(($start-2000))"
start_after="$(($start+$insertion+2000))"
if [[ $start_before -lt 1 ]]
	then
	start_before=1
	fi
samtools faidx $outputdir/redone/$index/$id"_read.fasta"
echo $sequence_name":"$start_before"-"$start_after
samtools faidx $outputdir/redone/$index/$id"_read.fasta" $sequence_name":"$start_before"-"$start_after | grep -v ">" >> $outputdir/redone/$index/$id"_1.fasta"

#map reads to template and run racon, repeat in total 3 times

COUNTER=1
while [  $COUNTER -lt 4 ]; do
	bwa index $outputdir/redone/$index/$id"_"$COUNTER.fasta
	bwa mem $outputdir/redone/$index/$id"_"$COUNTER.fasta $outputdir/redone/$index/$id.fq > $outputdir/redone/$index/$id"_"$COUNTER".sam"
	$racon $outputdir/redone/$index/$id.fq $outputdir/redone/$index/$id"_"$COUNTER".sam" $outputdir/redone/$index/$id"_"$COUNTER.fasta | cut -f1 -d' '> $outputdir/redone/$index/$id"_"$COUNTER"_new".fasta

	mv $outputdir/redone/$index/$id"_"$COUNTER"_new".fasta $outputdir/redone/$index/$id"_""$(($COUNTER + 1))"".fasta"
	let COUNTER=COUNTER+1
	done

echo ">"$id > $outputdir/redone/$index"/"$id".fasta"
grep -v ">" $outputdir/redone/$index/$id"_""$(($COUNTER))"".fasta" >> $outputdir/redone/$index"/"$id".fasta"


#map polished read to reference genome

bwa mem $reference $outputdir/redone/$index/$id".fasta" >  $outputdir"/redone/"$index/$id"_aligned.sam"

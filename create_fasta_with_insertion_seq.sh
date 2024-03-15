#!/bin/bash

## Script that takes insertion read, extracts the flanking sequence from reference genome, maps supporting reads to it
## and utilizes racon to polish the sequence.


line=$1		#line of vcf containing the insertion
outputdir=$2	#directory where reads_in_insertions
reference=$3	#reference genome file
racon=$4	#racon build

mkdir $outputdir"/fastq"

if [[ $line == \#* ]]; then
	continue
fi
id=$(echo $line | cut -f3 -d' ')
echo $id
rm $outputdir/fastq/$id.fq

for f in $outputdir/reads_in_insertions/*/$id.txt;	#extract supporting reads and makes a fq file
	do
	sample=$(echo ${f} | rev | cut -f1,2 -d '/' | rev | cut -f1 -d'/')
	grep -A3 -f $f <(samtools fastq $outputdir/$sample*"_rm.bam") | sed '/^--$/d' >> $outputdir/fastq/$id.fq
	done
sequence=$(echo $line | cut -f5 -d' ')
chr=$(echo $line | cut -f1 -d' ')
start=$(echo $line | cut -f2 -d' ')
printf ">" > $outputdir/fastq/$id"_1.fasta"
echo $id >> $outputdir/fastq/$id"_1.fasta"

#extract flanking sequence from refence genome

start_before="$(($start-2000))"
start_after="$(($start+2000))"


samtools faidx $reference $chr":"$start_before"-"$start | grep -v ">" >> $outputdir/fastq/$id"_1.fasta"
printf $sequence >> $outputdir/fastq/$id"_1.fasta"
samtools faidx $reference $chr":"$start"-"$start_after | grep -v ">" >> $outputdir/fastq/$id"_1.fasta"

#map reads to template and run racon, repeat in total 3 times


COUNTER=1
while [  $COUNTER -lt 4 ]; do
	bwa index $outputdir/fastq/$id"_"$COUNTER.fasta
	bwa mem $outputdir/fastq/$id"_"$COUNTER.fasta $outputdir/fastq/$id.fq > $outputdir/fastq/$id"_"$COUNTER".sam"
	$racon $outputdir/fastq/$id.fq $outputdir/fastq/$id"_"$COUNTER".sam" $outputdir/fastq/$id"_"$COUNTER.fasta | cut -f1 -d' '> $outputdir/fastq/$id"_"$COUNTER"_new".fasta

	mv $outputdir/fastq/$id"_"$COUNTER"_new".fasta $outputdir/fastq/$id"_""$(($COUNTER + 1))"".fasta"
	let COUNTER=COUNTER+1
	done

echo ">"$id > $outputdir/fastq/$id".fasta"
grep -v ">" $outputdir/fastq/$id"_""$(($COUNTER))"".fasta" >> $outputdir/fastq/$id".fasta"

#map polished read to reference genome

bwa mem $reference $outputdir/fastq/$id".fasta" > $outputdir"/fastq/"$id"_aligned.sam"

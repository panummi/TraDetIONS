#!/bin/bash

path=$1
output=$2
reference=$3
racon=$4
repeats=$5
first=$6
vcf=$7

while IFS= read -r line
	do
	if [[ $line == \#* ]]; then
		continue
	fi
	linepart=$(echo $line | cut -f1-7 -d' ')
	index=$(echo $1 | cut -f6 -d'/')
	        id=$(echo $line | cut -f3 -d' ')
        chr=$(echo $line | cut -f1 -d' ')
        start=$(echo $line | cut -f2 -d' ')
        [[ $chr =~ ^#.* ]] && continue

	if [[ $first == "T" ]]; then
		echo "H"
		mkdir -p $output"/annotation"
		bash $path/create_fasta_with_insertion_seq.sh "$linepart" $output/reads/ $reference $racon
		python $path/insertion_annotation.py -f $output"/reads/fastq/"$id".fasta" -s $output"/reads/fastq/"$id"_aligned.sam" -p $chr":"$start -o $output"/annotation/"$id.bed -r $repeats -g $reference
	else
		mkdir -p $output"/reads/redone/annotation"
		bash $path/create_fasta_with_read.sh "$linepart" $output/reads/ $reference fastq/ $racon
		python $path/insertion_annotation.py -f $output"/reads/redone/fastq/"$id".fasta" -s $output"/reads/redone/fastq/"$id"_aligned.sam" -p $chr":"$start -o $output"/reads/redone/annotation/"$id.bed -r $repeats -g $reference
	fi
done < $vcf

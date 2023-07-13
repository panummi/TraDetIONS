#!/bin/bash

#script that removes duplicates of reads
$1=dir

for f in $dir/*sorted.map-ont.bam; do
	sample=$(echo $f | cut -f1,2,3 -d'.')
	samtools rmdup -s $f $sample"_rm.bam"
done

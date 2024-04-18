#!/bin/bash

#script that removes duplicates of reads

dir=$1

echo $dir
for f in $dir/*bam; do
	sample=$(echo $f | cut -f1,2,3 -d'.')
	samtools rmdup -s $f $sample"_rm.bam"
	rm $f
done

#!/bin/bash

#samples=$("ls /mnt/cg38/Paivi/sv_analysis/reads/*bam")

for f in /mnt/cg38/Paivi/new_runs/pseudogenes/reads/*sorted.map-ont.bam; do
	sample=$(echo $f | cut -f1,2,3 -d'.')
	echo $sample

samtools rmdup -s $f $sample"_rm.bam"

done

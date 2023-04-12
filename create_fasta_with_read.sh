#!/bin/bash
#input="/mnt/cg38/Paivi/sniffles_data/somatic_filtering/germline_insertions_new_ids.vcf"
#input=$1
line=$1
outputdir=$2
reference=$3
index=$4
echo $input


#mkdir $outputdir"/redone/"$index


#while IFS= read -r line
#	do
	if [[ $line == \#* ]]; then
		continue
	fi
	id=$(echo $line | cut -f3 -d' ')
	echo $id
#	rm $outputdir/redone/$index/$id.fq
	rm $outputdir/redone/$index/$id.fq
	for f in $outputdir/reads_in_insertions/*/$id.txt;
		do
		sample=$(echo ${f} | rev | cut -f1,2 -d '/' | rev | cut -f1 -d'/')
		echo $f
#		samtools rmdup -s $outputdir/$sample".sorted.map-ont.bam" $outputdir/$sample".sorted.map-ont_rm.bam"
		echo  $outputdir/$sample".fastq"
		grep -A3 -f $f <(samtools fastq $outputdir/$sample".sorted.map-ont_rm.bam") | sed '/^--$/d' >> $outputdir/redone/$index/$id.fq
#		seqtk subseq $outputdir/$sample".fastq" $f >> $outputdir/redone/$index/$id.fq
#		grep -A3 -f $f <(samtools fastq $outputdir/$sample".sorted.map-ont_rm.bam")
		done;
	echo $outputdir/reads_in_insertions/*/$id".txt_largest.txt"
        for best in $outputdir/reads_in_insertions/*/$id".txt_largest.txt";
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
	echo "THIS IN HERE"
	sequence_name=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 2)
	echo $outputdir"/redone/fastq/"$id"_best_reads.txt"
	echo $sequence_name
	sequence=$(echo $line | cut -f5 -d' ')
	chr=$(echo $line | cut -f1 -d' ')
#        start=$(echo $line | cut -f2 -d' ')
	start=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 4)
	insertion=$(sort -k3,3n $outputdir/redone/$index/$id"_best_reads.txt" | head -n $half | tail -n 1 | cut -f 3)
#	echo $start
    printf ">" > $outputdir/redone/$index/$id"_1.fasta"
	echo $id >> $outputdir/redone/$index/$id"_1.fasta"
#	echo 
echo "HERE"        
#       start_before="$(($start-2000))"
#       start_after="$(($start+2000))"


#	samtools faidx $reference $chr":"$start_before"-"$start | grep -v ">" >> $outputdir/redone/fastq/$id"_1.fasta"
#	printf $sequence >> $outputdir/redone/fastq/$id"_1.fasta"
#	samtools faidx $reference $chr":"$start"-"$start_after | grep -v ">" >> $outputdir/redone/fastq/$id"_1.fasta"
#    grep -A1 $sequence_name $outputdir/fastq/$id.fq | tail -n1 >> $outputdir/redone/fastq/$id"_1.fasta"
    printf ">" > $outputdir/redone/$index/$id"_read.fasta"
    echo $sequence_name >> $outputdir/redone/$index/$id"_read.fasta"
#	grep -A1 $sequence_name $outputdir/fastq/$id.fq | tail -n1
    echo  $sequence_name $outputdir"/redone/"$index/$id".fq"
    grep -A1 $sequence_name $outputdir/redone/$index/$id.fq | tail -n1 >> $outputdir/redone/$index/$id"_read.fasta"
	
	start_before="$(($start-2000))"
    start_after="$(($start+$insertion+2000))"
	if [[ $start_before -lt 1 ]]
	then
	    start_before=1
	fi
#ECHO "här"
	samtools faidx $outputdir/redone/$index/$id"_read.fasta"
	echo $sequence_name":"$start_before"-"$start_after
	samtools faidx $outputdir/redone/$index/$id"_read.fasta" $sequence_name":"$start_before"-"$start_after | grep -v ">" >> $outputdir/redone/$index/$id"_1.fasta"

	COUNTER=1
        while [  $COUNTER -lt 4 ]; do
		bwa index $outputdir/redone/$index/$id"_"$COUNTER.fasta
		bwa mem $outputdir/redone/$index/$id"_"$COUNTER.fasta $outputdir/redone/$index/$id.fq > $outputdir/redone/$index/$id"_"$COUNTER".sam"
		/mnt/cg8/Päivi/software/racon/build/bin/racon $outputdir/redone/$index/$id.fq $outputdir/redone/$index/$id"_"$COUNTER".sam" $outputdir/redone/$index/$id"_"$COUNTER.fasta | cut -f1 -d' '> $outputdir/redone/$index/$id"_"$COUNTER"_new".fasta

		mv $outputdir/redone/$index/$id"_"$COUNTER"_new".fasta $outputdir/redone/$index/$id"_""$(($COUNTER + 1))"".fasta"
		let COUNTER=COUNTER+1
		done

	echo ">"$id > $outputdir/redone/$index/$id".fasta"
	grep -v ">" $outputdir/redone/$index/$id"_""$(($COUNTER))"".fasta" >> $outputdir/redone/$index/$id".fasta"

#	bwa mem /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_reference.fasta" /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_"$COUNTER".fasta" > "/mnt/cg38/Paivi/sniffles_data/reads/fastq/"$id"_aligned.sam"
	bwa mem /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta $outputdir/redone/$index/$id".fasta" >  $outputdir"/redone/"$index/$id"_aligned.sam"

#	coordinates=$(python get_insertion_sequence.py -i "/mnt/cg38/Paivi/sniffles_data/reads/fastq/"$id"_aligned.sam")
#	start_new = $($coordinates | cut -f 1 -d'-' | cut -f2 -d':')
#	echo $start_new
#	echo $coordinates
#	printf ">" >> "/mnt/cg38/Paivi/sniffles_data/reads/fastq/"$id"_aligned.sam"
#	echo $(($start_new+$start_before)) >> "/mnt/cg38/Paivi/sniffles_data/reads/fastq/"$id"_aligned.sam"
#	samtools faidx /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_"$COUNTER".fasta" $coordinates > /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_"$COUNTER"_insertion.fa"
	#bwa mem /mnt/cg8/Päivi/databases/hg38reps.html /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_"$COUNTER"_insertion.fa" > /mnt/cg38/Paivi/sniffles_data/reads/fastq/$id"_"$COUNTER"_insertion.sam"

#done < "$input"

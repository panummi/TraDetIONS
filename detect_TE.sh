#!/bin/bash

tradetionspath=$1
vcfinput=$2
repeatfasta=$3
samplefile=$4
outputdir=$5
referencegenome=$6
racondir=$7
TE=$8
origvcf=$9

mkdir -p $outputdir
mkdir -p $outputdir/reads

python $tradetionspath/survivor_retro_annotate.py -i $vcfinput -a $repeatfasta -o $outputdir/TE_output
python $tradetionspath/mark_variants.py -s $samplefile -v $outputdir/TE_output.vcf -o $outputdir/marked_variants.vcf
cat <(grep "#" $outputdir/marked_variants.vcf) <(awk '{print $1 "\t" $2 "\t" "TE_"NR "\t" $0}' <(grep -v "#" $outputdir/marked_variants.vcf) | cut -f1,2,3,7-) > $outputdir/preliminary_TE.vcf
python $tradetionspath/get_insertion_reads.py -i $outputdir/preliminary_TE.vcf -o $outputdir/reads/ -s $samplefile
bash $tradetionspath/rmdups.sh $outputdir/reads
bash $tradetionspath/polish_and_annotate.sh $tradetionspath $outputdir $referencegenome $racondir $repeatfasta "T" $outputdir/preliminary_TE.vcf
cat $outputdir/annotation/*.bed > $outputdir/annotation/annotation_first.bed
grep -v "#" $outputdir/preliminary_TE.vcf | cut -f3 > $outputdir/TE_ids.txt
bash $tradetionspath/annotation_check.sh $outputdir/annotation/annotation_first.bed $outputdir/preliminary_TE.vcf $outputdir/TE_ids.txt $TE > $outputdir/first_annotation.tsv
awk '{if ($6 != "Match") print $1}' $outputdir/first_annotation.tsv > $outputdir/redo_ids.txt
grep -w -f $outputdir/redo_ids.txt $outputdir/preliminary_TE.vcf > $outputdir/redo.vcf
bash $tradetionspath/polish_and_annotate.sh $tradetionspath $outputdir $referencegenome $racondir $repeatfasta "F" $outputdir/redo.vcf
cat $outputdir/reads/redone/annotation/*.bed > $outputdir/reads/redone/annotation/annotation_second.bed
bash $tradetionspath/annotation_check.sh $outputdir/reads/redone/annotation/annotation_second.bed $outputdir/preliminary_TE.vcf $outputdir/redo_ids.txt $TE > $outputdir/second_annotation.tsv

awk '{if ($6 == "Match" && $7 == "hallmarks") print $1}' $outputdir/first_annotation.tsv > $outputdir/good_ids_1.txt
awk '{if ($6 == "Match" && $7 == "hallmarks") print $1}' $outputdir/second_annotation.tsv > $outputdir/good_ids_2.txt

mkdir -p $outputdir"/TE"
cat <(grep "#" $outputdir/preliminary_TE.vcf ) <(grep -w -f <(cat $outputdir/good_ids_1.txt $outputdir/good_ids_2.txt) $outputdir/preliminary_TE.vcf) > $outputdir"/TE/TE.vcf"
while IFS= read -r line
	do
		cat $outputdir/reads/fastq/$line.fasta >> $outputdir"/TE/TE_insertions.fa"
		cat $outputdir/annotation/$line.bed >> $outputdir"/TE/TE_annotation.bed"
	done<$outputdir/good_ids_1.txt
while IFS= read -r line
        do
                cat $outputdir/reads/redone/fastq/$line.fasta >> $outputdir"/TE/TE_insertions.fa"
		cat $outputdir/reads/redone/annotation/$line.bed >> $outputdir"/TE/TE_annotation.bed"
        done<$outputdir/good_ids_2.txt

######somatic filtering

sed 's/;S=/\t/g' $outputdir"/TE/TE.vcf" | cut -f3,9 | cut -f1 -d';' | awk '{if ($2 == "1") print $1}' > $outputdir"/candidate_somatic.txt"

cat <(grep "#" $outputdir/preliminary_TE.vcf ) <(grep -w -f $outputdir"/candidate_somatic.txt" $outputdir/preliminary_TE.vcf) > $outputdir"/candidate_somatic.vcf"
bedtools window -w 10 -v -a $outputdir"/candidate_somatic.vcf" -b $origvcf > $outputdir"/somatic.vcf"

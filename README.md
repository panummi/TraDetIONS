# TraDetIONS

Transposon detection in Oxford Nanopore Sequencing data

Created by Päivi Nummi

Example data HG002 from https://human-pangenomics.s3.amazonaws.com/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/

Example run with this data:
        Download:
        aws s3 sync   --exclude='*.bam*' --exclude='*.fast5'    --no-sign-request s3://ont-open-data/gm24385_2020.11/flowcells/20201026_1644_2-E5-H5_PAG07162_d7f262d5/ $PWD/20201026_1644_2-E5-H5_PAG07162_d7f262d5/
        
        Map:
        minimap2 -x map-ont -a -2 --MD -R '@RG\tID:HG002\tSM:HG002\tPL:ONT' -t 27 reference_genome.fasta  $PWD//20201026_1644_2-E5-H5_PAG07162_d7f262d5/fastq_pass/*.fastq.gz|samtools sort -m 4G -@ 3 --reference reference_genome.fasta -O cram -o align/HG002.sorted.map-ont.cram
        
        Index:
        samtools index align/HG002.sorted.map-ont.cram
        
        SV calling:
        sniffles --vcf {output} --mapped_reads {input} --threads {threads} " --min_support 2 " --min_length  10 " --num_reads_report 3 " --min_seq_size 1000 --genotype  --cluster --report_seq
        
        SV merging:
        echo "HG002.sorted.map-ont.sniffles.vcf" > input.lst 
        SURVIVOR merge input.lst 50 -1 1 1 1 -1 survivor.vcf
        
        TraDetIONS:
        echo "HG002	N	-	1	HG002	HG002	align/HG002.sorted.map-ont.cram" > sample.txt
        bash detect_TE.sh tradetions_repo_path survivor.vcf data/hg38reps.fa sample.txt outputHG002 reference_genome.fasta racon_path data/TE_merge_subclass.tsv
        

TraDetIONS uses SVs called by Sniffles and merged with Survivor. Other dependencies are Racon and conda environment environment.yml.

INSTALL:

git clone https://github.com/panummi/TraDetIONS.git


EXAMPLE RUN:

prepairing svs:

sniffles --vcf {output} --mapped_reads {input}
        --threads {threads} "
        --min_support 2 "
        --min_length  10 "
        --num_reads_report 3 "
        --min_seq_size 1000
        --genotype  --cluster --report_seq
        
SURVIVOR merge input.lst 50 -1 1 1 1 -1 output.vcf


Get TEs:

bash detect_TE.sh tradetions_repo_path survivor.vcf data/hg38reps.fa data/sample.txt outputdir reference_genome.fasta racon_path data/TE_merge_subclass.tsv

Output will be in output/TE/*
and contain fasta file for insertions sequences, bed file for annotation of the sequences and vcf file for the variants.

TraDetIONS might overnight on HPC cluster

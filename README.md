# TraDetIONS

Transposon detection in Oxford Nanopore Sequencing data

Created by Päivi Nummi

Example data HG002 from https://human-pangenomics.s3.amazonaws.com/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/

TraDetIONS uses SVs called by Sniffles and merged with Survivor. Other dependencies are Racon and conda environment environment.yml.

The tools consist of following scripts, that perform actions:

EXAMPLE RUN:

prepairing svs:

sniffles --vcf {output} --mapped_reads {input}
        --threads {threads} "
        --min_support 2 "
        --min_length  10 "
        --num_reads_report 3 "
        --min_seq_size 1000
        --genotype  --cluster --report_seq



Get TEs:

bash detect_TE.sh tradetions_repo_path survivor.vcf /repepat/masker/data/hg38reps.fa data/sample.txt outputdir reference_genome.fasta racon_path TE_merge.tsv


Get orphan transductions and Pseudogenes:

bash detect_orphan_and_pseudogene.sh tradetions_repo_path survivor.vcf /repepat/masker/data/hg38reps.fa data/sample.txt polyA_outputdir reference_genome.fasta racon_path TE_merge.tsv A output/TE/TE.vcf data/Homo_sapiens.GRCh38.cdna.all.fa



# TraDetIONS
Transposon detection in Oxford Nanopore Sequencing data

Detects TE from nanopore sequencing

Use the following scripts for these functions
1. survivor_retro_annotate.py    recognizes TE sequence
2. mark_variants.py              recognizes potential somatic insertions
3. merge_insertoins.py            merges insertions
4. get_insertion_reads.py         gets reads supporting insertion
5. rmdups.sh                      removes duplicate reads 
6. create_fasta_with_insertion_seq.sh  polishes insertion with insertion sequence
7. create_fasta_with_read.sh  polishes insertion with representative read
8. insertion_annotation.py    annotates insertion

9. pseudogene_annotation.py    annotates pseudogenes

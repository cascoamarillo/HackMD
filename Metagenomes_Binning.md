# Concoct Binning
###### tags: `metagenomic`, `binning`

https://concoct.readthedocs.io/en/latest/usage.html

First get the alignment file (bam) from assembly and fastq reads with bwa
Sam file:
`bwa mem -t $threads assembly.fa reads_1 reads_2`

Sam to Bam:
`samtools view -@ 10 -bS NH_P-idba_contig_EukRep.sam > NH_P-idba_contig_EukRep.bam`

Index bam:
`samtools sort -@ 12 NH_P-idba_contig_EukRep.bam NH_P-idba_contig_EukRep.sort`

## concoct

Cut contigs into smaller parts:
`~/miniconda3/envs/concoct_env/bin/cut_up_fasta.py ../NH_P-idba_contig_EukRep.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa`

Generate table with coverage depth information per sample and subcontig:
`~/miniconda3/envs/concoct_env/bin/concoct_coverage_table.py contigs_10K.bed ../metawrap_binning/NH_P-idba_contig_EukRep.sort.bam > coverage_table.tsv`

Run concoct:
`~/miniconda3/envs/concoct_env/bin/concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/ -t 10`


Merge subcontig clustering into original contig clustering:
`~/miniconda3/envs/concoct_env/bin/merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv`

Extract bins as individual FASTA:
`mkdir concoct_output/fasta_bins`

`~/miniconda3/envs/concoct_env/bin/extract_fasta_bins.py ../NH_P-idba_contig_EukRep.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins`




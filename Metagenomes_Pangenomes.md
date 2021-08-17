# Pangenomes
###### tags: `metagenomes`, `pangenomes`


An anvi'o workflow for microbial pangenomics

https://merenlab.org/2016/11/08/pangenomics-v2/


##  anvi-gen-contigs-database

Generate the contigs (assembly) database for anvio:

https://merenlab.org/software/anvio/help/7/programs/anvi-gen-contigs-database/

Working directory:
`/groups/fr2020/metagenomes/Permafrost_RUS/Pangenome/`

First, remove extra headers in the fasta file for anvio

`sed '/^>/ s/ .*//' /groups/fr2020/metagenomes/Permafrost_RUS/PRJNA72925/active_layer/idba_ud_out/Active-layer_contig_Rotifera-GC50.fasta > Perma925_active.fasta`

USE the assembly after Rotifera and GC 50 filtering.

Then:
`anvi-gen-contigs-database -f Perma925_active.fasta -T 4 -o Perma925_active.db`

## anvi-run-hmms

Populate the contig database with HMM (Protein) hits:

https://merenlab.org/software/anvio/help/7/programs/anvi-run-hmms/

`anvi-run-hmms -c Perma925_active.db -H /groups/fr2020/metagenomes/Permafrost_RUS/Pangenome/Rotifera_4340 -T 6 --also-scan-trnas`

We are using a set of 4340 pfam protein models associated with Rotifera. Steps for this:
https://merenlab.org/software/anvio/help/7/artifacts/hmm-source/

-download the list of pfams from http://pfam.xfam.org/ncbiseq/398365647#tabview=tab4

Put it on a list, like in ~/bin/pfam/Rotifera_Pfam-list.txt:

```
PF00001
PF00002
PF00003
PF00004
PF00005
```

 And create anvi’o HMM sources from ad hoc PFAM accessions

`anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-file ~/bin/pfam/Rotifera_Pfam-list.txt -O ~/bin/pfam/Rotifera_4340/`

Make sure it has the same structure as in https://github.com/merenlab/anvio/tree/master/anvio/data/hmm

like change inside `kind.txt` to `singlecopy:eukarya`

## Create a file "external-genomes.txt":

```
name	contigs_db_path
Perma925_active	Perma925_active.db
Perma925_permafrost	Perma925_permafrost.db
Perma250_827	Perma250_827.db
```

## Create an anvio Storage Genomes

`anvi-gen-genomes-storage -e external-genomes.txt -o PERMAFROST-GENOMES.db`

## Running the pangenome analysis:

`anvi-pan-genome -g PERMAFROST-GENOMES.db --project-name "Permafrost_Rotifera" --output-dir PERMAFROST --num-threads 12 --minbit 0.5 --mcl-inflation 10 --use-ncbi-blast`

## Displaying the pan genome

`anvi-display-pan -g PERMAFROST-GENOMES.db -p PERMAFROST/Permafrost_Rotifera-PAN.db`

Then, in another terminal we need to open a ssh tunnel:
`ssh -L8080:localhost:8080 frodriguez@minnie`

Open Chrome and go to (http://localhost:8080/)

## Exporting graph
`inkscape -z -e file.png -w 1024 file.svg`
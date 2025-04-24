# Genome decontamination (*Anvi'o*)
###### tags: `Assembly`
## Workflow
1. Map raw reads to raw assembly
	*  `bowtie` or `minimap2`
	*  Input: Raw fasta & assembly fasta
	*  Output: Alignment `.sam` file
2. Convert alignment `.sam` file to `.bam` format
	* `samtools`
3. Inspect ONT mapping report
	* `samtools`
4. Annotate assembly taxonomies
	* `blastn`
5. Initially assess contamination
	* `blobtools`
	* Input: Alignment `.bam` files, raw assembly fasta & raw assembly taxonomic annotations
6. Create contigs database (CDB)
	* `anvi'o`
	* Holds basic information assembly scaffolds (k-mer frequency, GC-content, ...)
	* Additional steps for:
		* COG gene annotations
		* `centrifuge` taxonomy annotations
7. Profile `.bam` file to CDB
	* `anvi’o` Profiles describes statistics per scaffold in `.bam` file (average coverage, portion of each scaffold covered by at least one read, ...)
8. Merge `anvi’o` profiles
	* Computes hierarchical clustering of scaffolds
9. Visualize merged data on `anvi’o`'s interface
	* Allows identification of draft genome bins & removal of contaminants

## MBL Cluster Preparation
File paths:
*  Home directory (`~`): `/groups/fr2020/Isa/`
*  ONT (to decontaminate):
	*  Raw reads:
`~/ONT_raw_reads/ONT_2019_2021_all.fastq.gz`
	*  Assembly:
`~/pbjelly/Ds_jelly-assembly.fasta`
* Illumina (if needed):
	* Reads:
`~/groups/fr2020/Isa/Illumina_ds_paired_reads/`
	* Assembly:
`~/Ds_illumina_genome/Dst_b1v03.fasta`
	* Annotations:
`~/Ds_illumina_genome/Dst_b1v03.max_arth_b2g_droso_b2g_ctg-names.gff`

# [BlobTools](https://blobtoolkit.genomehubs.org/blobtools2/)
## Preparing `blobtools`
**Goal:** Re-run `blobtools` with lastest NCBI `nt` database (October 2025), and compare with coverage files.

### 1. Map reads to the assembly
#### Illumina
Mapping the Illumina reads to the `pbjelly` assembly with [`bowtie`](https://blobtoolkit.genomehubs.org/blobtools2/):
```
cd /groups/fr2020/Isa/

bowtie2 -p 12 -x Ds_jelly-assembly \
-1 Illumina_ds_paired_reads/trimmomatic/Darwinula_stevensoni_250_350_pass_paired_1_fixed.fastq.gz \
-2 Illumina_ds_paired_reads/trimmomatic/Darwinula_stevensoni_250_350_pass_paired_2_fixed.fastq.gz \
-S pbjelly/coverage/Ds_jelly-assembly_250_350.sam
```

#### ONT
Mapping the ONT reads to the `pbjelly` assembly with `minimap2`:
```
module load minimap2
cd pbjelly

minimap2 -ax map-ont Ds_jelly-assembly.mm ../ONT_raw_reads/ONT_2019_2021_all.fastq.gz -t 42 > coverage/Ds_jelly-assembly_ONT.sam
```

### 2. Convert alignment `.sam` to `.bam`
Convert SAM to BAM format:
```
cd coverage

samtools view Ds_jelly-assembly_250_350.sam -S -b > Ds_jelly-assembly_250_350.bam
samtools sort Ds_jelly-assembly_250_350.bam > Ds_jelly-assembly_250_350.sort.bam

samtools view Ds_jelly-assembly_ONT.sam -S -b > Ds_jelly-assembly_ONT.bam
samtools sort Ds_jelly-assembly_ONT.bam > Ds_jelly-assembly_ONT.sort.bam
```

### 3. Inspect ONT mapping report
View the mapping report of ONT reads (important for later):
```
samtools flagstat Ds_jelly-assembly_ONT.sort.bam

3588820 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
3385697 + 0 mapped (94.34%:-nan%)
```

### 4. Annotate taxonomies with `blastn`
Blast the raw assembly against NCBI's most recent *nt* database:
```
blastn -query ../Ds_jelly-assembly.fasta \
-db nt  \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-evalue 1e-25 \
-num_threads 20 \
-out ../blobtools/Ds_jelly-assembly-blastn.out
```
:::info
* `-outfmt '6 qseqid staxids bitscore std'`:
	* [Output format six](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) with 3 specified columns, and then `std`
	* `std`: Standard columns (n = 11)
	* Resulting table has 11 columns
		* `qseqid` Is both manually specified, and included in `std`
			* Therefore, it is included twice
:::
:::warning
The BLAST database on the MBL cluster is outdated. Therefore, this was done using RIT BLAST databases.
:::
### 5. Initially assess contamination via [`blobtools`](https://blobtoolkit.genomehubs.org/blobtools2/)
:::info
I have [`blobtools`](https://blobtoolkit.genomehubs.org/blobtools2/) installed in my home directory; it's not in bioware.
:::
#### Illumina
Run [`blobtools`](https://blobtoolkit.genomehubs.org/blobtools2/):

```
~/bin/blobtools/blobtools create \
-i ../Ds_jelly-assembly.fasta \
-b ../coverage/Ds_jelly-assembly_250_350.sort.bam \
-t Ds_jelly-assembly-blastn.out \
-o Ds_jelly-assembly_Ill-blastn-blobplot

~/bin/blobtools/blobtools view \
-i Ds_jelly-assembly-blastn.blobDB.json

~/bin/blobtools/blobtools plot -i Ds_jelly-assembly-blastn.blobDB.json
```

#### ONT
Run [`blobtools`](https://blobtoolkit.genomehubs.org/blobtools2/):
```
~/bin/blobtools/blobtools create \
-i ../Ds_jelly-assembly.fasta \
-b ../coverage/Ds_jelly-assembly_ONT.sort.bam \
-t Ds_jelly-assembly-blastn.out \
-o Ds_jelly-assembly-blastn_ONT-blobplot

~/bin/blobtools/blobtools view \
-i Ds_jelly-assembly-blastn.blobDB.json

~/bin/blobtools/blobtools plot -i Ds_jelly-assembly-blastn.blobDB.json
```
:::warning
**NOT SURE if we can combine with in one plot report (!!!)**
:::

**EVERYTHING** in: `/groups/fr2020/Isa/pbjelly/blobtools/`

Compare old report (old BLAST DB) in `/groups/fr2020/Isa/pbjelly/blobtools/oldDB/` with new ones:

Illumina:
`Ds_pbjelly.bestsum.Illumina.png`
`Ds_pbjelly.bestsum.Illumina.read_cov.png`

ONT:
`Ds_pbjelly.bestsum.ONT.png`
`Ds_pbjelly.bestsum.ONT.read_cov.png`

:::info
**++Result interpretation++**
* Illumina:
	* Little (negligible) contamination
	* 8.20 % did not map
* Illumina + old NCBI `nt` database:
	* Suggested 11.65% of mapped reads **(11.64% increase)** came from conaminants
	* (Did not affect mapping)
* ONT + latest NCBI `nt` database:
	* Some contamination -> Requires cleanup
		* Mostly Proteobacteria
	* 49.85 % did not map -> Reverse reads? No idea (FR). But I get (samtools flagstat) 94.34% of ONT reads mapped to the Ds pbjellt assembly with minimap2. **We need to investigate this.**
:::

# ANVIO
### 6. Create contigs database
Create CDB:
```
module load anvio/8-conda
cd ../anvio

anvi-gen-contigs-database -f ../Ds_jelly-assembly.fasta \
-o contigs.db -T 24 -n Ds_jelly -g 1
```
:::info
* `-f`: The FASTA file that contains reference sequences you mapped your samples against.
* `-n`: Name of the project.
* `-g`: Translation table to use
:::
`anvio 8` Output clarifies used parameters:
```
##uses k-mer = 4
##--split-length 20000 (default)
##--prodigal-translation-table 1 (in PRODIGAL -g 1)

### CONTIGS DB CREATE REPORT
===============================================
Split Length .................................: 20,000
K-mer size ...................................: 4
Skip gene calling? ...........................: False
External gene calls provided? ................: False
Ignoring internal stop codons? ...............: False
Splitting pays attention to gene calls? ......: True
Contigs with at least one gene call ..........: 1235 of 1243 (99.4%)
Contigs database .............................: A new database, contigs.db, has been created.
Number of contigs ............................: 1,243
Number of splits .............................: 20,912
Total number of nucleotides ..................: 427,252,629
Gene calling step skipped ....................: False
Splits broke genes (non-mindful mode) ........: False
Desired split length (what the user wanted) ..: 20,000
Average split length (what anvi'o gave back) .: 20,433
```

"Populate `contigs.db` with Hidden Markov Model hits":
```
anvi-run-hmms -c contigs.db --num-threads 20
```

#### Additional CDB annotations
Protocol reference can be found [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC4824900/)

##### Annotation via COGs
Annotate genes in `contigs-db` with functions using NCBI’s Clusters of Orthologus Groups (COGs) database:
```
anvi-run-ncbi-cogs -c contigs.db \
--cog-data-dir /blastdb/cogs -T 24

anvi-get-dna-sequences-for-gene-calls -c contigs.db \
-o gene-calls.fa
```

*"How many genes did Anvi'o find on your contigs?"*:
```
egrep ">" gene-calls.fa | wc -l
```
```
468511
```

##### Annotation via `centrifuge`
Protocol reference can be found [here](https://merenlab.org/2016/06/18/importing-taxonomy/).
Classify microbes in metagenomic sequences using `centrifuge`:
```
export PATH=/bioware/centrifuge-1.0.3-beta:$PATH

export CENTRIFUGE_BASE="/groups/fr2020/Isa/pbjelly/anvio/CENTRIFUGE"
cd $CENTRIFUGE_BASE

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p+h+v.tar.gz
tar -zxvf p+h+v.tar.gz && rm -rf p+h+v.tar.g

centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v gene-calls.fa \
-S centrifuge_hits.tsv
```

##### Import annotations to CDB
Import taxonomies into `anvi'o`:
```
anvi-import-taxonomy -c contigs.db \
-i centrifuge_report.tsv centrifuge_hits.tsv \
-p centrifuge
```

### 7. Profile `.bam` to CDB
Profile `.bam` files relative to `anvi'o` CDB:
```
anvi-profile -i ../coverage/Ds_jelly-assembly_250_350.sort.bam \
--output-dir ./illumina_profile \
--sample-name illumina_250_350PE \
-c contigs.db --min-contig-length 800 --num-threads 24

anvi-profile -i ../coverage/Ds_jelly-assembly_ONT.sort.bam \
--output-dir ./ONT_profile \
--sample-name ONT \
-c contigs.db --min-contig-length 800 --num-threads 24
```

#### 8. Merge profiles
Merge the different `anvi'o` profiles:
```
anvi-merge illumina_profile/PROFILE.db ONT_profile/PROFILE.db \
-o Ds_illumina-ONT_merge -c contigs.db \
-S Ds_pbjelly --enforce-hierarchical-clustering
```

Interactive (fun!) mode:
```
anvi-interactive -p Ds_illumina-ONT_merge/PROFILE.db \
-c contigs.db --server-only
```
:::info
* For use of `anvi-interactive` from a PC outside the MBL network, connect to the servers using Windows Powershell and the following command:
```
ssh -L8123:localhost:8080 -J ${USERNAME}@evol5.mbl.edu ${USERNAME}@minnie
```
:::

:::info
++**Report is found at:**++
`/groups/fr2020/Isa/pbjelly/anvio/bacteria/Ds_illumina-ONT_merge/SUMMARY_default/`
:::

# NOW WHAT?
Goods news: there are not many contaminants ~16 Megabases (~250 contigs)

Bad news: we have to remove those.

Anvio gives you the Fasta file with the Two Bins (bacteria vs Ds) in the reprot (bin_by_bin). But I think it would be ideal to combine both reports (blobtools and anvio), like a big table with contigs and taxonomic assignation and extract our filtered Darwinula ONLY contigs for future analysis. :) 

### Combine `blobtools` and `anvi'o` reports
**Goal:** One table, listing all contigs and their taxonomies (see below)

|Contig|tax_blobtools_Ill|tax_blobtools_ONT|bin_anvio|
|-|-|-|-|
| Contig 0     | Arthropoda     |  Arthropoda   |Ds|
|...|...|...|...|
|Contig1242|Proteobacteria|Proteobacteria|Bacteria|

Input files:
* `blobtools/Ds_jelly-assembly-blastn_ONT.Ds_jelly-assembly-blastn_ONT-blobplot.blobDB.table.txt`
* `blobtools/Ds_jelly-assembly-blastn-blobplot.blobDB.table.txt`
* `anvio/bacteria/Ds_illumina-ONT_merge/SUMMARY_default/bin_by_bin/Bin_bacteria/Bin_bacteria-contigs.fa`
* `anvio/bacteria/Ds_illumina-ONT_merge/SUMMARY_default/bin_by_bin/Bin_Ds/Bin_Ds-contigs.fa`
:::warning
**What about `anvio/Ds_illumina-ONT_merge/SUMMARY_default/bin_by_bin/...`?**
:::

Processing:
* `blobtools` files:
	1. Remove `#`-containing lines (`grep`)
	2. Extract columns of interest (`awk`) as `blobtools_contigs`
```
grep -v "#" blobtools/Ds_jelly-assembly-blastn-blobplot.blobDB.table.txt \
| awk '{print $1, $6}' \
> blobtools/blobtools_comb-report-prep_Ill.txt

grep -v "#" blobtools/Ds_jelly-assembly-blastn_ONT.Ds_jelly-assembly-blastn_ONT-blobplot.blobDB.table.txt \
| awk '{print $6}' \ 
> blobtools/blobtools_comb-report-prep_ONT.txt

paste blobtools/blobtools_comb-report-prep_Ill.txt blobtools/blobtools_comb-report-prep_ONT.txt | awk '{print $1, $2, $4}' | sort > blobtools/blobtools_comb-report-prep_ONT-Ill.txt
```
:::warning
`Contig0` misses from `anvio_comb-report-prep.txt`
:::
* `anvi'o` files:
	1. Extract headers (`grep`)
	2. Remove `>` (`sed`)
	3. Export lists of contigs as `anvio_contigs_Ds` and  `anvio_contigs_Bacteria`
```
grep ">" anvio/bacteria/Ds_illumina-ONT_merge/SUMMARY_default/bin_by_bin/Bin_Ds/Bin_Ds-contigs.fa | sed 's/>//g' | awk '{print $1, "Ds"}' > anvio/anvio_contigs_Ds.txt

grep ">" anvio/bacteria/Ds_illumina-ONT_merge/SUMMARY_default/bin_by_bin/Bin_bacteria/Bin_bacteria-contigs.fa | sed 's/>//g' | awk '{print $1, "Bacteria"}'> anvio/anvio_contigs_Bacteria.txt

cat anvio/anvio_contigs_Ds.txt anvio/anvio_contigs_Bacteria.txt | sort > anvio/anvio_comb-report-prep.txt

join blobtools/blobtools_comb-report-prep_ONT-Ill.txt anvio/anvio_comb-report-prep.txt | awk 'BEGIN{print "contig", "blob_Illumina", "blob_ONT", "anvio"}1' > blob_anvio_report.txt
```
	5. Create a new, two-column file in a for-loop (see table below) as `anvio_contigs`:
		* This file lists the contig (in the 1^st^ column), and the `anvi'o` taxonomy (in the 2^nd^column, added through a conditional statement comparing to the `anvio_contigs_...` files)
```
touch anvio_comb-report-prep_Ds_Bacteria.txt
```
	6. Horizontally concatenate `blobtools_contigs` and `anvio_contigs` (`paste`)
	7. Add header line (`cat`)

Table from step 4:
|contig|bin_anvio|
|-|-|
|Contig1| Bacteria|
|...|...|
|Contig1242| Ds|

### Extract *Darwinula* contigs
1. Unzip reads:
```
gunzip /groups/fr2020/Isa/ONT_raw_reads/ONT_2019_2021_all.fastq.gz
```
2. Filter summary file to *D. stevensoni* contigs (`awk`)
3. Extract *D. stevensoni* contigs by name  (`asub`)

### Assemble *Cardinium* genome
1. Filter summary file to *Cardinium* contigs (`awk`)
	* Filter to Proteobacteria, or something more specific?
2. Extract *Cardinium* contigs by name from the summary file  (`asub`)
3. Map raw ONT reads to extracted *Cardinium* contigs (`minimap2`)
4. Extract raw ONT reads which mapped to *Cardinium* contigs from `.sam` file (`awk`)
5. Assemble *Cardinium* genome (`flye`)
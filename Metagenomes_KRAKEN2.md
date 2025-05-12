# KRAKEN2
###### tags: `metagenomic`, `KRAKEN`

A nice tutorial on taxonomic investigation:
https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html
##  Run kraken2
Manual of kraken2 https://ccb.jhu.edu/software/kraken2/
`module load kraken2`

`kraken2  --threads 10 --report assembly.txt assembly.fa > assembly.kraken`

this will use the standard DB loaded with kraken "/blastdb/kraken2/kraken2-db-20200325", which is good for bacteria, virus and human sequences.


Specific database:
`kraken2 --db /groups/arklab/KRAKEN2/minikraken_8GB_20200312 --threads 10 --report assembly.txt assembly.fa > assembly.kraken`

this DB is a small representation of bacteria/virus sequences

This classification may take a while, depending on how many sequences we are going to classify. The resulting content of the file evol1.kraken looks similar to the following example:

```
C       7001326F:121:CBVVLANXX:1:1105:2240:12640        816     251     816:9 171549:5 816:5 171549:3 2:2 816:5 171549:4 816:34 171549:8 816:4 171549:2 816:10 A:35 816:10 171549:2 816:4 171549:8 816:34 171549:4 816:5 2:2 171549:3 816:5 171549:5 816:9
C       7001326F:121:CBVVLANXX:1:1105:3487:12536        1339337 202     1339337:67 A:35 1339337:66
U       7001326F:121:CBVVLANXX:1:1105:5188:12504        0       251     0:91 A:35 0:91
U       7001326F:121:CBVVLANXX:1:1105:11030:12689       0       251     0:91 A:35 0:91
U       7001326F:121:CBVVLANXX:1:1105:7157:12806        0       206     0:69 A:35 0:68
```

Each sequence classified by Kraken2 results in a single line of output. Output lines contain five tab-delimited fields; from left to right, they are:

    C/U: one letter code indicating that the sequence was either classified or unclassified.

    The sequence ID, obtained from the FASTA/FASTQ header.

    The taxonomy ID Kraken2 used to label the sequence; this is 0 if the sequence is unclassified and otherwise should be the NCBI Taxonomy identifier.

    The length of the sequence in bp.

    A space-delimited list indicating the lowest common ancestor (in the taxonomic tree) mapping of each k-mer in the sequence. For example, 562:13 561:4 A:31 0:1 562:3 would indicate that:

        the first 13 k-mers mapped to taxonomy ID #562

        the next 4 k-mers mapped to taxonomy ID #561

        the next 31 k-mers contained an ambiguous nucleotide

        the next k-mer was not in the database

        the last 3 k-mers mapped to taxonomy ID #562


Rotifera database:

`kraken2 --db /groups/arklab/KRAKEN2/Rotifera --threads 10 assembly.fa > SCN_BF_Rotifera_DB.kraken`


#  Investigate taxa

We can use the webpage NCBI TaxIdentifier to quickly get the names to the taxonomy identifier. However, this is impractical as we are dealing potentially with many sequences. Kraken2 has some scripts that help us understand our results better.

Because we used the Kraken2 switch --report FILE, we have got also a sample-wide report of all taxa found. This is much better to get an overview what was found.

The first few lines of an example report are shown below.

```
83.56  514312  514312  U       0       unclassified
16.44  101180  0       R       1       root
16.44  101180  0       R1      131567    cellular organisms
16.44  101180  2775    D       2           Bacteria
13.99  86114   1       D1      1783270       FCB group
13.99  86112   0       D2      68336           Bacteroidetes/Chlorobi group
13.99  86103   8       P       976               Bacteroidetes
13.94  85798   2       C       200643              Bacteroidia
13.94  85789   19      O       171549                Bacteroidales
13.87  85392   0       F       815                     Bacteroidaceae
```

The output of kraken-report is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

    Percentage of reads covered by the clade rooted at this taxon

    Number of reads covered by the clade rooted at this taxon

    Number of reads assigned directly to this taxon

    A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply “-“.

    NCBI Taxonomy ID

    The indented scientific name

Now, we use the tool ktImportTaxonomy from the Krona tools to create the html web-page. We first need build a two column file (read_id<tab>tax_id) as input to the ktImportTaxonomy tool. We will do this by cutting the columns out of the Kraken2 file:

```
cd kraken dir
cat assembly.kraken | cut -f 2,3 > assembly.kraken.krona
ktImportTaxonomy assembly.kraken.krona -o assembly.kraken.krona.html
firefox taxonomy.krona.html
```

# Build custom DB

Download assembly/genome fasta files

Rename headers
`>Sp_ctg_1`

Get tax ID from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

Add tax ID to headers:

`~/bin/bioawk/bioawk -c fastx '{print ">"$name"|kraken:taxid|104782\n"$seq}' Avaga.fa >Avaga_tax104782.fa`


Concatenate fasta files
cat *.fa > Rotifera_tax.fa

`kraken2-build --threads 12 --no-masking --add-to-library Rotifera_tax.fa --db Rotifera`

Also, download the taxonomy information:

    kraken2-build --download-taxonomy --db $DBNAME

Finally, build the database

    kraken2-build --build --db $DBNAME
    
### Protein database:
`kraken2-build --protein --threads 12 --no-masking --add-to-library /groups/fr2020/RepeatPeps.lib_taxid.fa --db RepeatPeps`

## UPDATE 2021

We can run KRAKEN with reads ONLY too:
Run kranken2 on merge fastq files:

Bioproject PRJNA680161 (Permafrost)


`module load kraken2`

Decompress .fastq.gz files:

`gunzip *.fastq.gz`

Concatenate the libraries *_1_paired.fastq & *_2_paired.fastq:

`cat 1*_1_paired.fastq 2*_1_paired.fastq 3*_1_paired.fastq > _1_paired.fastq`

Do the same for *_2_paired

Run kraken2 with --paired
`
kraken2 --paired --db /groups/arklab/KRAKEN2/Rotifera --threads 12 *_1_paired.fastq *_2_paired.fastq > PRJNA680161_Rotifera.kraken`

When finished:

`cat *_Rotifera.kraken | cut -f 2,3 > *_Rotifera.kraken.krona`

Import module KRONA:
    `module load kronatools/`
    
`ktImportTaxonomy *_Rotifera.kraken.krona -o PRJNA680161_Rotifera.kraken.krona.html`

> Now (May 2025), it seems the taxonomy file should be directed in the prompt:
    `ktImportTaxonomy *_paired_Rotifera.kraken.krona -o *_paired_Rotifera.kraken.krona.html -tax /blastdb/kronatools-2.7/`
    
`firefox taxonomy.krona.html`
or simply open it with Filezilla and your browser.

## Extract sequences (assembly) matching a specific database.....

Let's say we want to extract the sequences (contigs) that somehow match the Rotifera_DB with kraken. We need to generate the kraken.krona.html and select on the group of interest (Rotifera). Then, on the top right corner it says "count" with the number of sequences belonging to that group. If we click on that number (596 on the pic) it will show a list with the contig IDs.

![](https://i.imgur.com/60rXVRp.png)

We copy that list and paste it in a text file in the directory (you can use filezilla for this).
The text file with the list will be something like this:

$ more tabulated-IDs.txt
ID_1
ID_2
ID_3
...
ID_N

then:

`seqkit grep --pattern-file Tabulated-IDs.txt Input-fasta.fa > IDs-fasta.fa`

Count number of sequences inside a fasta file:

`grep -c "^>" Perma827_Rotifera.fa`

If we have a large portion of reads with high GC content, meaning bacteria, we can get rid of the by:

`~/bin/prinseq-lite-0.20.4/prinseq-lite.pl -max_gc 50 -fasta input.fa`

Will output two files: *_good (with the seqs within the limit, ie. less than 50%GC) and *_bad the rest of seqs.


### Extract sequences matching a specific database...another way:

After running Kraken, Kraken2, or KrakenUniq, users may use the extract_kraken_reads.py program to extract the FASTA or FASTQ reads classified as a specific taxonomy ID. 
https://ccb.jhu.edu/software/krakentools/index.shtml?t=extractreads

`extract_kraken_reads.py`

In:
```
/groups/fr2020/bin/KrakenTools/
```

```
usage: extract_kraken_reads.py [-h] -k KRAKEN_FILE -s SEQ_FILE1
                               [-s2 SEQ_FILE2] -t TAXID [TAXID ...] -o
                               OUTPUT_FILE [-o2 OUTPUT_FILE2] [--append]
                               [--noappend] [--max MAX_READS] [-r REPORT_FILE]
                               [--include-parents] [--include-children]
                               [--exclude] [--fastq-output]
extract_kraken_reads.py: error: the following arguments are required: -k, -s/-s1/-1/-U, -t/--taxid, -o/--output
```

# TE assignation with kraken

New Rotifera TE database with RepBase classification is in
 `/groups/arklab/KRAKEN2/RotiferaTE_Repbase`
 
 So you can run:
 
 `kraken2 --db /groups/arklab/KRAKEN2/RotiferaTE_Repbase --threads 12 Assembly.fa > Assembly_RotiferaTE_Repbase.kraken`
 
 or
 
 `kraken2 --paired --db PASSED_R1.fastq PASSED_R2.fastq > Passes_RotiferaTE_Repbase.kraken`
 
 Then
 
 `cat *RotiferaTE_Repbase.kraken | cut -f 2,3 > *RotiferaTE_Repbase.kraken.krona`
 
 And
 
 `ktImportTaxonomy -tax /groups/arklab/KRAKEN2/kronatools/taxonomy/ *RotiferaTE_Repbase.kraken.krona -o *RotiferaTE_Repbase.kraken.krona.html`
 
 **IMPORTANT**: use `-tax /groups/arklab/KRAKEN2/kronatools/taxonomy/` it is where the taxonomy.tab file contains the hierarchy of TEs to show.
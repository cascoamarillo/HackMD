# TE annotation

###### tags: `TE`, `annotation`, `Repeat Masker`


Once we have a nice valid TE custom library of our new species, we can take a closer look at how it is distributed (across contigs/cosmids) in the genome. Along with the proportions (%), a way to look at this is generating annotation files like [.gff](https://m.ensembl.org/info/website/upload/gff3.html)

```
##gff-version 3
ctg123 . mRNA            1300  9000  .  +  .  ID=mrna0001;Name=sonichedgehog
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mrna0001
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mrna0001
ctg123 . exon            3000  3902  .  +  .  ID=exon00003;Parent=mrna0001
ctg123 . exon            5000  5500  .  +  .  ID=exon00004;Parent=mrna0001
ctg123 . exon            7000  9000  .  +  .  ID=exon00005;Parent=mrna0001
```

Of course there are versions and another variants (.bed, genebank, etc). I like gff3 btw.

Repeat Masker allows you to, first map you TElib in your genome and export this as .gff file format. We need to specify the `-gff` option. Like:

`RepeatMasker -dir RM -pa 6 -a -low -no_is -e ncbi -gff -lib /groups/arklab/Repbase/invrep25.12.ref.fasta contig.fa`

`-dir RM` is to better organize the output in folders

I would like to do this for:
1. The custom TE library
2. The RepeatMasker/Repbase library

Custom TE sould be `Dscosmid_denovoLibTEs_filtered_MCL.fa.classified`, cos the extra #TEtype in each fasta header will by use by RepeatMasker
```
>RIX-comp_MCL1_Dscosmid-B-G12-Map3#LINE/CR1 
AGGCTTTGAGTCGCTGTAGAAGTCCGATGGCTAGGGATGGAGTACATTTT
TCGAAAGAGGGAGCGGTAGCAGTAGGATTGGCTATCATGAAGGAAAATGA
GCCTTTTTTAGGATTGTAGATGGGGGGGCAAGGGATCACAAGCAAGGGAC
AGCCACTACGAGCAGCTTTAGACCGCCGCGTACTGGTGTGCAGGAAAGGC
```
The Repbase is just not using -lib option. It will pick `/groups/arklab/bin/RepeatMasker/Libraries/RepeatMasker.lib`

Let's do first the RM/Repbase lib (no -lib option). Just type:

`RepeatMasker -dir RM -pa 6 -a -low -no_is -e ncbi -gff cosmid.fa`

## Repeat Protein masker


Uses a protein TE library:
`/groups/arklab/bin/RepeatMasker/Libraries`

ALWAYS symlink the input fasta file (assembly.fa) to the working dir; otherways it's going to place output files where the fasta file is located!

`RepeatProteinMask -engine ncbi -pvalue 0.001 -noLowSimple contig.fa`

This will create two files:
.annot: table with TE annotation
.masked: input fasta file masked (regions with TE are N's). Not useful for us.
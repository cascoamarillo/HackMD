# Fasta Input manipulation

## GC content
To check GC% distribution in any given fasta file (assembly) we can run the script:

`/groups/fr2020/bin/average_GC_hist_FASTA.sh` as (in your working dir with the fasta file):

`/groups/fr2020/bin/average_GC_hist_FASTA.sh infile.fa 1`

I will create an output GC.png with an histogram. Always, rename the output for further details "infile_CG.png"
*The number "1" dictates how to bin/group the data...not really important now.



## Fasta/Assembly Summary

We are using the stats script from BBTools
https://jgi.doe.gov/data-and-tools/bbtools/

In the folder with your fasta file:

`/bioware/bbtools-38.84/stats.sh in=contig.fa > contig.summary.txt`


## Extract fasta seqs by IDs:

####Extract IDs from assembly:
You nee a text file (pattern file) with a list of IDs (>ID) tabulated:

```
$ more abulated-IDs.txt
ID_1
ID_2
ID_3
...
ID_N
```
then:

`seqkit grep --pattern-file Tabulated-IDs.txt  Input-fasta.fa > IDs-fasta.fa`

## Count number of sequences inside a fasta file:

```
grep -c "^>" Perma827_Rotifera.fa
```

## More about seqkit tool in
https://bioinf.shenwei.me/seqkit/


## Split fasta file based on GC%

`~/bin/prinseq-lite-0.20.4/prinseq-lite.pl -max_gc 50 -fasta input.fa`

Will output two files: *_good (with the seqs within the limit) and *_bad the rest of seqs.

## Fasta size selection

First, let check the size distribution of a given multifasta file. like:

`cat file.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'`

It will output in the terminal, or you can pipe it to a file `> file-sizes.txt`

Then select by size with seqkit tool:

fitlering sequences with minimum 100 aa/nt length i.e > 100,

`seqkit seq -m 100 test.fa`

For filtering sequences with maximum of 1000 aa/nt length<=1000

 `seqkit seq -M 1000 test.fa`

For filtering sequences between maxiumum size 200 and minimum size of 100:

 `seqkit seq -m 100 -M 200 test.fa`



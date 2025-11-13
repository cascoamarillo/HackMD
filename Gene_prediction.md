# Gene prediction
###### tags: `gene`, `prediction`, `annotation`


Since we have a small genome (cosmids), we are going to run Augustus web server:

http://bioinf.uni-greifswald.de/webaugustus/

We pick the first option: AUGUSTUS training submission"

## AUGUSTUS training submission"

Fill all the fields

![](https://i.imgur.com/Y2Jgkku.png)

-Your email

-Upload the cosmid fasta file

-They ask for one of this protein-related files....we are going to cheat a little: we are using the protein annotation from Illumina in ostracods:
`/groups/arklab/ostracods/illumina/Dst_b1v03.prot.fa`

-Yep, we are not working with human data.

Everything should be OK. They send the results to the email...check it for spam

Then you need to save the .tar.gz archive 

## Augustus prediction submission

![](https://i.imgur.com/GkJP1xE.png)

-Upload AUGUSTUS species parameters (.tar.gz archive)

-Genome file: cosmids fasta file. BUT:

You are going to do this twice:
    -One with Dscosmid.fa
    -another with masked file: Dscosmid.fa.masked in `/groups/arklab/ostracods/cosmids/RedoTEdenovo/divergencePlot_N3`

-No need to upload cDNA file or Hints file.

-No UTR prediction

![](https://i.imgur.com/bcPvtlT.png)


Wait for results


# Gene prediction: BRAKER2
**Used for A. vaga L1 strain (hybrid assembly with PBJelly)**

```
braker.pl --species=As11 --genome=genome.fasta \
--bam=file1.bam,file2.bam \ --bam=/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L1/Av11_RNA_L1.sort.bam,/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L2/Av11_RNA_L2.sort.bam
--prot_seq=/users/frodriguez/Avgenome/Genoscope/Adineta_vaga_v2.0.annot.pep_edit.fa \
--prg=(gth|exonerate|spaln)
--cores=
--gff3
--UTR=on
--softmasking
--AUGUSTUS_CONFIG_PATH=/users/frodriguez/bin/Augustus/config
--GENEMARK_PATH=/bioware/genemark-et-4.32/gmes_petap
```


### Some scripts to be included the env
```
export PATH="/users/frodriguez/bin/Augustus/scripts:$PATH"
export PATH="/users/frodriguez/bin/Augustus/bin:$PATH"
export BAMTOOLS ="/bioware/bamtools-2.4.0/bin:$PATH"
```

### For gth aligner
`export PATH="/users/frodriguez/bin/gth-1.7.1-Linux_x86_64-64bit/bin:$PATH"`


### Augustus
```
export AUGUSTUS_CONFIG_PATH=/users/frodriguez/bin/Augustus/config
export AUGUSTUS_BIN_PATH=/users/frodriguez/bin/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/users/frodriguez/bin/Augustus/scripts
```


### Some additional things to be included
```
setenv $BSSMDIR         "/users/frodriguez/bin/gth-1.7.1-Linux_x86_64-64bit/bin/bssm"
setenv $GTHDATADIR      "/users/frodriguez/bin/gth-1.7.1-Linux_x86_64-64bit/bin/gthdata"
setenv $GENEMARK_PATH      "/bioware/genemark_suite_linux_64-20160428"
```

## Different BRAKER2 runs to check parameter/options/outputs 
### Simple run to see if it works

`braker.pl --species=As11 --genome=/users/frodriguez/Avgenome/Av11-2/PBJelly/te/RM/softmasking/Av11_pb_v2.soft.mask.fa --bam=/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L1/Av11_RNA_L1.sort.bam,/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L2/Av11_RNA_L2.sort.bam --prot_seq=/users/frodriguez/Avgenome/Genoscope/Adineta_vaga_v2.0.annot.pep_edit.fa --prg=gth --cores=24 --gff3 --UTR=on --softmasking --AUGUSTUS_CONFIG_PATH=/users/frodriguez/bin/Augustus/config --GENEMARK_PATH=/bioware/genemark-et-4.32/gmes_petap`


### BAM with UTRs

`braker.pl --species=As11-2 --genome=/users/frodriguez/Avgenome/Av11-2/PBJelly/te/RM/softmasking/Av11_pb_v2.soft.mask.fa --bam=/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L1/accepted_hits.bam,/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L2/accepted_hits.bam --cores=24 --UTR=on --softmasking --AUGUSTUS_CONFIG_PATH=/users/frodriguez/bin/Augustus/config --GENEMARK_PATH=/bioware/genemark-et-4.32/gmes_petap`

#### Output: Wed Jan  2 13:29:12 2019: Logfile: /users/frodriguez/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-2/braker.log!

### BAM + protein

`braker.pl --species=As11-prot --genome=/users/frodriguez/Avgenome/Av11-2/PBJelly/te/RM/softmasking/Av11_pb_v2.soft.mask.fa --bam=/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L1/accepted_hits.bam,/users/frodriguez/illumina/working/As11-RNAseq_180828/tophat/Av11_RNA_L2/accepted_hits.bam --prot_seq=/users/frodriguez/Avgenome/Genoscope/Adineta_vaga_v2.0.annot.pep_edit.fa --prg=gth --cores=24 --softmasking --AUGUSTUS_CONFIG_PATH=/users/frodriguez/bin/Augustus/config --GENEMARK_PATH=/bioware/genemark-et-4.32/gmes_petap`
#### Output:  2 18:22:39 2019: Logfile: /users/frodriguez/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-prot/braker.log!


### This was my final scripts to get the working output:

#######USE ~/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-2/ annotations (with UTRs)
####

###Get gff3 from gff file (convert_augustus_to_gff3.py) & copy UTRs from gff/gtf files. convert_augustus_to_gff3.py DO NOT carry UTRs!!

###sort files 

`~/bin/gff3sort/gff3sort.pl ~/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-2/gff/As11.augustus.hints_braker.nosort.gff3 > ~/Avgenome/Av11-2/PBJelly/BRAKER/As11.augustus.hints_braker.gff3`

`~/bin/gff3sort/gff3sort.pl ~/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-2/augustus.hints_utr.gtf > ~/Avgenome/Av11-2/PBJelly/BRAKER/As11.augustus.hints_braker.gtf`

`~/bin/gff3sort/gff3sort.pl ~/Avgenome/Av11-2/PBJelly/BRAKER/braker/As11-2/augustus.hints_utr.s.gff > ~/Avgenome/Av11-2/PBJelly/BRAKER/As11.augustus.hints_braker.gff`

### Get some statistics

cd ~/bin/eval-2.2.8/

`./get_general_stats.pl ~/Avgenome/Av11-2/PBJelly/BRAKER/As11.augustus.hints_braker.gtf > ~/Avgenome/Av11-2/PBJelly/BRAKER/As11.augustus.hints_braker.gtf.stats`


### Kepp only transcript 1 (.t1), with patterns.txt (IDs to keep) with list of transcripts.t1

```
more patterns.txt
g1.t1
g10.t1
g100.t1
g1000.t1
g10000.t1
g10001.t1
g10002.t1
g10003.t1
```


`awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' As11.augustus.hints_braker.codingseq | grep -Ff patterns.txt - | tr "\t" "\n" > As11.augustus.hints_braker.codingseq.t1.fa`


# Gene prediction: BRAKER3

....
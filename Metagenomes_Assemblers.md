# IDBA / Megahit Assembler
###### tags: `metagenomic`, `assembler`

https://github.com/loneknightpy/idba

Not 100% sure, but I think we need to load all metawrap enviroment to use idba
`module load metawrap`

IDBA series assemblers accept FASTA format reads. FASTQ format reads can be converted by fq2fa program in the package. SO we need to conver files like this:

fasq.gz > fastq > fasta 

First, we have compressed fastq file (.gz). We need to unzip them:

`gunzip SRR*.fastq.gz`

Since we have normally more than one (1 & 2 paired), using the same terminal window, just add "&" at the end of the command:
`gunzip SRR*_1.fastq.gz &`
hit Enter
`gunzip SRR*_2.fastq.gz &`

*What to do if you want to reverse the process, compress a file: use gzip

`gzip SRR*.fastq`

it will create a .fastq.gz file and delete the original .fastq

IMPORTANT: since fastq file are quite big, taking a lot of space, if we are not going to use then any further, ALWAY compress them back into .fastq.gz

OK, going back to the assembler. IDBA series assemblers accept FASTA format reads. FASTQ format reads can be converted by fq2fa program in the package. SO we need to conver files like this:

`fq2fa --merge SRR*_1.fastq SRR*_2.fastq SRR_merge.fa`

AND one more final step, since some libraries (Illumina) are actually made from several runs (fastq files), we need to merge/concatenate all of them into one single file: in our case it will be SCN_BF

cat SRR01*.fa SRR02*.fa ....SRR10*.fa > SCN_BF.fa

## IDBA-UD - Iterative de Bruijn Graph Assembler for sequencing data with highly uneven depth:

*What is a Bruijn Graph Assembler??? One simple tutorial:
https://towardsdatascience.com/genome-assembly-using-de-bruijn-graphs-69570efcc270

inside dir: `/groups/fr2020/metagenomes/Permafrost_RUS/PRJNA680161`

Run the assembly in minnie, with screen session

`idba_ud -r Permafrost680161.fa --num_threads 12 -o idba_ud_out`

## Megahit assembler

This assembler work with the metawrap env:

`module load metawrap`

Basic usage:
https://github.com/voutcn/megahit
https://www.metagenomics.wiki/tools/assembly/megahit
```
megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
```

It works basically with fasta/fastq PE data. With the input used for idba assembler (merge fasta reads):

`megahit --12 3.4.fa --num-cpu-threads 10 --out-dir megahit_assembly`
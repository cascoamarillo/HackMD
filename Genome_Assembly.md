---
title: Genome Assembly
tags: [Assembly]

---

# Lab12 BIOL-530-630

## Intro: Get your computer ready

We will need at least two pieces of software:
the bash terminal and a program to transfer files from and to your computer.


what is a file transfer software? check [here](https://en.wikipedia.org/wiki/FileZilla)

## filezilla

download the appropiate version (Client) [here](https://filezilla-project.org/)


# Connect to oedipus
...now the fun begings:
Open the terminal and type

`ssh -X user@oedipus.bioinformatics.rit.edu `

`password: type here`

#### The 'X' in ssh protocol is to initiate the display. Eg. in MUMMER, when you plot the graph in the terminal a display will pop up on your screen. However, my MacOS doens't like it (LOL).

# Don't remember your password? Come to see me!

## Workspace in oedipus.bioinformatics.rit.edu 

### Class directory
/mnt/classes/biol_530_630

### Lab datasets
/mnt/classes/biol_530_630/data/

### Class Software
/mnt/classes/biol_530_630/bioware/

### YOUR directory
/mnt/classes/biol_530_630/home/xyz1234/

### General Software
/mnt/sde_dir/software/


# For the Lab, you can either create symlinks to your folder with the reads
`ln -s existing_source_file optional_symbolic_link`

copy the entire Lab folder (`cp -R`) 
In your home dir:
`cp -R /mnt/classes/biol_530_630/data/Lab12 /mnt/biol_530_630/home/xyz1234/`

or copy the tar file 
 `cp /mnt/biol_530_630/data/Lab12.tar.gz /mnt/biol_530_630/home/xyz1234/`

Once you finish all your analysis, and want a copy of them in your local computer:

`scp -r user@oedipus.bioinformatics.rit.edu:/mnt/classes/biol_530_630/home/xyz1234/Lab12/ ~/your_local_directory`

`password: type here`
## Just make sure you don't delete/modify anything in
/mnt/biol_530_630/data/Lab12/
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

:::info
Also, make sure you are in the right working directory, and the inputs either are in your work dir or write the full path to them. Same for outputs and scripts in general.
:::

## Running heavy stuff

:::info
## screen

While working in a server, to avoid stoping your run because you loose your connection (wifi...), you can use **screen**. That way, if you lost connection or close the terminal the processstill will be running in the server.

Just log into your user name, and before running anything, invoke `screen`. You'll  see some additional infor about `screen`. Then, hit enter, and a blank terminal will show up, which has your previous environment. Run your script/command prompt, and (if no error shows up) you can close your terminal, turn off your computer, etc. When you log back in `oedipus`, just type `screen -dr`. Your previous terminal screen session/s will show up.

https://linuxize.com/post/how-to-use-linux-screen/
:::

:::info
## algo

Since we are running heavy demanding processes, whith several classes and students, we will implement the use of a special script for rosources allocation:
`algo`
Just type algo in front of your command prompt
Example: `algo python3 hello_world.py`
That way we will work within our memory allocation parameters on the server.
:::
## Genome Assembly comparison

In Lab12 we will run two genome assemblers:

A) Celera assembler (v7): overlap-layout-consensus

B) SPADES (v4): de Bruijn graph

Short reads (100 bp Paired End) are provided from Staphylococcus aureus subsp. aureus USA300 Genome sequencing, which it is a  whole genome shotgun sequencing project of genomic DNA paired-end library 'Solexa-8293' (aka Illumina). More info: 
https://www.ncbi.nlm.nih.gov/sra/SRR022868

Once we have both assemblies, we will compare them with the reference genome (Genome Comparison - Lecture6/Lab6!) with the genome aligner MUMMER.

### SPAdes assembler
https://ablab.github.io/spades/index.html

Spades was originally developed for de novo assembly of genome sequencing data produced for cultivated microbial isolates and for single-cell genomic DNA sequencing. SPAdes is a universal A-Bruijn assembler in the sense that it uses k-mers only for building the initial de Bruijn graph.


Check SPAdes is working on your side:
```
/mnt/classes/biol_530_630/bioware/SPAdes-4.0.0-Linux/bin/spades.py
```
And check all parameter needed and options. Everything looking good? This is a good moment to invoke "screen" if you haven't yet.

Go ahead and run:
```
/mnt/classes/biol_530_630/bioware/SPAdes-4.0.0-Linux/bin/spades.py -1 frag_1.fastq -2 frag_2.fastq -o Saureus_spadesPE -t 4

```

Output (assembly) would be located in 
/Lab12/Saureus_spadesPE/contigs.fasta

### Celera Assembler
https://sourceforge.net/projects/wgs-assembler

`/mnt/classes/biol_530_630/bioware/wgs-7.0/Linux-amd64/bin/runCA -options`

First, we need to process the sort reads.

`/mnt/classes/biol_530_630/bioware/wgs-7.0/Linux-amd64/bin/fastqToCA -libraryname illuminaPE -technology illumina -reads frag_1.fastq -mates frag_2.fastq -insertsize 180 20 > Reads12.frg`

Then, run Celera Assembler

`/mnt/classes/biol_530_630/bioware/wgs-7.0/Linux-amd64/bin/runCA -p Saureus_celera -d runCA_mates12_threads Reads12.frg -s specFile`

Output (assembly) is located in 
/Lab12/runCA_mates12/9-terminator/Saureus_celera.ctg.fasta

### Compare assemblies
MUMmer is a system for rapidly aligning entire genomes.
https://mummer4.github.io/manual/manual.html

MUMmer is installed in `/mnt/sde_dir/software/mummer/mummer-4.0.0rc1/`

We can call it by the full path. There are serveral tools in the MUMmer package. We will use two:

1. `nucmer` is for the all-vs-all comparison of nucleotide sequences contained in multi-FastA data files.
2. `mummerplot `is a perl script that generates gnuplot scripts and data collections for plotting with the gnuplot utility.

**Simple use case (similar assemblies:re-sequencing)**

Given a file containing a single reference sequence (Saureus_ref-genome.fasta) in FASTA format and another file containing multiple sequences in FastA format (Spades_contigs.fasta) type the following at the command line:

`/nucmer  -p <prefix>  ref.seq  qry.seq`

To produce the following files:

`<prefix>.delta`

eg.

`/mnt/sde_dir/software/mummer/mummer-4.0.0rc1/nucmer -p Sa_Celera Saureus_ref-genome.fasta Saureus_celera.ctg.fasta`


To see a simple gnuplot output, run the perl script mummerplot on the output files. Or explore the \<prefix>.[fr]plot file to see the data collection.

You can edit the <prefix>.gp file that is created to change colors, line thicknesses, etc.
    
`./mummerplot  -p <prefix>  <prefix>.delta`

There are some problems (bugs) with the last script. Not sure if it is my end (Mac with XQuartz as display), or more general with the script/configuration.
    
The option `--png` should output a PNG file in your working directory:

`./mummerplot  -p <prefix> --png <prefix>.delta`
    
eg.
    
`/mnt/sde_dir/software/mummer/mummer-4.0.0rc1/mummerplot -p Sa_Celera --png  Sa_Celera.delta`


## Lab12 Activity/Discussion

### Discussion12
Collect all your steps and parameters, and save them into your MarkDown (MD) file (Lab12+Name.md), and export it to pdf (if you can). What is the output of SPAdes//Celera like? How big is the genome assembled? Compared to the reference? How many contigs compose the SPAdes/Celera contigs  versus the reference file?
Tip: there are a thousand ways to do this. One simple way in the terminal:

`grep -v ">" sequence.fasta | grep -E -o "G|C|T|A|N" | wc -l`

`grep ">" -c Saureus_spades.fasta`


### Activity12
Create genome comparison plots, organize them into a single PDF file, and give each plot a descriptive title. Explain what each plot illustrates and the insights it provides. Share your reflections on the results. :)

## Genome Assembly comparison

Do you want to do MORE!?

C) Velvet: de Bruijn graph
Here is another de Bruijn assembler
https://github.com/dzerbino/velvet 

First, velveth helps you contruct the dataset and tell the system what each sequence file represents.

`/mnt/sde_dir/software/velvet/velveth Saureus_velvet 31 -shortPaired -fastq frag_1.fastq frag_2.fastq`

where "31' is the kmer length

Then, velvetg is the core of Velvet where the de Bruijn graph is built then manipulated.

`/mnt/sde_dir/software/velvet/velvetg Saureus_velvet -exp_cov auto -cov_cutoff auto`

Assembly results in 
/Saureus_velvet/contigs.fa
# Analysis I: Get data. SRA & FastQC
---
Edited from Dr. Peredo

Analysis pipeline I: SRA & FastQC 


###### tags: `fastq-dump`, `fastqc`, `trimmomatic`


## Update from July 23, 2020

Get Fastqc report for any .fastq file in a directory:
(do the screen sesion first)
in 
cd /groups/fr2020/metagenomes/bioreactor/
```

ls
post-QC_report  pre-QC_report  SRR3901696_1.fastq  SRR3901696_2.fastq
[frodriguez@flicker SRR3901696]$ 
[frodriguez@flicker SRR3901696]$ fastqc SRR3901696_1.fastq SRR3901696_2.fastq
```

Get .fastq files back to the zipped version:
(do the screen sesion first)
```
ls
post-QC_report  pre-QC_report  SRR3901700_1.fastq  SRR3901700_2.fastq

[frodriguez@minnie SRR3901700]$ gzip SRR3901700_1.fastq &
[1] 49105

[frodriguez@minnie SRR3901700]$ gzip SRR3901700_2.fastq &
[2] 49330
[frodriguez@minnie SRR3901700]$ 
```


## Getting your data ready. Types of files. Quality filtering
Supose er are using Illumina technology. Your data will be provided as fastq files. These files are huge, containing millions of reads, and are not designed for human reading.
Here, we will briefly cover some aspect of how to first interact with these files.

Illumina sequencing steps
https://www.illumina.com/documents/products/techspotlights/techspotlight_sequencing.pdf


### FastQ files
The results of sequencing will be provided to you as fastq files, FASTQ is a text file format (human readable) that provides 4 lines of data per sequence.
> Machine_Sequence_identifier/1-2
> The sequence
> Comments
> Quality scores

In the case of paired-end reads, sequencing results may be stored either in one FASTQ file (alternating) or in two different FASTQ files. Paired-end reads may have sequence identifiers ended by "/1" and "/2" respectively.

Example FASTQ entry for one Illumina read:
```
@EAS20_8_6_1_3_1914/1 ---identifier
CGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACAC
TTTGCTATGCCATAGCATTTTTATCCATAAGATT sequence
+ comments
HHHHHHHHHFHGGHHHHHHHHHHHHHHHHHHHHEHHHHHHHHHHHHHHGHHHGHHHGHIH
HHHHHHHHHHHHHHGCHHHHFHHHHHHHGGGCFHBFBCCF quality
```
[What is a quality score?](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf)

[What are those funny letters in the score?](https://bioinformaticsworkbook.org/introduction/fastqquality-score-encoding.html#gsc.tab=0)

Generally a FASTQ file is stored in files with the suffix .fq or .fastq 

We use compressing programs such as Gzip to reduce the file size. Compression is indicated by the suffix .gz or .gzip.

```
gzip file.3.fastq &
gzip file.2.fastq &

```
decompress 

```
gunzip file.fastq.gz
```

```
or
```

```
gzip -d file.fastq.gz

```




### Fasta files
FASTA it is a text file format in which is sequence is characterized by two lines, 

> one the identifier starting with >, 
> the other the sequence itself. 

No data on read quality is included, but the size of the file is considerable smaller than that of a fastq.
```
>EAS20_8_6_1_3_1914_1 ---identifier
CGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACAC
TTTGCTATGCCATAGCATTTTTATCCATAAGATT sequence
```

## Getting data from NCBI short read archive (SRA)
1. **Go to NCBI https://www.ncbi.nlm.nih.gov/ and select SRA in the database menu. Also, you can search in the project database.**

2. **You can search samples of interest by description, e.g. human microbiome or by usinga  SRA sample number (usually these number will appear in the papers as data are usually required to be deposit in a public database prior publication)**

3. example. Search in the search bar for 'salt marsh amplicon'. Select sample of interest, e.g. 16SRNA-transcript-saltmarshsediment

![](https://i.imgur.com/61LRPrN.png)

![](https://i.imgur.com/mFp1KeV.png)



    a. What is this sample from?
    b. What sequencing technology was used?
    c. When was it published?
    d. Does it belong to a bigger study?

4. **Open the Project accession in a separate tab.**

![](https://i.imgur.com/Dgr2pAO.png)




5. **Click on the SRA experiments link (the number 53).** A new tab will open. You can now send the results to the **run selector** (link in center top). This will provide you access to each sample within the project and also to the metadata of the experiment.
Download the metadata and accession list files. these are text files, separated by commas. You can import the files into excel (search into how to import data, key words for office Data import, text, tabulated commas). to avoid confussion add the project name to the file name

e.g. 

SRR_Acc_List_PRJNA610907.txt
SraRunTable_PRJNA610907.txt



![](https://i.imgur.com/JULPIaU.png)


6. **Click in run SRR11277344. You will navigate to a new menu with lots of info about the sample.**

However, you cannot process or analyze the full scope of the information contained in this file just by looking at it. To really analyze this data set, you will need to download it to your computer.

### Download one SRA file using fastdump
We are going to be downloading the data from the project PRJNA610907

On the terminal navigate to the designed folder in /groups/spart2020


make folder for the raw data

```
cd /groups/fr2020/rawData
mkdir SRAdump_PRJNA610907
```
this will create a folder for you to download your data.
```
cd SRAdump_PRJNA610907
```
navigate to the folder.
we will use fastq-dump to download the data.

check the data characteristics. this is a set of files with [paired reads](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html)
so we want to separate each read in a separate file 1/ and 2/ (we use the split files flag. we will also skip any technical reads that might be present)

[SRR11277344](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11277344)

https://edwards.sdsu.edu/research/fastq-dump/

https://ncbi.github.io/sra-tools/fastq-dump.html 

```
module purge

module load bioware
module load fastqc 

fastq-dump --split-files --gzip --read-filter pass --skip-technical ERR011348
ls

```
**if this does not work**

try
```
module purge

module load bioware
module unload sratoolkit
module load sratoolkit/2.10.7
module load fastqc 
```
The first time you use fastq-dump is going to fail and tell you to configure the program. run this command
```
vdb-config --interactive
```
a blue screen will appear. click save and exit. then enter. the blue screen will go away. now fastq-dump will work. 

now you should be ready to go
```
fastq-dump --split-files --gzip --read-filter pass --skip-technical SRR11277344
ls
```
or Metawrap test libraries

`fastq-dump --split-files --gzip --read-filter pass --skip-technical ERR011347`


`fastq-dump --split-files --gzip --read-filter pass --skip-technical ERR011348`

`fastq-dump --split-files --gzip --read-filter pass --skip-technical ERR011349`


you should see two files.

SRR11277344_pass_1.fastq.gz  SRR11277344_pass_2.fastq.gz



### Checking the QC using fastqc

Now let’s have a look to the quality of the reads. We will be using fastqc, first have a look to
the options. Then, have a look to the data in the read 1 & 2.
```
fastqc -h
fastqc SRR11277344_pass_1.fastq.gz
ls 

firefox SRR11277344_pass_1_fastqc.html 
```

With fastqc you can run more than one file at a time and use option -t (threads) for multiple CPUs:

`fastqc -t 6 SRR11277344_pass_1.fastq.gz SRR11277344_pass_2.fastq.gz`

Open the .html file and have a look to the report.


No sequence has been flagged as bad quality for this file. If that was the case, we can use programs such as trimmomatic to remove low quality reads from our data sets.
in the grey box basic info e.g. number of reads in the file, length.
good quality data stays in the green. (we selected reads that pass the QC, but still it is important to check)

![](https://i.imgur.com/wAMKkIu.png)
 close the firefox window to return to the terminal
 
 lets check read 2. what you see?
 
 ```
 
fastqc SRR11277344_pass_2.fastq.gz
firefox SRR11277344_pass_2_fastqc.html 

```

### trimming reads
You can use the following command to trim your reads. this will remove any potential  (remember to change the names of files to those you are using!)
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```
module load trimmomatic 
java -jar /bioware/trimmomatic-0.36/trimmomatic-0.36.jar PE SRR11277344_pass_1.fastq.gz SRR11277344_pass_2.fastq.gz SRR11277344_pass_paired_1.fastq.gz SRR11277344_pass_unpaired_1.fastq.gz  SRR11277344_pass_paired_2.fastq.gz  SRR11277344_pass_unpaired_2.fastq.gz  ILLUMINACLIP:/bioware/trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

ls

```
It will generate FOUR output files: 1&2 paired reads, and 1&2 unpaired reads
*It is very important the order of output files!
SRR*_1_paired.fastq.gz 
SRR*_1_unpaired.fastq.gz 
SRR*_2_paired.fastq.gz
SRR*_2_unpaired.fastq.gz

in this case, the files are identical to those we were using (we downloaded reads that had already been fitlered for QC.)
to avoid taking to much space in the disk, lets remove the initial set of reads we downloaded

```
rm *_pass_1.fastq.gz
rm *_pass_2.fastq.gz
rm *_unpaired_*

```
Remember that `rm` will remove files FOREVER. To check if you have the right list of files replace `rm` by `echo`. It will show you the list of files without removing them.
```
echo *_pass_1.fastq.gz
echo *_pass_2.fastq.gz
echo *_unpaired_*

```

Of course, if you were analyzing a relatively long set of samples, you would not be doing this one by one!

## Alternative to trimming reads: fastp

In /groups/fr2020/bin/fastp is installed another fastq trimming tool. We are going to use it for Bioproject PRJNA680161, since it contains Nextera adaptor sequences too (along with Illumina). Go to the working dir:

`/groups/fr2020/bin/fastp/fastp -i SRR*_pass_1.fastq.gz -I SRR*_pass_2.fastq.gz -o SRR*_pass_1_clean.fastq.gz -O SRR*_pass_2_clean.fastq.gz -w 2 -j SRR*.json -h SRR*.html`

-w number of CPUs
-j & -h name of report



## Download MULTIPLE SRA files, fast-dump, fastq.
We will create a **for loop**. [more here](https://www.tutorialspoint.com/unix/unix-shell-loops.htm)
All we need is a text file containing the list of files you want to process. Here we will cover an example, using the same dataset as before. In the terminal, get out of the current folder and make a new one.

do you remember the file SRR_Acc_List_PRJNA610907.txt you downloaded?
copy it to the folder. (you can fdrag and drop using filezilla or youcan create a file using the terminal)We will start at 9. 

```
cat > SRR_Acc_List_PRJNA610907.txt
```
this creates the file, let's populate it. 

Enter your document's text. open the txt file in your computer, select all (crt +a) and copy (crt +c) now copy all the text in the terminal **(crt + shift +v)** note: copy and past require crt +shift. 
```
SRR11277344
SRR11277345
SRR11277346
SRR11277347
SRR11277348
SRR11277349
SRR11277350
SRR11277351
SRR11277352
SRR11277353
SRR11277354
SRR11277355
SRR11277356
SRR11277357
SRR11277358
SRR11277359
SRR11277360
SRR11277361
SRR11277362
SRR11277363
SRR11277364
SRR11277365
SRR11277366
SRR11277367
SRR11277368
SRR11277369
SRR11277370
SRR11277371
SRR11277372
SRR11277373
SRR11277374
SRR11277375
SRR11277376
SRR11277377
SRR11277378
SRR11277379
SRR11277380
SRR11277381
SRR11277382
SRR11277383
SRR11277384
SRR11277385
SRR11277386
SRR11277387
SRR11277388
SRR11277389
SRR11277390
SRR11277391
SRR11277392
SRR11277393
SRR11277394
SRR11277395
SRR11277396
```
#Press enter to insert and end line.exit using Ctrl + Z
and check!
```
ls -l SRR_Acc_List_PRJNA610907.txt
head SRR_Acc_List_PRJNA610907.txt

```

### Loop it!

now let's run all the commands we just used in loop.
```

module purge
module load bioware
module unload sratoolkit
module load sratoolkit/2.10.7
module load trimmomatic

for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastq-dump --split-files --read-filter pass --gzip --skip-technical $file; done

for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastqc ${file}_pass_1.fastq.gz; done

for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastqc ${file}_pass_2.fastq.gz; done

for file in `cat SRR_Acc_List_PRJNA610907.txt`; do java -jar /bioware/trimmomatic-0.36/trimmomatic-0.36.jar PE ${file}_pass_1.fastq.gz ${file}_pass_2.fastq.gz ${file}_pass_paired_1.fastq.gz ${file}_pass_unpaired_1.fastq.gz  ${file}_pass_paired_2.fastq.gz  ${file}_pass_unpaired_2.fastq.gz  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36; done
rm *_pass_1.fastq.gz
rm *_pass_2.fastq.gz
rm *_unpaired_*


```
### write a script
you can write all these commands in a small scrip. 

generate a text file (using your note pad or the command line)

here instructions to write it in command line

type in terminal 

```
cat > sra_download_PRJNA610907.sh

```
now copy and past in the terminal

```
 #!/bin/bash
 module purge
module load bioware
module unload sratoolkit
module load sratoolkit/2.10.7
module load trimmomatic
 
 for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastq-dump --split-files --read-filter pass --gzip --skip-technical $file; done
 
 for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastqc ${file}_pass_1.fastq.gz; done
 
 for file in `cat SRR_Acc_List_PRJNA610907.txt`; do fastqc ${file}_pass_2.fastq.gz; done
 
 for file in `cat SRR_Acc_List_PRJNA610907.txt`; do java -jar /bioware/trimmomatic-0.36/trimmomatic-0.36.jar PE ${file}_pass_1.fastq.gz ${file}_pass_2.fastq.gz ${file}_pass_paired_1.fastq.gz ${file}_pass_unpaired_1.fastq.gz  ${file}_pass_paired_2.fastq.gz  ${file}_pass_unpaired_2.fastq.gz  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36; done
 rm *_pass_1.fastq.gz
 rm *_pass_2.fastq.gz
 rm *_unpaired_*
 
```
press enter and crt +Z


 lets activate the script and run it!
 
 ```
 chmod +x sra_download_PRJNA610907.sh
 
 ./sra_download_PRJNA610907.sh
 
 ```
 
 more [here](https://ryanstutorials.net/bash-scripting-tutorial/bash-script.php)
 
 # Wget for EBI ENA fastq files

Instead of fastq-dump, sometimes you have a URL link to your data, just type:
`wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz`

For getting the link, go to 
https://www.ebi.ac.uk/ena/browser/search

search the project/run, and copy the link location:
![](https://i.imgur.com/qq6Fdgq.png)




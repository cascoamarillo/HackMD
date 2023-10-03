# TE denovo annotation

###### tags: `TE`, `annotation`, `REPET`

Let's have some practice with REPET and creating TE diverge plots

## Get your genome (assembly) of interest

We are goin to test it in a small rotifer genome Brachionus koreanus. You can find it in NCBI
https://www.ncbi.nlm.nih.gov/

Select Taxonomy and search "Brachionus koreanus"

![](https://i.imgur.com/50roTCo.png)

Upper left corner: click on Genome link and we'll get some info about this rotifer assembly

![](https://i.imgur.com/3tP1o4Z.png)

Now, time to download the fasta file with the genome; we can download in our system or in the terminal (once you are in the working folder) type:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/177/125/GCA_009177125.1_ASM917712v1/GCA_009177125.1_ASM917712v1_genomic.fna.gz
```
For working folder, in /groups/fr2020/ type

`mkdir Bkoreanus`

Great, we have it saved in our folder!

## REPET environment

Although we have REPET packeage installed and working, we need to set a few things before.

### copy TEdenovo_template
In your working dir: /groups/fr2020/Bkoreanus

`cp -R /groups/arklab/bin/REPET/TEdenovo_template ./`

And now, take a look at the files in it, specially README!!

We also have to look at /groups/fr2020/REPET_README

# This is very important

You'll have to modify your bash_profile, what is that?? it's like an initial configuration for your terminal.
In REPET_README, it says:
```
#------------------------------------------------------------------------------

This directory contains the main REPET program (modified to work in the JBPC's infrastructure) and various files needed to run it.

The "TEdenovo_template" directory provides a template for the TEdenovo pipeline; copy it to use it. The same applies for "TEannot_template" for the TEannot pipeline and "divergencePlot_template" for generating divergence plots.

The bin directory contains several files:
	* The RepBase databases; they are symlinked to by the project templates.
	* The sqldb script; it's used to access the MySQL database.
	* The createRepeatLandscape_UnknownOnTop.pl script; it was copied from /bioware/repeatmasker-4.0.6/util/createRepeatLandscape.pl and modified to place unknowns on top of the graphs generated.

You can use any FASTA file in either pipeline, provided that:
	* There are no more than 15 characters in the name of the file (and the project name)
	* There exists a line break at least once every 60 bases
	* There are no illegal characters in the FASTA headers

#------------------------------------------------------------------------------

###Add these lines to your ~/.bash_profile:

export REPET_PATH="/groups/arklab/bin/REPET/REPET_jbpc"
export tools="/groups/arklab/bin"
export libcpp="$tools/libcppunit-1.12.0/lib"
export LD_LIBRARY_PATH="$libcpp:$LD_LIBRARY_PATH"
export recon="$tools/RECON-1.08/scripts"
export piler="$tools/piler"
export gt="$tools/genometools/bin"
export censor="$tools/censor-4.2.29/bin"
export mreps="$tools/mreps"
export rs="$tools/RepeatScout-1"
export twoBitInfo="$tools/UCSCuserApps/bin"
export modeler="$tools/RepeatModeler"
export rm="$tools/RepeatMasker"
export PATH="$rm:$modeler:$twoBitInfo:$rs:$mreps:$censor:$gt:$piler:$recon:$LD_LIBRARY_PATH:$tools:$REPET_PATH/bin:$PATH"

export GT_NO_FLOCK=yes

#------------------------------------------------------------------------------

Dependencies:

Censor:		girinst.org/downloads/software/censor
CrossMatch:	phrap.org/consed/consed.html#howToGet
	* This utility is delivered alongside Phrap
GenomeTools:	genometools.org/pub
CppUnit:	sourceforge.net/projects/cppunit/files/cppunit/1.12.0
	* Must be version 1.12.0 for some reason; newer or older versions won't work
Mreps:		mreps.univ-mlv.fr/howto.html
Piler:		drive5.com/piler
	* Download the source files
	* In Makefile, remove the --static flag
	* Compile the object files with 'make'
	* Assemble the binary with 'g++ -o piler [list of .o files]'
Recon:		repeatmasker.org/RECON-1.08.tar.gz
RepeatModeler:	repeatmasker.org/RepeatModeler.html
RepeatScout:	bix.ucsd.edu/repeatscout
RMBlast:	repeatmasker.org/RMBlast.html
TwoBit:		hgdownload.soe.ucsc.edu/admin/exe
	* This utility is part of the UCSC userApps package
Repbase:	girinst.org/repbase
	* Symlink the RepBase sequences in the top-level project directories
ESL-Shuffle:	hmmer.org/download.html
	* This utility is part of the HMMER package
TRF:		tandem.bu.edu/trf/trf.download.html


The dependencies listed below should already be installed on the JBPC servers (or at least, I asked Rich to install them):

Text-Soundex:	search.cpan.org/~rjbs/Text-Soundex-3.05/Soundex.pm
RepeatMasker:	repeatmasker.org/RMDownload.html

#------------------------------------------------------------------------------

Changes to various program files:

Edited /groups/arklab/bin/REPET/REPET_jbpc/bin/TEdenovo.py to insert the line:
	* sys.path.append("/groups/arklab/bin/REPET/REPET_jbpc")
directly above this line:
	* from commons.tools.PASTEClassifier import PASTEClassifier

Edited /groups/arklab/bin/RECON-1.08/scripts/recon.pl to insert:
	* "/groups/arklab/bin/RECON-1.08/bin"
on the third line between the double quotes.

Edited /groups/arklab/bin/REPET/REPET_jbpc/bin/GFF3Maker.py to change instances of the text "Identity" to "identity" (necessary for proper formatting of GFF3 files).

#------------------------------------------------------------------------------

Brandon M. Le (Brandon_Le@brown.edu), August 2016
```

### All that is done, but this: Add these lines to your ~/.bash_profile:
Which is a file hidden in your system, you need to tell to show your hidden files.

```
#Add these lines to your ~/.bash_profile:

export REPET_PATH=/groups/arklab/bin/REPET/REPET_jbpc
export tools=/groups/arklab/bin

export libcpp=$tools/libcppunit-1.12.0/lib
export LD_LIBRARY_PATH="$libcpp:$LD_LIBRARY_PATH"
export recon=$tools/RECON-1.08/scripts
export piler=$tools/piler
export gt=$tools/genometools/bin
export censor=$tools/censor-4.2.29/bin
export mreps=$tools/mreps
export rs=$tools/RepeatScout-1
export twoBitInfo=$tools/UCSCuserApps/bin
export modeler=$tools/RepeatModeler
export rm=$tools/RepeatMasker
export PATH="$rm:$modeler:$twoBitInfo:$rs:$mreps:$censor:$gt:$piler:$recon:$LD_LIBRARY_PATH:$tools:$REPET_PATH/bin:$PATH"

export GT_NO_FLOCK=yes

```

We can do it with a text editor and Filezilla help.

Then, it's time to modify the template, fasta file and see if everithing works!

#Read the README file from TEannotation_template:

```
#==============================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPORTANT!!!! TEdenovo.py WORKS with phyton2:

module load python/2.7.15-201901021238

#==============================================================================
#MYSQL
mysql -u arklab_repet -p -h jbpcdb.bpcservers.private
pwd: Gt0qw93#u3tgsgge
mysql> use arklab_repet;
mysql> show tables;
#==============================================================================

This directory serves as a template for running the TEdenovo pipeline.

	1. Copy this directory to another location.

	2. Replace the example "DmelChr4.fa" fasta file with the one you want to run.

	3. Edit TEdenovo.cfg and ./run to include the appropriate project names and working directory.

	4. Execute the following command: "nohup ./run > output.log &"

The output is a MCL-filtered library of TEs located in {NAME}_Blaster_GrpRecPil_Struct_Map_TEclassif_Filtered_MCL. The console output is redirected to output.log.

#==============================================================================

On the FASTA files you use:
	* There is a 15-character limit on the name of the file (and the project name)
	* Make sure there is a line break at least once every 60 bases
	* Make sure there are no illegal characters in the FASTA headers

##LINE BREAK
seqkit seq -w 60 Sp.fasta > Sp.fna

##RENAME fasta headers
awk '/^>/{print ">Sp_ctg" ++i; next}{print}' < Sp.fna > Sp.fst

##ALL CAPITAL LETTERS
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' Sp.fst > Sp.fa

#==============================================================================

Brandon M. Le (Brandon_Le@brown.edu), August 2016
*Edited by Fernando Rodriguez, June 2020
```

##Again, we should work on the dir:
In your working dir: /groups/fr2020/Bkoreanus

`cd /groups/fr2020/Bkoreanus/TEdenovo_template`

Rename TEdenovo_template
.
.
.
help?

### Replace the example "DmelChr4.fa" fasta file with the one you want to run.

Simply, deleted DmelChr4.fa
rm DmelChr4.fa

### Let's get a new fasta name (Spoiler alert: this is very important)
Bkoreanus.fa?
Better not to use funny symbols (*^$@/-_)

### On the FASTA files you should (let's do this every time):
	* There is a 15-character limit on the name of the file (and the project name)
	* Make sure there is a line break at least once every 60 bases
	* Make sure there are no illegal characters in the FASTA headers



##LINE BREAK
```
seqkit seq -w 60 Sp.fasta > Sp.fna
```

##RENAME fasta headers (Two initials plus "_ctg#" seems to work well everytime!).

```
awk '/^>/{print ">Sp_ctg" ++i; next}{print}' < Sp.fna > Sp.fst
```

##ALL CAPITAL LETTERS

```
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' Sp.fst > Sp.fa

```
!!REPET use the fasta file with the .fa extention!!

Make sure file.fa is in TEdenovo directory, you can symlink it (ln -s)

## Edit TEdenovo.cfg and ./run

To include the appropriate project names and working directory

### TEdenovo.cfg

Open TEdenovo.cfg

```
[project]
project_name: DmelChr4
project_dir: /groups/.../TEdenovo
```
Project name is the name of the FASTA file without the extension (.fa)
Project dir, simply check it with `pwd` in the terminal and paste it.

Is your genome polyploid? How many Ns? Check this line out:

`minNbSeqPerGroup: N+1`

### run

The name of the FASTA file you are working with
export PROJECT_NAME=DmelChr4

## Run the script

Better run this in the clusters (cluster5, cricket, grendel)
`status`

	General Use Clusters:

    * cluster5             CentOS 7.4  2 @ 40*995   
    * cricket              CentOS 7.4  1 @ 40*230  
    * grendel              CentOS 7.4  8@8*96     
                                      3@12*128   

`ssh cricket`

`clusterize ./run -n 20`

-n parameter id for the number of nodes (how fast)

For checking your job `qstat`

Alternatively (this probably will take longer), if run it in one of the servers 

```
status
    * arthur               CentOS 7.4    8*96  
    * domino               CentOS 7.4   24*96 
    * flicker              CentOS 7.4   24*250
    * jake                 CentOS 7.4    4*96  
    * luna                 CentOS 7.4    4*96  
    * minnie               CentOS 7.4   32*1000 
    * rocket               CentOS 7.4   12*128 
    * taiga                CentOS 7.4   12*128
    * tern                 CentOS 7.4   12*128
```
`ssh tern`

Open up a screen terminal

`screen`

Check you are in the working dir (TEdenovo)

`./run > output.log`
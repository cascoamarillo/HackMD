# TE divergence Plots
###### tags: `TE`, `divergence`, `plot`
Go to your working dir
`cd /groups/fr2020/...`

Copy the template folder divergencePlot_template in 
/groups/arklab/bin/REPET

`cp -R /groups/arklab/bin/REPET/divergencePlot_template`

- [ ] Take a look at the folder and its contents

- [ ] Read README:


#==============================================================================

This directory contains the scripts and files necessary to output a divergence landscape plot. RepeatMasker's RepBase poorly categorizes animals that aren't mammals or flies (as it appears to me), so using a TEdenovo library is highly suggested:

	1. Copy this folder to a new location.

	2. Symlink the fasta file for the organism's genome to this directory.
`ln -s /groups/fr2020/metagenomes/assembly.fa ./`
`ln -s /TEdenovo_project_directory/{NAME}_Blaster_GrpRecPil_Struct_Map_TEclassif_Filtered_MCL/{NAME}_denovoLibTEs_filtered_MCL.fa ./`


	3. After running TEdenovo on the genome you want to generate a plot for, take the file located in:
	/TEdenovo_project_directory/{NAME}_Blaster_GrpRecPil_Struct_Map_TEclassif_Filtered_MCL/{NAME}_denovoLibTEs_filtered_MCL.fa
	and symlink it to this new location.

	4. Edit "run" to modify the variable names as needed.

	5. Execute the following command: "./run > output.log"

The landscape divergence plot should be labeled {NAME}.html.

#==============================================================================

Notes

The script "extractUnknowns.py" is used to pull out all of the consensuses labeled as "unknown" into a separate file. It also generates a "known" file that contains everything else.

In /groups/arklab/bin/RepeatMasker, there is a script called "createRepeatLandscape_UnknownOnTop.pl" that is a copy of "createRepeatLandscape.pl", with the only difference being that the unknown consensuses are plotted on top of the other results in the divergence plots.

#==============================================================================

Brandon M. Le (Brandon_Le@brown.edu), 2016


OK, step one is done, let's create a symlink (step2) to the fasta file:

ln -s /groups/fr2020/Bkoreanus/Bkoreanus.fa ./

Step 3: Look for the TEdenovo fasta file as described (once TEdenovo is finish) and create a symlink in /groups/fr2020/Bkoreanus/divergencePlot_template

Step: edit `run`:

## Modify these variables as needed
file="DmelChr4"
TElib="DmelChr4_denovoLibTEs_filtered_MCL.fa"

instead of DmelChr, it should be Bkoreanus

We are going to use one of the servers (no clusters) for this run. Then modifie the option "-pa" in
`RepeatMasker -pa 12 -a -low -no_is -e ncbi -lib ${TElib}.classified ${file}.fa`

with the numer of nodes available, but not all of then! It's going to run overnight anyway (so use SCREEN)

Save `run` file

Before run anythig we need to check every commad in run works on your terminal (bash_profile issue). You need to check:

`RepeatClassifier`
`faToTwoBit`
`RepeatMasker`

### Everything ready?

Check which machine are you going to use
`status`
etc...

This step (divergePlot) is better (faster) in the servers:

Generally Available Bioinformatics Servers 

    * arthur               CentOS 7.4    8*96  
    * domino               CentOS 7.4   24*96 
    * flicker              CentOS 7.4   24*250
    * jake                 CentOS 7.4    4*96  
    * luna                 CentOS 7.4    4*96  
    * minnie               CentOS 7.4   32*1000 
    * rocket               CentOS 7.4   12*128 
    * taiga                CentOS 7.4   12*128
    * tern                 CentOS 7.4   12*128


`ssh tern`
Just remener, adjust -pa in ./run

`RepeatMasker -pa 12 -a -low -no_is -e ncbi -lib ${TElib}.classified ${file}.fa`

to the number of processors each machine has (processors*Memory)

Use a screen session
`screen`

and run it
`./run`

If you use one of the clusters
	General Use Clusters:

    * cluster5             CentOS 7.4  2 @ 40*995   
    * cricket              CentOS 7.4  1 @ 40*230  
    * grendel              CentOS 7.4  8@8*96     
                                      3@12*128   
```
ssh cricket
clusterize -log output.log ./run
```



## EXTRA July 9

I would like to check two more divergence plots:
1. Use the combination of cosmid + illumina TEdenovo
`cat Dscosmid_denovoLibTEs_filtered_MCL.fa Dstevenill_denovoLibTEs_filtered_MCL.fa > newname.fa`

For building the plot within the cosmid seqs (Dscosmid.fa) in the divergencePlot_Template 
`customTElib="newname.fa"`

2. And use another custom TE library in 
`/groups/arklab/ostracods/divergencePlot_BL/customLibrary/ostracods_denovoLibTEs_filtered.fa.classified`

copy it in your new divergencePlot_Template folder and:

Run divergencePlot_Template script, but use this script instead: ./run

```
# Genererates a landscape plot through RepeatMasker (RM)

# Modify these variables as needed
file="ostracods"
customTElib="customLibrary/ostracods_denovoLibTEs_filtered.fa.classified"

# Run RM: query against database (Repbase in this scenario) with -a (alignment output .align)
/groups/arklab/bin/RepeatMasker/RepeatMasker -pa 12 -a -low -no_is -lib ${customTElib} ${file}.fa

# Run script calcDivergenceFromAlign.pl using output file (.align) from RM
perl /groups/arklab/bin/RepeatMasker/util/calcDivergenceFromAlign.pl -s ${file}.divsum ${file}.fa.align

# Output file ${file}.divsum which contains Kimura divergence is used to generate the landscape plot
perl /groups/arklab/bin/RepeatMasker/util/createRepeatLandscape_UnknownOnTop.pl -div ${file}.divsum -twoBit ${file}.2bit > ${file}.html
```
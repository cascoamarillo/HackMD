# basic unix commands
###### tags: `basic commands`, `unix`


mv, cp, ls, cd, rm… not familiar?

there are many many tutorials and cheat sheets.

This is a good [cheat sheet](https://files.fosswire.com/2007/08/fwunixref.pdf)

tutorials
http://www.ee.surrey.ac.uk/Teaching/Unix/

https://people.ischool.berkeley.edu/~kevin/unix-tutorial/toc.html

https://www.doc.ic.ac.uk/~wjk/UnixIntro/ #this one has exercises

you should be familiar with:
`mkdir directory` (creates a new directory)
`cd` navigates to a directory
`ls` list the the folders or files in a directory
`ls -l` list the the folders or files with specs
`rm file` removes a file
`rm -r directory` remove directories and their contents recursively
`mv` moves a file
`pwd` actual directory path
`cd ~` home dir
`cd ..` go back one step in the dir path
`cd ../..` go back two steps in the dir path

for all these commands, typing the flag -h or --help usually list the options of usage.

# screen

REALLY important command!
all about the program [here](https://kb.iu.edu/d/acuy)
OR `screen --help` in your bash terminal

Basically, what would happen in you are running something in and the battery laptop dies? You'll get all your new data lost or corrupted (probably). Since we'll run all in the clusters; using the screen you can turn off your computer...the process will keep running remotely. Yeah!

screen # opens a protected window. the processes wont die if you close the terminal

crt +d+a #DeAttach close the terminal screen window, return to regular terminal.

screen -r # return to screen window. if several screen windows are open, there will be a list select the one of interest

IMPORTANT: exit for close the screen session. Once you finish running your script and checkingin your files in your screen session, you have to terminate it by typing `exit`. Or the screen session will be hanging there, until the System Adimistrator will get mad at us!

## nohup
An alternative to screen, where you get a new screen terminal, is using the `nohup` command. Just type
`nohup COMMAND [ARGS]`
like
`nohup ./run > output.log &`
and you can keep working on your terminal while the process is running in the background. But for longer (overnigth) runs, better use screen.

# use the tab in your keyboard!

It will fill up info for you so you dont have to type us much. if two files share a name, it will stop. type a little more and use the tab again. Or press it 3 times and will show all posibilities.

:)
# how to make a text file in terminal

type in terminal

`cat > fileName.txt`

now write (or paste) in the terminal

 `text of your document!`
 

press enter and crt +Z

use head to check contain.

`head fileName.txt`

This can be done by other ways too, yo can find yours!

# how to write a simple script

you can write all these commands in a small scrip.

generate a text file (using your note pad or the command line)

here instructions to write it in command line

type in terminal

`cat > sra_download_PRJNA610907.sh`

now copy and past in the terminal
This is just an example using SRA data, which doesn't have to be our case.
 ```
#!/bin/bash
 module purge
 module load bioware
 module load fastqc 
 module load trimmomatic
```
```
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

`chmod u+x sra_download_PRJNA610907.sh`

`./sra_download_PRJNA610907.sh`

Do you want to know [more](https://ryanstutorials.net/bash-scripting-tutorial/bash-script.php)?

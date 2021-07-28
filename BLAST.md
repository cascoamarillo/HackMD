# BLAST

## Tutorials

Blast hand-book:
https://www.ncbi.nlm.nih.gov/books/NBK279690/

Nice & more friendly tutorial:
https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/


In sum... blastn, blastp, tblasnx??

![](https://i.imgur.com/AwYH5r1.png)

TRY running blastn -help to see them for the blastn program.
Here is a summary of a few parameters that are most commonly used for blast:

    -query <fasta file>
        The name (or path) of the FASTA-formatted file to search for as query sequences.
    -db
        database name
    -evalue <real number>
        Only HSPs with E values smaller than this should be reported. For example: -evalue 0.001 or -evalue 1e-6.
    -outfmt <integer>
        How to format the output. The default, 0, provides a human-readable (but not programmatically parseable) text file. The values **6 and 7 produce tab-separated rows and columns in a text file, with 7 providing explanatory comment lines**. Similarly, a value of 10 produces comma-separated output; 11 produces a format that can later be quickly turned into any other with another program called blast_formatter. Options 6, 7, and 10 can be highly configured in terms of what columns are shown.
    -max_target_seqs <integer>
        When the output format is 6, 7, or 10 for each query sequence, only report HSPs for the first <integer> different subject sequences.
    -max_hsps <integer>
        For each query/target pair, only report the best <integer> HSPs.
    -out <output file>
        Write the output to <output file> as opposed to the default of standard output.         
    -num_threads <integer>
        Number of CPUs
        
### Output types we use:
```
-outfmt 7
    table style with line comments
-outfmt 0
    Blast hit descriptions (names) and aligment
```

I normally like to have both; run first with "7" and then with "0", and use different names as output.

### Databases
We used different database based on the search type:
blastn
    `-db nt`
blastp, blastx
    `-db nr`

## Examples

### blastn
`blastn -query Perma827_Rotifera_GC50.fasta -db nt -evalue 0.001 -outfmt 7 -max_target_seqs 10 -num_threads 6 -out Perma827_Rotifera_GC50.blasn-nt7.txt`


### blastx
`blastx -query Perma827_Rotifera_GC50.fasta -db nr -evalue 0.001 -outfmt 0 -max_target_seqs 10 -num_threads 6 -out Perma827_Rotifera_GC50.blasx-nr0.txt`


## XML file output

Another useful output format is .xml (-outfmt 5), which produce the resourceful .xml that can be loaded into other applications, like Blast-Viewer or EPOS, java based GUI programs:


https://github.com/pgdurand/BlastViewer

https://bio.informatik.uni-jena.de/epos/downloads/


### Example for BlastViewer
```
blastn -query Perma827_Rotifera_GC50.fasta -db nt -evalue 0.001 -outfmt 5 -max_target_seqs 10 -num_threads 6 -out Perma827_Rotifera_GC50.blasn-nt.xml
```

Then 

`java -jar blastviewer-5.4.1.jar -in "./data/Perma827_Rotifera_GC50.blasn-nt.xml"`

## Limiting a Search by Taxonomy

To limit a search by taxonomy, you would need access to the taxid or taxid list of a specifc organism or species. 

To retrieve that information, you would need to install EDirect Utilities on your machine and update your bashprofile accordingly.

### Downloading EDirect
You can download EDirect using this link and following the steps under "Installation:" https://www.ncbi.nlm.nih.gov/books/NBK179288/ "

Here is a summary of the process:
1. Download the [EDirect installer](https://www.ncbi.nlm.nih.gov/books/NBK179288/bin/install-edirect.sh) to retrieve a .sh script
2. Using FileZilla or any other software, create a copy of the .sh script you just downloaded in your home directory on MBL servers. (You can reach your home directory by running the `cd` command on the terminal.) 
3. In the directory you just added the file in, run this command: `source ./install-edirect.sh` 
4. After it is done installing, you will be asked if you want to automatically add a path to this directory in your bash profile. You can type `y` to add automatically, or you can manually run the following command
`echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile` 

### Using taxids in the blast command
To limit a search by taxonomy, you can follow the steps outlined in this [website](https://www.ncbi.nlm.nih.gov/books/NBK569846/).
#### Limiting a BLAST search with a high-level taxonomic node
1. Get the taxid of a pecies from the NCBI website. (Rotifera: 10190)
2. Run `get_species_taxids.sh -t 10190 > 10190.txids` which will produce a file of all Rotifer taxids
3. Edit your blast call by adding the -taxidlist. 
```
blastn -query contigs.fa -taxidlist 10190.txids -db nt -evalue 0.001 -outfmt 5 -max_target_seqs 5 -num_threads 6 -out contigs.blasn-nt5Rotifera.xml
```
#### Limiting a BLAST search with a species-level taxonomic node

For this, you can just use the taxid of the species. 
```
$ blastn –db nt –query QUERY –taxids 9606 –outfmt 7 –out OUTPUT.tab
```






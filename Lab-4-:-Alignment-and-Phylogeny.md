Today you will align gene sequences from 10 species and build a gene tree. The evolutionary histories of individual genes can be quite different than the evolutionary history of species, though, so your resulting gene tree may or may not match the true species tree. If this type of gene tree/species tree discordance interests you, be sure to check out Coalescent Theory:
![](https://biologos.org/files/resources/dnavariant.png)

There will be one major difference in lab today, though -- you will not be given the exact commands to type on the command line. I will point you to each program and the manual page for you to run. When in doubt, use defaults to begin with, and then google.

Log onto an interactive node before you start.

## Download cds and peptide sequences
I have placed two fasta files in /home/student/binf4550/data/04.Phylogeny

Make a new directory (with mkdir) and copy them into your current directory. 

    interesting_seqs.cds
    interesting_seqs.pep

These two fasta files are paired, meaning that the cds and peptides correspond to each other. The length of each cds sequence should be 3 times the length of every peptide sequence. How would you use grep to quickly look at the fasta headers in each of these files?

It makes sense to double-check your files sometimes. I often use the head command to look at the first few lines of multiple files. In the cds file, the Cow sequence starts with an ATG (start codon) and the in the peptide file, Cow starts with an M (Methionine start codon). Looking good there!

![](http://i.imgur.com/BhtJW8n.png)

OK, your first task is to take the cow peptide sequence and use BLAST to identify its closest homolog and annotation. Which flavor of BLAST would you most reasonably use for this scenario?
    blastn, blastp, blastx, tblastn, tblastx ?

## Perform multiple sequence alignments with a variety of programs

1) Aligning sequences with MUSCLE:

Muscle is a good, all-around multiple sequence aligner.

    /usr/local/muscle/latest/bin/muscle -h
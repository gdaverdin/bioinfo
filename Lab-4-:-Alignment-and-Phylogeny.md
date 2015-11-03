Today you will align sequences from 10 species and build a phylogeny. There will be one major difference in lab today, though -- you will not be given the exact commands to type on the command line. Instead, that is your homework. I will point you to each program and the manual page for you to run. When in doubt, use defaults to begin with, and then google.

## Download cds and peptide sequences
I have placed two fasta files in /home/student/binf4550/data/04.Phylogeny

Make a new directory (with mkdir) and copy them into your current directory. 

    interesting_seqs.cds
    interesting_seqs.pep

These two fasta files are paired, meaning that the cds and peptides correspond to each other. The length of each cds sequence should be 3 times the length of every peptide sequence. How would you use grep to quickly look at the fasta headers in each of these files?

It makes sense to double-check your files sometimes. I often use the head command to look at the first few lines of multiple files.

![](http://imgur.com/lOBgIQ6)

## Perform multiple sequence alignments with a variety of programs

1) Aligning sequences with MUSCLE:

    /usr/local/muscle/latest/bin/muscle -h
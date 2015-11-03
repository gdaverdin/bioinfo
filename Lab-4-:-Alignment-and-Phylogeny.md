Today you will align gene sequences from 10 species and build a gene tree. 

The evolutionary histories of individual genes can be quite different than the evolutionary history of species, though, so your resulting gene tree may or may not match the true species tree. If this type of gene tree/species tree discordance interests you, be sure to check out Coalescent Theory:
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

OK, your first task is to take the cow peptide sequence and use BLAST to identify its closest homolog and annotation. Use NCBI, Uniprot, whatever you prefer. Which flavor of BLAST would you most reasonably use for this scenario?

    blastn, blastp, blastx, tblastn, tblastx ?

Once you know the gene ID, you can continue.

## Perform multiple sequence alignments with a variety of programs

Multiple sequence alignment is computationally difficult. This has been a developing mathematical and empirical field over the last 20+ years. As with much informatics, there is no single program that will work perfectly for all applications. Some aligners are better with peptides, some better with short sequences, some better with long, some are more "gappy" than others (aka they tend to open and extend gaps more freely).  

1) Aligning sequences with MUSCLE:

Muscle is a good, all-around multiple sequence aligner. Check out the manual page and learn how to use Muscle with default settings. Generate a multiple sequence alignment for the peptide file.

    /usr/local/muscle/latest/bin/muscle -h

2) Aligning sequences with Prank

Prank uses a different alignment algorithm with different penalties for opening and extending gaps.

    /usr/local/prank/latest/prank -h

3) Aligning sequences with PASTA

PASTA is my current favorite of the bunch. It builds on concepts developed in SATE to break large alignments down into small problems and generates Hidden Markov Models for alignment backbones. Specifically, PASTA utilizes transitivity to quickly improve alignments (ie., if A -> B and B -> C, then A must align to C). 

https://github.com/smirarab/pasta

    python2.7 /usr/local/pasta/latest/run_pasta.py -h

What obvious differences did you notice in the time it takes to run all three programs? Remember that this is a small alignment. Time matters, sometimes. More importantly, did the three programs produce different alignments? What can you say about the propensity for each aligner to open a gap? When I want to quickly look at small alignments with a nice color scheme, I keep this page bookmarked: http://www.ebi.ac.uk/Tools/msa/mview/


## Build a peptide gene tree in RAxML

Let's use our PASTA alignment (pastajob.marker001.interesting_seqs.pep.aln) to build a gene tree. RAxML is my recommended choice for building maximum likelihood (ML) trees. Here we will use the GTR (General Time Reversible) model of protein substitution with the gamma model of rate heterogeneity (because not all species evolve at the same rate!), running 100 rapid bootstraps. We will use 2 threads.

run_raxml.sh
    
    #!/bin/bash
    /usr/local/raxml/latest/raxmlHPC-PTHREADS -T 2 -f a -x 12345 -m PROTGAMMAGTR -s pastajob.marker001.interesting_seqs.pep.aln -# 100 -p 2 -n COX2alignment

and submit it, making sure to request two threads.
    
    qsub -q rcc-30d -pe thread 2 ./run_raxml.sh


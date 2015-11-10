So far in this course, we have been working on one analysis at a time. What if we resequence two, ten, or hundreds of individuals, and want to run the same analysis on all individuals? Since we don't want to individually type commands, then wait, then type again -- we can make the computer work for us while we do something more fun. To do this, we can use loops. Today we will revisit the human brain and adrenal RNAseq from the transcriptomics lab and work on an automation pipeline to save us time. 

The most basic concept of a loop is that we "loop over" multiple things and do something.

## The for loop

Here is a basic skeleton example of a for loop in bash.

    #!/bin/bash
    for i in *.fastq
    do
    SOMETHING INTERESTING
    done
    

Let's start by running FASTQC on the brain and adrenal RNAseq fastq files. These are single end reads (not paired end). LOG ONTO AN INTERACTIVE NODE USING QLOGIN. Copy the fastqs, chr19.fa, and the gtf gene models into your own directory. I have placed the data here:

    aharkess@zcluster:/home/student/binf4550/data/06.Loops/human_rnaseq$ ls -lh
    total 88M
    -rw-rw---- 1 aharkess jlmlab 7.9M Oct 19 13:56 Adrenal.fq
    -rw-rw---- 1 aharkess jlmlab 6.0M Oct 19 13:56 Brain.fq
    -rw-rw---- 1 aharkess jlmlab  58M Oct 19 13:56 chr19.fa
    -rw-rw---- 1 aharkess jlmlab 5.9M Oct 19 13:57 USCS_hg19_chr19.genes.gtf

To run FASTQC on both of our fastq files in one fell swoop, we can write a loop. Let me show you an example first, then we can break it down afterwards. I called mine "batch_fastqc.sh".

    #!/bin/bash
    for i in *.fastq
    do
    echo "Beep boop. I'm a computer running fastqc on file $i"
    /usr/local/fastqc/latest/fastqc $i
    done

What this loop does in "pseudocode": For every file "i" that ends in *.fastq in my current directory, print the phrase "Beep boop. I'm a computer running fastqc on file $i" (see how I referenced the current fastq file by using $i ?). Then run fastqc using $i as the input. Then move on to the next file that ends with *.fastq in my directory, and do it all over again.

To run this code on the interactive node, we can type:

    sh ./batch_fastqc.sh

Or we can submit it to the queue using 

    qsub -q rcc-30d ./batch_fastqc.sh


**Weird. This code isn't producing any output. Find and fix the error.**

## Write your own loop

Spend a few minutes to write your own loop to count the number of reads in each fastq file. I would like the output to look like this:

    aharkess@compute-14-7:/home/student/binf4550/data/06.Loops/human_rnaseq$ sh ../count_reads.sh
    The read count of file Adrenal.fq is:
    50121
    The read count of file Brain.fq is:
    37992

## Using command substitution

We can also insert the output from commands into a bash script using the backticks

    `command`

Here's an example:

    #!/bin/bash
    for i in *.fq
    do
    echo "The read count of file $i is:`grep -c "^@ERR" $i`"
    done

## Practical exercise: Count the number of sequences in a directory filled with fasta files

Now you're going to write your own loop again. This time, you're going to write a loop to go over every fasta file in a directory, then print out a two column file with 

Fasta_name     Number_of_sequences

Start out by copying a directory of multiple sequence alignments for the pineapple genome from /home/student/binf4550/data/06.Loops/pineapple_alignments/



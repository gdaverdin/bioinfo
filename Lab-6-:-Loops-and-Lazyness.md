So far in this course, we have been working on one analysis at a time. What if we resequence two, ten, or hundreds of individuals, and want to run the same analysis on all individuals? Since we don't want to individually type commands, then wait, then type again -- we can make the computer work for us while we do something more fun. To do this, we can use loops.  

The most basic concept of a loop is that we "loop over" multiple things and do something. We will talk about two types of loops today, the "for" loop and the "while" loop. We will implement both types in the bash language in unix. 

## The for loop

Here is a basic skeleton example of a for loop in bash.

    #!/bin/bash
    for i in *.fastq
    do
    SOMETHING INTERESTING
    done
    

Let's start by running FASTQC on a set of FASTQ files. 

## The while loop
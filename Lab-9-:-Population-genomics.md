Diversity is the spice of life. 

Today you will be using the Human 1000 genome project data from a set of individuals to characterize diversity along chromosome 1. We will calculate allele frequency, measures of pairwise diversity (pi), signatures of selective sweeps (Tajima's D), transition/transversion (Ts/Tv) rate, and other various population metrics. 

Luckily, once you have a VCF file of high quality SNPs, running these analyses is easy. Like most informatics, the difficult part is knowing the underlying biology and mathematics, then interpreting the results. We will be using VCFtools intensively today. Spend some time now to read the manual carefully and check out all the available options and analyses it can run. 

# Copy data

I have placed a gzip-compressed vcf file of 2,504 human genome SNP annotations here:
    
    /home/student/binf4550/data/09.PopulationGenomics/1.1-4000000.ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

This file is a "slice", only including positions 1-4000000 on chromosome 1. 

# Calculate pairwise nucleotide diversity (pi)

Let's start by using VCFtools to calculate pi, or pairwise nucleotide diversity. You can think of pi as "If you grab two random alleles out of a hat, what's the probability that they are different?"

    /usr/local/vcftools/latest/bin/vcftools --gzvcf 1.1-4000000.ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --site-pi --chr 1 --out chr1_pi.out 

Great. Look at the head of your *sites.pi file. It has 3 columns (CHROM, POS, PI). Download that data to your desktop and let's plot it using R. 

Let's just plot this in base R real quick.

    setwd('/wherever/you/put/your/data')
    dat <- read.table('pi_test.sites.pi', header=TRUE) # or whatever you named your file.
    plot(dat$PI~dat$POS)

![](http://i.imgur.com/2GI24Gb.jpg)

X axis is position, Y axis is Pi. Well that's a pretty ugly and hard to interpret plot. Can you make a better one?

Google how to calculate the mean Pi for all individuals in R. The answer is 
    
    [1] 0.03279469

Hint: You can reference a specific column from a dataframe using the dollar sign, like this

    dat$PI

# Calculate Fst 

Fst is a measure of population differentiation. High Fst indicates that two populations are fixed for allelic differences. 

![](http://www.nature.com/nrg/journal/v5/n8/images/nrg1401-i1.jpg)

Read the VCFtools manual to calculate Fst in 1kb sliding windows (--fst-window-size). Plot the data in R.

# Challenge: Look for evidence of a selective sweep using Tajima's D

From wikipedia:

> A negative Tajima's D signifies an excess of low frequency polymorphisms relative to expectation, indicating population size expansion (e.g., after a bottleneck or a selective sweep) and/or purifying selection. A positive Tajima's D signifies low levels of both low and high frequency polymorphisms, indicating a decrease in population size and/or balancing selection. 

My output looks like this:

    CHROM	BIN_START	N_SNPS	TajimaD
    1	10000	11	-0.301026
    1	11000	3	-0.0088863
    1	12000	0	0
    1	13000	21	-1.55401
    1	14000	14	-0.381699
    1	15000	25	-0.941925
    1	16000	17	-1.97147
    1	17000	12	-1.84959
    1	18000	2	-0.512236

Plot the results and look for  
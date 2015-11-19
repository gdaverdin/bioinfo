Diversity is the spice of life. 

Today you will be using the Human 1000 genome project data from a set of individuals to characterize diversity along chromosome 1. We will calculate allele frequency, measures of pairwise diversity (pi), signatures of selective sweeps (Tajima's D), transition/transversion (Ts/Tv) rate, and other various population metrics. 

Luckily, once you have a VCF file of high quality SNPs, running these analyses is easy. Like most informatics, the difficult part is knowing the underlying biology and mathematics, then interpreting the results. 

# Copy data

I have placed a gzip-compressed vcf file of 2,504 human genome SNP annotations here:
    
    /home/student/binf4550/data/09.PopulationGenomics/1.1-4000000.ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

This file is a "slice", only including positions 1-4000000 on chromosome 1. 

# Calculate pairwise nucleotide diversity (pi)

Let's start by using VCFtools to calculate pi, or pairwise nucleotide diversity. You can think of pi as "If you grab two random alleles out of a hat, what's the probability that they are different?"

    /usr/local/vcftools/latest/bin/vcftools --gzvcf 1.1-4000000.ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --site-pi --chr 1 --out chr1_pi.out

# Calculate Fst 

# Calculate Ts/Tv


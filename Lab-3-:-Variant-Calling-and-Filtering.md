You have just been hired as a bioinformatician in at the Fred Hutchinson Cancer Research Institute. Your first task is to analyze BRCA1 amplicon resequence for a 45 year old female individual that has developed breast cancer. Breast cancer is currently the most common type of cancer in female humans, so clinical diagnostics must be fast and accurate. Your job is to identify any mutations and characterize their potential outcomes.

You will be using the GATK pipeline and best practices to align your reads and call SNPs. Since this dataset is so small, for the sake of time, we can run this entire analysis on the interactive node without creating submission scripts. Just enter the commands, do not create a submission script and use qsub. 

## Download amplicon reads
I have simulated paired-end reads from a breast cancer-positive individual and placed them in our lab directory. They called Brca1Reads_0.1.fastq and Brca1Reads_0.2.fastq. Get on an interactive node and copy them into your directory on the cluster. 

    /home/student/binf4550/data/03.VariantCalling

## Aligning shotgun reads to the human genome with bwa-mem
For this section, I am following this page of the best practices guide : https://www.broadinstitute.org/gatk/guide/article?id=2799

GATK is very particular. It requires that each sequencing library have metadata attached to it. These metadata will end up in our .bam alignments, which is nifty for a lot of downstream reasons. First we have to compose a short string of this metadata, called the read group identifier, in the following format:

    @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 
    where the \t stands for the tab character.

    @RG = an indicator that this line is a read group identifier line
    ID = a group ID
    SM = sample ID
    PL = platform (Illumina, 454, Ion torrent, Pacbio)
    LB = library number (can have multiple libraries for a given individual)
    PU = platform unit (for big facilities that have multiple sequencers)

We can just use the default @RG string above. 

Then we have to index the reference genome (in this case, were just looking at chromosome 17). 

    /usr/local/bwa/latest/bwa index chr17.fa

Now we can run bwa-mem. BWA is the Burrow-Wheelers Aligner written by alignment guru Heng Li. bwa-mem is a workhorse in the alignment world; it balances speed with alignment sensitivity and accuracy.

Here's what the manual tells us to run:

    bwa mem -M -R ’<read group info>’ -p reference.fa raw_reads.fq > aligned_reads.sam

But we have to use the path to bwa on our cluster. Remember, on the zcluster, programs are kept in /usr/local/:

    /usr/local/bwa/latest/bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' chr17.fa Brca1Reads_0.1.fastq Brca1Reads_0.2.fastq > Brca1Reads_aligned.raw.sam

## Convert SAM -> BAM, sort, and mark PCR duplicates

SAM is big. Converting to BAM right away saves us some some valuable hard drive space. Then we need to sort every read by its alignment coordinate relative to the reference.



## Calling SNPs and indels with the GATK Unified Haplotyper



## Using SnpEff to annotate variants
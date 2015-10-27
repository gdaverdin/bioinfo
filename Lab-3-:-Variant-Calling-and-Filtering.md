You have just been hired as a bioinformatician in at the Fred Hutchinson Cancer Research Institute. I have simulated DPCK gene amplicon resequence for a 45 year old female individual that has developed pancreatic and breast cancer. Your job is to identify any mutations in this individuals' genes that may be causal.

You will be using the GATK pipeline and best practices to align your reads and call SNPs. Since this dataset is so small, for the sake of time, we can run this entire analysis on the interactive node without creating submission scripts. Just enter the commands, do not create a submission script and use qsub. 

## Download amplicon reads
I have simulated paired-end reads from a breast cancer-positive individual and placed them in our lab directory. They called Brca1Reads_0.1.fastq and Brca1Reads_0.2.fastq. I have also placed a fasta file for the human genome chromosome 17 (chr17.fa) and a .gtf of the UCSC gene models for this chromosome (Hg19.Chr17.UCSC-3.gtf) Get on an interactive node and copy this entire directory into your directory on the cluster. The directory is located here:

    /home/student/binf4550/data/03.VariantCalling

If you don't know how to copy the contents of an entire folder into your directory, try to google it before you ask for help. 

## Aligning shotgun reads to the human genome with bwa-mem
For this section, I am following this page of the best practices guide : https://www.broadinstitute.org/gatk/guide/article?id=2799

Check out the GACRC page for BWA, as well. https://wiki.gacrc.uga.edu/wiki/Burrows-Wheeler_Aligner_%28BWA%29

GATK is very particular. It requires that each sequencing library have metadata attached to it. These metadata will end up in our .bam alignments, which is nifty for a lot of downstream reasons. First we have to compose a short string of this metadata, called the read group identifier, in the following format:

    @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 
    where the \t stands for the tab character.

    @RG = an indicator that this line is a read group identifier line
    ID = a group ID
    SM = sample ID
    PL = platform (Illumina, 454, Ion torrent, Pacbio)
    LB = library number (can have multiple libraries for a given individual)
    PU = platform unit (for big facilities that have multiple sequencers)

We can just use the default @RG string above. This Readgroup (@RG) string becomes particularly important when we have multiple samples to align and analyze. 

Then we have to index the reference genome (in this case, were just looking at chromosome 17). 

    /usr/local/bwa/latest/bwa index chr17.fa

Now we can run bwa-mem. BWA is the Burrows-Wheeler Aligner written by alignment guru Heng Li. bwa-mem is a workhorse in the alignment world; it balances speed with alignment sensitivity and accuracy.

Here's what the manual tells us to run:

    bwa mem -M -R ’<read group info>’ -p reference.fa raw_reads.fq > aligned_reads.sam

But we have to use the path to bwa on our cluster. Remember, on the zcluster, programs are kept in /usr/local/:

    /usr/local/bwa/latest/bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' chr17.fa Brca1Reads_0.1.fastq Brca1Reads_0.2.fastq > Brca1Reads_aligned.raw.sam

## Convert SAM -> BAM, sort, and mark PCR duplicates

SAM is big. Converting to BAM right away saves us some some valuable hard drive space. Then we need to sort every read by its alignment coordinate relative to the reference so that it is easily searchable. Afterwards, we want to filter out reads derived from PCR duplicates. These are not biologically informative reads, so we don't want them to contribute as depth evidence for a SNP.

Picard is an incredible set of tools to manipulate sequence and analysis formats. Let's see what all it can do

    java -jar /usr/local/picard/latest/dist/picard.jar -h

Now let's **sort** our .sam file using picard and convert it to bam.

    java -jar /usr/local/picard/latest/dist/picard.jar SortSam INPUT=Brca1Reads_aligned.raw.sam OUTPUT=Brca1Reads_aligned.sorted.bam SORT_ORDER=coordinate

And then **mark duplicate reads** so they don't get counted during SNP calling.

    java -jar /usr/local/picard/latest/dist/picard.jar MarkDuplicates INPUT=Brca1Reads_aligned.sorted.bam OUTPUT=Brca1Reads_aligned.sorted.dedup.bam METRICS_FILE=metrics.txt

Now we have to **index** our sorted, duplicate-marked .bam alignment file.

    java -jar /usr/local/picard/latest/dist/picard.jar BuildBamIndex INPUT=Brca1Reads_aligned.sorted.dedup.bam

## Perform local re-alignments around putative indel sites

Almost ready. Lastly we need to realign indel sites. Remember that we need to do small local alignments here to correct for alignment issues around gap opening and extension sites. This will ensure that we get true indel calls instead of spurious false-positive SNPs . Index the genome real quick, once to create a dictionary for GATK, and the other to generate a regular index using samtools faidx:

    java -jar /usr/local/picard/latest/dist/picard.jar CreateSequenceDictionary REFERENCE=chr17.fa OUTPUT=chr17.dict
    samtools faidx chr17.fa

Then create a list of putative indel sites that need to be realigned.

    java -jar /usr/local/gatk/latest/GenomeAnalysisTK.jar -T RealignerTargetCreator -R chr17.fa -I Brca1Reads_aligned.sorted.dedup.bam -o targets_for_realignment.list

Then do the realignment.

    java -jar /usr/local/gatk/latest/GenomeAnalysisTK.jar -T IndelRealigner -R chr17.fa -I Brca1Reads_aligned.sorted.dedup.bam -targetIntervals targets_for_realignment.list -o Brca1Reads_aligned.sorted.dedup.realigned.bam


## Calling SNPs and indels with the GATK Unified Haplotyper

    java -jar /usr/local/gatk/latest/GenomeAnalysisTK.jar -T HaplotypeCaller -R chr17.fa -I Brca1Reads_aligned.sorted.dedup.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o Brca_raw_variants.vcf

Remember that this is an unfiltered set of SNPs. We can apply some hard filters to our variants for depth, call quality. This is a difficult thing to teach -- you'll need to read about this yourself and cater variant filtration to your own experiment (https://www.broadinstitute.org/gatk/guide/article?id=2806). 

## Download and visualize

You have used IGV to visualize your .bam files before. Download the .bam, .bam.bai, chr17.fa, the .gtf file, and your .vcf file. You can drag and drop the .bam and .vcf files into the IGV viewer to overlay them. Here is an example of what mine looks like:

![](http://i.imgur.com/ha95pnD.png)

Boy, downloading and visualizing SNPs to figure out if they're in a coding region is kind of a pain. I wonder if we could do this computationally instead of by eye...

## Homework:

1) Are there any coding SNPs or indels? At what position? Use 1-based coordinates. 


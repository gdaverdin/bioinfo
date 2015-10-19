Today's lab will be a "choose your own adventure" based on your personal research needs. You can pick between 

* building a _de novo_ transcriptome assembly, annotation, and quantification, or
* aligning RNAseq reads from adrenal and brain tissue to the human genome, and finding differential expression between two conditions


## _De novo_ transcriptome assembly

You have formed a collaboration with an orchid biologist who studies the "bee orchid", _Ophrys apifera_, which is an incredible example of mimicry. The bee orchid mimics the coloration and scent of a bee species that pollinates it. Male bees will attempt to copulate with the orchid, confusingly thinking that it is a female bee. The male bees unsuspectingly collect pollen on their bodies, move on to another flower, attempt to copulate it, and ultimately pollinate it. Your collaborator wants to identify candidate genes that are potentially responsible for this coloration and scent variation, so he generated a small amount of RNAseq data to assemble and annotate expressed transcripts.

![](https://thmcf.files.wordpress.com/2013/06/bee-orchid-imc-3702.jpg) 

For the sake of time, you do not have to assemble a transcriptome today if you would rather focus on the downstream steps. Here is a walkthrough of how I assembled this one, though, if you would like to follow along and learn.

1) Download paired-end 75nt RNAseq reads from SRA

    /usr/local/sra/latest/bin/fastq-dump --split-files SRR609403

2) Run FASTQC and check the quality of the reads

    /usr/local/fastqc/latest/fastqc SRR609403_1.fastq SRR609403_2.fastq

3) Run Trinity to _de novo_ assemble the transcriptome

Conveniently, Trinity has a lot of built-in features that will save you time. For instance, it can run Trimmomatic to clean your RNAseq data of adapter contamination and poor quality bases at the 3' ends of reads.

    #!/bin/bash
    export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
    export PATH=/usr/local/gmap-gsnap/latest/bin/:${PATH}
    time /usr/local/trinity/2.0.6/Trinity --seqType fq --max_memory 10G --left SRR609403_1.fastq --right SRR609403_2.fastq --CPU 2 --trimmomatic --full_cleanup

and submit this job to the rcc-30d queue, requesting 2 threads (because we asked Trinity to use 2 with the --CPU flag), and telling the cluster that were going to use at least 10G of RAM on a node. This will bump us onto a node with 48Gb of RAM. 

    qsub -q rcc-30d -cwd -pe thread 2 -l mem_total=20G ./run_trinity.sh


## Reference-based alignment and differential expression

You have isolated RNA and generating sequencing libraries from brain and adrenal gland tissues to understand what genes are differentially expressed between these two tissue types. For speed, we will be working with just a small subset of data that maps to a small region of chromosome 19. I have placed a link to the RNAseq data, the human genome chromosome 19, and the human genome gene models here:

http://jlmwiki.plantbio.uga.edu/~aharkess/human_rnaseq.tar.bz2

Use wget to download this data. It is compressed using bzip2 (hence, the .bz2 ending) in a tarball (.tar). To decompress it, you can use the command:

    tar -xjvf human_rnaseq.tar.bz2

    -x = extract
    -j = file has been compressed with bzip
    -v = be verbose. tell me what you're extracting as you do it
    -f = the file to extract

First, you need to align RNAseq reads to the genome using Tophat2. Tophat2 requires that the reference genome be indexed for quick searching. We can index the genome using "bowtie2-build". Make sure you're on an interactive node before you begin. Normally I would submit this command using a submission script to the queue, but since this is so small, it won't take more than a minute or two. Look at the help page for bowtie2-build:

    /usr/local/bowtie2/latest/bin/bowtie2-build -h

    /usr/local/bowtie2/latest/bin/bowtie2-build chr19.fa chr19.fa

Now that we have an indexed genome, we can run Tophat2 to align our reads to the genome. Remember that Tophat2 is splice-aware, meaning that it will "break" a read to map across exon/intron boundaries when aligning to a genome.

First we check the GACRC page to see if we need to load some special libraries (hint: we do). Read the manual to understand which arguments go where. 

https://wiki.gacrc.uga.edu/wiki/Tophat

We're going to align both the Brain and the Adrenal tissues separately. Go ahead and make a submission script, call it whatever you like (I named mine "run_tophat_brainadrenal.sh), and format it like this:

    #!/bin/bash
    export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
    /usr/local/tophat/latest/bin/tophat2 -o adrenal_tophat chr19.fa Adrenal_1.fq Adrenal_2.fq 
    /usr/local/tophat/latest/bin/tophat2 -o brain_tophat chr19.fa Brain_1.fq Brain_2.fq 

Submit it to the rcc-30d queue using 

    qsub -q rcc-30d -cwd ./run_tophat_brainadrenal.sh

Aligned reads are output in a format called SAM/BAM. SAM (Sequence Alignment/Map) is a uniform and accepted format to output read alignment locations and quality scores. BAM is a compressed (binary) format of SAM, which will be smaller in size. Read up on SAM/BAM format and how data is stored here : http://genome.sph.umich.edu/wiki/SAM

Let's look at the end of one of these alignment files. SAMtools is a critically important set of tools that let you view and manipulate SAM/BAM alignment files. "samtools view" is a command that can convert .bam files to .sam very quickly, so that they are human readable. See how I piped the output into the "tail" command, just so I could see the last few lines of it?

    samtools view adrenal_tophat/accepted_hits.bam | tail

Now we can use these two files of aligned reads, one per condition, to identify differentially expressed genes using cuffdiff. 

    /usr/local/cufflinks/latest/bin/cuffdiff USCS_hg19_chr19.genes.gtf brain_tophat/accepted_hits.bam adrenal_tophat/accepted_hits.bam

Cuffdiff produces a variety of output:

    -rw-rw---- 1 aharkess jlmlab 114K Oct 19 15:54 var_model.info
    -rw-rw---- 1 aharkess jlmlab 190K Oct 19 15:56 isoform_exp.diff
    -rw-rw---- 1 aharkess jlmlab 135K Oct 19 15:56 tss_group_exp.diff
    -rw-rw---- 1 aharkess jlmlab 115K Oct 19 15:56 gene_exp.diff
    -rw-rw---- 1 aharkess jlmlab 139K Oct 19 15:56 cds_exp.diff
    -rw-rw---- 1 aharkess jlmlab 134K Oct 19 15:56 splicing.diff
    -rw-rw---- 1 aharkess jlmlab 114K Oct 19 15:56 promoters.diff
    -rw-rw---- 1 aharkess jlmlab  96K Oct 19 15:56 cds.diff
    -rw-rw---- 1 aharkess jlmlab 213K Oct 19 15:56 isoforms.fpkm_tracking
    -rw-rw---- 1 aharkess jlmlab 147K Oct 19 15:56 tss_groups.fpkm_tracking
    -rw-rw---- 1 aharkess jlmlab 153K Oct 19 15:56 cds.fpkm_tracking
    -rw-rw---- 1 aharkess jlmlab 128K Oct 19 15:56 genes.fpkm_tracking
    -rw-rw---- 1 aharkess jlmlab  84K Oct 19 15:56 isoforms.count_tracking
    -rw-rw---- 1 aharkess jlmlab  57K Oct 19 15:56 tss_groups.count_tracking
    -rw-rw---- 1 aharkess jlmlab  58K Oct 19 15:56 cds.count_tracking
    -rw-rw---- 1 aharkess jlmlab  47K Oct 19 15:56 genes.count_tracking
    -rw-rw---- 1 aharkess jlmlab 147K Oct 19 15:56 isoforms.read_group_tracking
    -rw-rw---- 1 aharkess jlmlab  99K Oct 19 15:56 tss_groups.read_group_tracking
    -rw-rw---- 1 aharkess jlmlab  98K Oct 19 15:56 cds.read_group_tracking
    -rw-rw---- 1 aharkess jlmlab  81K Oct 19 15:56 genes.read_group_tracking
    -rw-rw---- 1 aharkess jlmlab  207 Oct 19 15:56 read_groups.info
    -rw-rw---- 1 aharkess jlmlab  206 Oct 19 15:56 run.info
    -rw-rw---- 1 aharkess jlmlab   53 Oct 19 15:56 bias_params.info


A full explanation of each of these files can be found in the manual (http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/#cuffdiff-output-files). Differentially expressed genes at a False Discovery Rate (FDR) < 0.05 are found in the "gene_exp.diff" file. The last column in this file, titled "significance" is a yes/no boolean where "yes" means the gene is differentially expressed, and "no" means the gene is NOT differentially expressed between the two conditions. 

Use grep to identify the differentially expressed gene. (Hint: you are searching for a match to the word "yes". Use Google to learn how to use grep to extract lines that have a match for some string)


## Homework (due 10/27)

1) What gene is differentially expressed between brain and adrenal tissue?

2) What is the fundamental flaw in this reference-based differential expression analysis, and how would you have redesigned the experiment? Be specific.


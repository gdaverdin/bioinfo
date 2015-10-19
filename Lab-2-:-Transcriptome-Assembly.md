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


First, you need to align RNAseq reads to the genome using Tophat2. Tophat2 requires that the reference genome be indexed for quick searching. We can index the genome using "bowtie2-build". Make sure you're on an interactive node before you begin. Normally I would submit this command using a submission script to the queue, but since this is so small, it won't take more than a minute or two. Look at the help page for bowtie2-build:

    /usr/local/bowtie2/latest/bin/bowtie2-build -h

    /usr/local/bowtie2/latest/bin/bowtie2-build chr19.fa chr19.fa

Now that we have an indexed genome, we can run Tophat2 to align our reads to the genome. Remember that Tophat2 is splice-aware, meaning that it will "break" a read to map across exon/intron boundaries when aligning to a genome.

First we check the GACRC page to see if we need to load some special libraries (hint: we do). Read the manual to understand which arguments go where. 

https://wiki.gacrc.uga.edu/wiki/Tophat

We're going to align both the Brain and the Adrenal tissues separately. 

    export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
    /usr/local/tophat/latest/bin/tophat2 -o adrenal_tophat chr19.fa Adrenal_1.fq Adrenal_2.fq 
    /usr/local/tophat/latest/bin/tophat2 -o brain_tophat chr19.fa Brain_1.fq Brain_2.fq 

Aligned reads are output in a format called SAM/BAM. SAM (Sequence Alignment/Map) is a uniform and accepted format to output read alignment locations and quality scores. BAM is a compressed (binary) format of SAM, which will be smaller in size. Read up on SAM/BAM format and how data is stored here : http://genome.sph.umich.edu/wiki/SAM

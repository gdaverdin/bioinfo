# Lab 1 : Bacterial genome assembly

There has been a deadly outbreak of a new strain of the Vibrio cholerae bacterium, known to cause cholera. 

A typical V. cholerae genome is organized into two circular chromosomes with a total length of about 4Mbp (4 million base pairs) with 3,885 annotated genes. 

### Adventure 1: De novo assemble the cholera genome 

Massive amounts of Illumina, Pacbio, 454, Sanger, and other data are stored in the Sequence Read Archive (SRA). Whenever you publish a paper that generates sequence data, you should always submit it to a public repository like SRA. The SRA developers maintain a set of tools to quickly let you download sequence data submitted to the SRA.
	First, download a subset of paired end Illumina whole genome shotgun V. cholerae reads generated at the Center for Disease Control (CDC) using the SRA “fastq-dump” program. 

    /usr/local/sra/latest/bin/fastq-dump --split-files ERR632095

When downloading paired-end data, the --split-files flag separates the forward and reverse reads into two fastq files. If you don't use that flag, you get a single interleaved fastq file.

Let’s first look at the first few lines of our fastq file using “head”. Look at the headers, spot the quality scores. Remember that each read occupies 4 lines of a fastq file.

    head SRR787507.fastq

Let's first check out the quality of this data using FASTQC. 
Hybrid Illumina + Pacbio genome assembly using SPAdes
	SPAdes is a 


https://wiki.gacrc.uga.edu/wiki/SPAdes

    #!/bin/bash
    export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
    python2.7 /usr/local/spades/latest/bin/spades.py --pe1-1 ERR632095_1.fastq --pe1-2 ERR632095_2.fastq --pacbio Vcholerae_ElTor.pacbio.fastq --threads 2 -m 12 -o Vcholera_spades

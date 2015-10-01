There has been a deadly outbreak of a new strain of the Vibrio cholerae bacterium, known to cause cholera. You have been contracted by the Center for Disease Control (CDC) to _de novo_ assemble and annotate this novel, rapidly evolving strain that is causing mass panic as the infection spreads. The CDC has generated whole genome shotgun data using Illumina paired-end 50nt reads (PE50) and longer Pacbio reads, but a government shutdown has forced all bioinformaticians to turn off laboratory computers.

All you know is that a typical _V. cholerae_ genome is organized into two circular chromosomes with a total length of about 4Mbp (4 million base pairs) with ~3,800 annotated genes. 

![source: Wikipedia](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9d/Cholera_bacteria_SEM.jpg/240px-Cholera_bacteria_SEM.jpg)

### Downloading data from SRA
Massive amounts of Illumina, Pacbio, 454, Sanger, and other data are stored in the Sequence Read Archive (SRA). Whenever you publish a paper that generates sequence data, you should always submit it to a public repository like SRA so that it safely remains in the public domain. The SRA developers maintain a set of tools to quickly let a user download sequence data that is archived in the SRA. First, download a subset of paired end 50nt Illumina whole genome shotgun V. cholerae reads generated at the Center for Disease Control using the SRA “fastq-dump” program. Go ahead and log onto an interactive node, too, since we don't want to crush the head node.

    qlogin
    mkdir Vcholerae_genome   # make a new directory for this project
    cd Vcholerae_genome/
    /usr/local/sra/latest/bin/fastq-dump --split-files ERR632095

When downloading paired-end data, the --split-files flag separates the forward and reverse reads into two fastq files. If you don't use that flag, you get a single interleaved fastq file.

Let’s first look at the first few lines of our fastq file using “head”. Look at the headers, spot the quality scores. Remember that each read occupies 4 lines of a fastq file. To save space, the SRA automatically reformats and shortens the read IDs of submitted sequences (@ERR632095.1, @ERR632095.2, @ERR632095.3, ...)

    head ERR632095_1.fastq

### Check the overall quality of the Illumina reads

Now let's look at these reads more objectively using FASTQC. We are looking for any widespread issues of adapter contamination or poor sequence quality. There will usually be some, but less is better.

https://wiki.gacrc.uga.edu/wiki/FastQC

You can run FASTQC on this small dataset on the interactive node as opposed to writing a submission script and submitting a job. It doesn't matter for really small jobs like this.

    /usr/local/fastqc/latest/fastqc ERR632095_1.fastq ERR632095_2.fastq

Download the two html output files (ERR632095_1_fastqc.html and ERR632095_2_fastqc.html) and explore them. Do you see evidence of adapter read-through? A high percentage of poor quality reads?

![](http://i.imgur.com/wUPJKGT.png)

### Hybrid Illumina + Pacbio genome assembly using SPAdes

SPAdes is a cutting-edge genome assembler that specializes in leveraging multiple data types to assemble smaller genomes. Here we will use both Pacbio and Illumina data together to build an assembly for the cholera genome. 

https://wiki.gacrc.uga.edu/wiki/SPAdes

Just like with the Illumina data, we want to first check the quality of the Pacbio data. I have combined SRA data from several SMRT cells into one fastq file for you; copy this fastq file into your working directory. Run FASTQC again on this single fastq file of Pacbio reads.

    /home/student/binf4550/data/01.GenomeAssembly/Vcholerae_ElTor.pacbio.fastq

In your command-line text editor of choice (nano, pico, vim, emacs, etc), create a bash shell submission script. I named mine "run_spades.sh". The *.sh ending is not necessary, but it is good practice, since we are writing this simple script in the bash language. Practice and become proficient at editing text on the command line!

The GACRC will include an example submission script for every program installed on the cluster. There are always some differences in the way programs run and what they require in terms of dependencies. If I had not read the GACRC wiki page for SPAdes, I would NOT have remembered to include the "export LD_LIBRARY_PATH" line and the program would fail immediately. Always, always check the GACRC wiki page first. 

run_spades.sh:

    #!/bin/bash
    export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
    python2.7 /usr/local/spades/latest/bin/spades.py --pe1-1 ERR632095_1.fastq --pe1-2 ERR632095_2.fastq --pacbio Vcholerae_ElTor.pacbio.fastq --threads 2 -m 12 -o Vcholera_spades

Then submit your job by typing this on the command line:

    qsub -q rcc-30d -pe thread 2 ./run_spades.sh

qsub is a command that submits a job to the cluster to be prioritized and run on a node. There are a number of flags you can use with qsub.

    -q = the queue you want to use. rcc-30d is a high volume, all-purpose set of nodes
    -pe thread 2 = specificies that we want to multithread and reserve 2 cores
    ./run_spades.sh = the bash submission script we want to submit. This packet of information contains all the code and information we want to run on the cluster. 
 

This assembly took me 42 minutes to finish. Before you leave lab, you should get this assembly started. 

To calculate some simple and quick statistics on the scaffolds.fasta file (N50, Total length, GC%) I have placed a perl script on my lab server. To fetch it, you can use the wget command. wget is a nifty way to download files from the internet onto the cluster. 

    wget http://jlmwiki.plantbio.uga.edu/~aharkess/calculate_N50.pl

You can look at the code (note: this is ugly code) using less or cat, then execute it with

    perl calculate_N50.pl scaffolds.fasta

### Bonus Points: Bacterial genome annotation using Glimmer

An assembled genome isn't very valuable to us without a set of gene annotations, though. To identify where the genes are in your assembly, we can use the BASys webserver or an NCBI portal that runs Glimmer. In short, Glimmer utilizes Hidden Markov Models (HMMs) from a "training" set of genes to ideally pick the positions of all genes in a query genome. When your assembly finishes, download the scaffolds.fasta file and start an annotation run on the BASys server using default settings.

# Homework (due XX)

1) Assuming a 4.0 megabase (Mb) _V. cholerae_ genome, calculate the approximate coverage of Illumina data that you downloaded and used in the assembly (example: nearly 22X coverage).

2) Run the calculate_N50.pl script. What is your assembled scaffold N50? 
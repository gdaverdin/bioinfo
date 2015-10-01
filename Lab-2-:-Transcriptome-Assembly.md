Today's lab will be a "choose your own adventure" based on your personal research needs. You can pick between 

* building a _de novo_ transcriptome assembly, annotation, and quantification, or
* aligning RNAseq reads to the _Arabidopsis thaliana_ genome, and finding differential expression between two conditions


## _De novo_ transcriptome assembly

You have formed a collaboration with an orchid biologist who studies the "bee orchid", _Ophrys apifera_, which is an incredible example of mimicry. The bee orchid mimics the coloration and scent of a bee species that pollinates it. Male bees will attempt to copulate with the orchid, confusingly thinking that it is a female bee. The male bees unsuspectingly collect pollen on their bodies, move on to another flower, attempt to copulate it, and ultimately pollinate it.
![](https://thmcf.files.wordpress.com/2013/06/bee-orchid-imc-3702.jpg) 

For the sake of time, you will not assemble a transcriptome today. Here is a walkthrough of how I assembled one for you, though.

* Download paired-end 75nt RNAseq reads from SRA


    /usr/local/sra/latest/bin/fastq-dump --split-files SRR609403


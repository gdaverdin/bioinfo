Today we will use R to visualize RNAseq-based gene expression using heatmaps. R is a statistical computing language that is open source and free (unlike SAS). Today you will be generating a heatmap of TMM-normalized FPKM values from an RNAseq differential expression paper I published [(Harkess et al. 2015)](https://www.dropbox.com/s/0s1j2hhhx88n7qs/Harkess_et_al-2015-New_Phytologist-2.pdf?dl=0). 

You will not be using the cluster today. Instead, open R on your computers. 

# Install packages in R

Base or "vanilla" R is plenty powerful and useful, but the benefit of having an open source platform is that anyone can write packages that have specific, improved functions in them. R has a few built-in ways of downloading and installing packages, which makes this extra useful.

My favorite package for plotting data is ggplot2. It is written by Hadley Wickham, a notable R guru. He also develops and maintains a few other packages, most known by RStudio (an R GUI). We'll use two packages today, ggplot2 and reshape. ggplot2 is the ultimate package for developing publication quality figures. reshape is a very nifty package for "reshaping" data into different structures.

    install.packages("ggplot2")
    install.packages("reshape")

# Download data

I have left a copy of the TMM-normalized FPKMs for a set of genes here. Open up your terminal, change directories to the desktop, and make a directory called "heatmaps". cd into this directory, then use wget to download the data to that directory.

    http://jlmwiki.plantbio.uga.edu/~aharkess/Lim_Italian.counts.matrix.TMM_normalized.DEgenes.FPKM

# Read the data into R

    setwd('

    dat <- read.table('Lim_Italian.counts.matrix.TMM_normalized.DEgenes.FPKM', header=TRUE)
    dat.m <- melt(dat)

    p <- ggplot(dat.m, aes(x=variable, y=value)) + geom_boxplot() + coord_flip()

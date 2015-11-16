Today we will use R to visualize RNAseq-based gene expression using heatmaps. R is a statistical computing language that is open source and free (unlike SAS). Today you will be generating a heatmap of TMM-normalized FPKM values from an RNAseq differential expression paper I published [(Harkess et al. 2015)](https://www.dropbox.com/s/0s1j2hhhx88n7qs/Harkess_et_al-2015-New_Phytologist-2.pdf?dl=0). 

You will not be using the cluster today. Instead, open R on your computers. 

# Install packages in R

Base or "vanilla" R is plenty powerful and useful, but the benefit of having an open source platform is that anyone can write packages that have specific, improved functions in them. R has a few built-in ways of downloading and installing packages, which makes this extra useful.

My favorite package for plotting data is ggplot2. It is written by Hadley Wickham, a notable R guru. He also develops and maintains a few other packages, most known by RStudio (an R GUI). 

    install.packages("ggplot2")

# Download data

I have left a copy of the TMM-normalized FPKMs for a set of genes here. 

    http://jlmwiki.plantbio.uga.edu/~aharkess/Lim_Italian.counts.matrix.TMM_normalized.DEgenes.FPKM

# 
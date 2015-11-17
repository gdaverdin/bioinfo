Today we will use R to visualize RNAseq-based gene expression using heatmaps. R is a statistical computing language that is open source and free (unlike SAS). Today you will be generating a heatmap of TMM-normalized FPKM values from an RNAseq differential expression paper I published [(Harkess et al. 2015)](https://www.dropbox.com/s/0s1j2hhhx88n7qs/Harkess_et_al-2015-New_Phytologist-2.pdf?dl=0). 

You will not be using the cluster today. Instead, open R on your computers. 

One caveat you should keep in your mind -- R uses 1-based numbering, meaning that it starts counting at 1 instead of 0. Makes a lot of sense when you're doing math...!

# Open R

All of the lab computers have R installed on them. Open it up, and you'll find an R terminal window opens. If you press command+N, you will open up a blank script that you can enter commands in (and later save them to a file). You can either enter commands directly into the terminal window, or use command+enter to run one or more lines of your script.  

# Install packages in R

Base or "vanilla" R is plenty powerful and useful, but the benefit of having an open source platform is that anyone can write packages that have specific, improved functions in them. R has a few built-in ways of downloading and installing packages, which makes this extra useful.

My favorite package for plotting data is ggplot2. It is written by Hadley Wickham, a notable R guru. He also develops and maintains a few other packages, most known by RStudio (an R GUI). We'll use twothree packages today, ggplot2, reshape, and gplot. ggplot2 is the ultimate package for developing publication quality figures. reshape is a very nifty package for "reshaping" data into different structures. gplot is an additional package which I use for producing heatmaps of data frames. 

    install.packages("ggplot2")
    install.packages("reshape")
    install.packages("gplot")

# Download data

I have left a copy of the TMM-normalized FPKMs for a set of genes here. Open up your terminal, change directories to the desktop, and make a directory called "heatmaps". cd into this directory, then use wget to download the data to that directory.

    http://jlmwiki.plantbio.uga.edu/~aharkess/Lim_Italian.counts.matrix.TMM_normalized.DEgenes.FPKM

# Read the data into R

Well. R is a new language for many of you, and you won't learn its intricacies in one day. But, like many informatic things, it's best to have a project in mind and jump right in. I'll comment the code so you can see what its doing. Remember that comments (lines starting with #) don't get executed. To bring up the help/manual page for a function, use a question mark in front of the command, like this:

    ?read.table

Alright, let's start by reading our data in. R can parse data tables using the read.table() function. Since our data table has a header line with the library names, I tell R to name each column by that header using the "header=TRUE" boolean. We also assign the output of the read.table command to save to a variable called "dat". See how I used the <- arrow to create a new variable?

    dat <- read.table('Lim_Italian.counts.matrix.TMM_normalized.DEgenes.FPKM', header=TRUE)
    
Now let's look at our data. Just like in unix, let's use head to look at the first few lines.

    head(dat)

# Manipulate and reshape the data

The melt() function in the reshape package is infinitely useful. It can take a data frame and guess the best way to reformat it into a more simple structure. Let's use the melt function, and save the output to a new variable called dat.m

dat.m <- melt(dat)
head(dat.m)

# 
    p <- ggplot(dat.m, aes(x=variable, y=value)) + geom_boxplot() + coord_flip()

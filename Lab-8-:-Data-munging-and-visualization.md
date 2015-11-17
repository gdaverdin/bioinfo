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
    install.packages("gplots")

and load them by using the library() command

    library(ggplot2)
    library(reshape)
    library(gplots)

# Download data

I have left a copy of the TMM-normalized FPKMs for a set of genes here. Open up your terminal, change directories to the desktop, and make a directory called "heatmaps". cd into this directory, then use wget to download the data to that directory.

    http://jlmwiki.plantbio.uga.edu/~aharkess/TMM_normalized_FPKM_matrix.txt

# Read the data into R

Well. R is a new language for many of you, and you won't learn its intricacies in one day. But, like many informatic things, it's best to have a project in mind and jump right in. I'll comment the code so you can see what its doing. Remember that comments (lines starting with #) don't get executed. To bring up the help/manual page for a function, use a question mark in front of the command, like this:

    ?read.table

Alright, let's start by reading our data in. R can parse data tables using the read.table() function. Since our data table has a header line with the library names, I tell R to name each column by that header using the "header=TRUE" boolean. We also assign the output of the read.table command to save to a variable called "dat". See how I used the <- arrow to create a new variable?

First, set your working directory to the heatmaps directory you created and downloaded the data to.

    setwd('/path/to/heatmaps/folder')

    dat <- read.table('TMM_normalized_FPKM_matrix.txt', header=TRUE)
    
Now let's look at our data. Just like in unix, let's use head to look at the first few lines.

    head(dat)

Let's figure out how many genes are in this set by calculating the dimensions of our data frame. The output is "rows   columns"

    dim(dat)
    [1] 570   9
# Manipulate and reshape the data

The melt() function in the reshape package is infinitely useful. It can take a data frame and guess the best way to reformat it into a more simple structure. Let's use the melt function, and save the output to a new variable called dat.m

    dat.m <- melt(dat)
    head(dat.m)

# Plot the reshaped data using ggplot2 to generate a boxplot and heatmap for each sample

    p <- ggplot(dat.m, aes(x=variable, y=value)) + geom_boxplot() + coord_flip()
    p

Kinda hard to see anything with such a large Y axis. Let's limit this to genes with FPKMs < 100
   
    p + ylim(0,100)

![](http://i.imgur.com/dfTkN2V.png)

Now let's make a heatmap. It helps to define your own color palette instead of using defaults, though. 

heatmap.2 requires that we convert our data frame into a matrix. This luckily is easy.

    dat <- as.matrix(dat)

    my_palette <- colorRampPalette(c("green","black","red"))(n = 1000)
    heatmap.2(dat, dendrogram="both", trace="none", scale="row", density.info="none", col=my_palette)

![](http://i.imgur.com/ovAz600.png)
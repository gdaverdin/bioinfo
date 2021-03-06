For your final project for this class, you have quite a bit of freedom, but the foundation must including generating a genome assembly, transcriptome assembly, variant identification or data wrangling and visualization. This will culminate in a short presentation (10-15min) describing 

* your experimental design and question,
* your data source and analysis pipeline, 
* any pitfalls you encountered
* a summary of your results and conclusions

You must also turn in a minimum 6 page paper (double spaced) detailing these points. The paper should be formatted like a journal article: Introduction, Methods, Results, Discussion. **Undergraduates are exempt from the written paper.**

If you would like, you can work with a partner on the analysis and interpretation of your data, and present jointly for the final project. If you work in a group, you must clearly articulate a specific analysis that you were in charge of. You must **independently** write your own 6 page papers. Please ask me if you have any questions regarding potential plagiarism. 

# Potential Projects:

* Project 1: Use SPAdes and Velvet to assemble Illumina reads for an _E. coli_ genome that produces Shiga toxin. This sample was isolated from faeces of diarrhetic cattle and sequenced using Miseq PE300 reads. Annotate your genome with the NCBI Glimmer portal, and use Bedtools to extract a multifasta file of every gene location. Use BLAST to compare your gene annotations to that of another E. coli strain of your choice and plot a distribution of the % similarity for each best hit. http://www.ncbi.nlm.nih.gov/sra/SRX1368955[accn]

* Project 2: Assemble the Salmonella enterica subsp. enterica serovar Typhimurium genome, a highly pathogenic serovar that was used as a biological weapon in the [1984 Rajneeshee bioterror attack](https://en.wikipedia.org/wiki/1984_Rajneeshee_bioterror_attack). These are Miseq PE250 reads. Annotate your genome with the NCBI Glimmer portal, then compare the quality of your de novo assembly to another Salmonella genome assembly using [Quast](http://bioinf.spbau.ru/quast).
http://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=2018078 

* Project 3: Use [QIIME](https://wiki.gacrc.uga.edu/wiki/Qiime) to characterize alpha and beta diversity in a heterotrophic bacterial community assembly in an activated sludge wastewater treatment plant.  http://www.ncbi.nlm.nih.gov/sra/SRX669500[accn]

* Project 4: _De novo_ assemble the bee orchid transcriptome, _Ophrys apifera_ (SRA id = SRR609403) using Trinity. Compare your assembled transcripts against Uniprot genes to identify how many are putative full length assemblies by following the [Trinity manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts). 

* Project 5: Use R and the [rWBclimate package](https://ropensci.org/tutorials/rwbclimate_tutorial.html) to scrape model data from the World Bank Climate Data to visualize the temperature changes in Africa that are predicted under various climate models for the next 40 years. 



If you have your own data, you are encouraged to use that. We are here to help!
Mothur_to_PhyloSeq_to_Differential_Abundance_Tools
==================
An R-script providing guided import from Mothur data to PhyloSeq followed by profiling for differential abundance according to a given factor. This script was developed for application to identifying OTUs enriched in a SIP-DNA experiment.

Note: This is one of the author's first attempts at writing a full analysis pipeline AND, he admits that it is not written in a legible or efficient manner. A large proportion of the script could have been written as a loop, and it is not adequately commented. Hopefully the latter will be remedied shortly.


Block00 - Installs all required libraries/packages  (Universal)

Block01 - Guides user to import all necessary files for doing analysis with PhyloSeq (for manipulating and performing ecological analyses of community data). The scripts also creates extensive back-up .rds files, since it is generally slow to import large scale sequencing datasets.	(Universal)

Block02 - Guides user through determining differential abundance analysis using DESeq2, Limma-voom, or a custom-made relative abundance method. (Universal, but some of the language is geared towards stable isotope probing, i.e. the main factor of interest is "Enrichment", but this can be substituted for whatever you are interested in doing differential analyses on)

Block03 - Not automated - Just some validation steps I used to compare the three methods of identifying OTUs indentified as differentially abundant.

Block05 - A script to generate a list of fasta headers with taxonomic information for downstream submission or further sequence-based analysis (note: I've included the Python script "grab_ALL_READS_from_mothur_list_file.py" I use for recovering sequence information for each OTU. There is usage information included within it)

Feel free to e-mail the author with any questions. It has had two data sets run through it and is reasonably well debugged. The author is novice-to-intermediate, so some portions of code may be slower or could use improvement.


The actual PhloSeq project GitHub can be found at the following address.

Visit : https://github.com/joey711/phyloseq
==================

The code was used in the following publication(s)...

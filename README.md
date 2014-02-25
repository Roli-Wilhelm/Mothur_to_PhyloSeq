Mothur_to_PhyloSeq
==================
An R-script providing guided import from Mothur data to PhyloSeq

This code was used to simplify the import of "count data" output from the sequence processing pipeline Mothur.


Block00 - Installs all required libraries/packages  (Universal)

Block01 - Guides user to import all necessary files for doing analysis with PhyloSeq (for manipulating and performing ecological analyses of community data). The scripts also creates extensive back-up .rds files, since it is generally slow to import large scale sequencing datasets.	(Universal)

Block02 - Guides user through determining differential expression according to some treatment factor using DESeq2 (for examining differential expression) (Universal, but some of the language is geared towards stable isotope probing, i.e. the main factor of interest is "Enrichment", but this can be substituted for whatever you are interested in doing differential analyses on)

Feel free to e-mail the author with any questions. It has had two data sets run through it and is reasonably well debugged. The author is novice-to-intermediate, so some portions of code may be slower or could use improvement.


The actual PhloSeq project GitHub can be found at the following address.

Visit : https://github.com/joey711/phyloseq
==================

The code was used in the following publication(s)...

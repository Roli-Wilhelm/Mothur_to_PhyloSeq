library(phyloseq)
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

PATH <- getwd()
dir.create(paste(PATH, "\\Output", sep=""))
dir.create(paste(PATH, "\\Output\\BackUps", sep=""))

writeLines("\nThe minimum files required to create a phyloSeq object from mothur are : \".list\", \".groups\", \".taxonomy\" and a \".csv\" indicating which factors correspond to your samples.\n\nWe will begin by inputting all of this information.")

writeLines("\nPlease provide a one word descriptor of the library you are about to process.\nExample : \"Bacteria\"")
sign_save<-scan(n=1, what = character())

###Importing file into phyloSeq
writeLines("\nChoose your \".list\" file.")
mothlist<-file.choose()

writeLines("\nChoose your \".groups\" file.")
mothgroup<-file.choose()

writeLines("\nDo you wish to import a phylogenetic tree file (\".tre\") file?")
switch(menu(c("Yes", "No")), tree<-("Yes"), tree<-("No"))

if (tree == "Yes") {

  writeLines("\nChoose your phylogenetic tree file \".tre\" file. If you have none, ENTER: None")
  mothtree<-file.choose()

}

writeLines("\n\nThese are the following cut-offs detected in your \".list\" file:")
print(show_mothur_list_cutoffs(mothlist))
writeLines("\nEnter the dissimilarity cut-off you'd like to base your analysis on. Note: this must be present in your \".list\" file.\nExample : 0.05\n")
cutoff<-scan(n=1, what = character())

writeLines("\n\nPlease note: The following step may take upwards of an hour to import your \".list\" file. Grab a coffee; do something else. If you have already done this step, and have had to re-start the script, please answer accordingly.\n\n Have you specified then imported (long step) your \".list\" and \".groups\" files?")
switch(menu(c("Yes", "No")), redo<-("Yes"), redo<-("No"))

##Option of skipping

if (redo == "No") {
  
###Creating "OTU table" phyloSeq object
  if (tree == "No") {
    x <- import_mothur(mothlist, mothgroup, mothur_tree_file=NULL, cutoff)
    x
    saveRDS(x, file=paste("Output\\BackUps\\OTU_Table_", sign_save, "_", cutoff,".rds", sep=""))

  } else {
    x <- import_mothur(mothlist, mothgroup, mothtree, cutoff)
    x
    saveRDS(x, file=paste("Output\\BackUps\\OTU_Table_wTree_", sign_save, "_", cutoff,".rds", sep=""))
  }
} else {
  
  x<-readRDS(paste("Output\\BackUps\\OTU_Table_", sign_save, "_", cutoff,".rds", sep=""))
  
}

###Create Taxonomy table
writeLines("\nTo input a taxonomy table from Mothur, you must first manually change all tabs to \";\", because R cannot read in two different delimiters. Do this on a linux system with sed.\n Example : sed -i 's/\t/;/g' yourtax.tax.\nYou will also have to remove \";\" found at the end of your line with sed -i 's/;$//g'.\n\nAgain, if you are re-running this script, you are given the opportunity to skip importing your taxonomy file. Note: You would not have to re-import your taxonomy file if you have imported your OTU table at a new cutoff.")
writeLines("\nWould you like to skip importing your taxonomy table and use the one that exists in your \"Output\" folder?")
switch(menu(c("Yes", "No")), redo<-("Yes"), redo<-("No"))

if (redo == "No") {
  
  tax <- as.matrix(read.table(file.choose(), sep = ";", header = FALSE))
  tax <- tax[,colSums(is.na(tax))<nrow(tax)] #Remove accidental NA-filled columns
  saveRDS(tax, file=paste("Output\\BackUps\\Tax_Table_", sign_save, "_", cutoff,".rds", sep=""))
  
  ###Subset Taxonomy data to matching read names 
  #This step is necessary because the otu table now contains a sub-set of OTUs based on whichever OTU cut-off you used. 
  #It will not be in synch with the newly input taxonomy table and this will cause errors.
  OTU_overlap<-intersect(tax[,1], taxa_names(x))

  #Create empty matrix for containing intersection
  tax_adj <- matrix(ncol = ncol(tax), nrow = length(OTU_overlap))

  for (i in 1:length(OTU_overlap)) {
  
    tax_adj[i,] <- tax[grep(paste("^",OTU_overlap[i],"$",sep=""), tax[,1]),]
  
  }

  #Order matrices according to read names in the phyloseq objects (had to do this b/c of issues of losing the phyloseq class with otu_table and phy_tree present and "import_mothur_tree" not working independently.)
  Neworder<- vector(mode="numeric", length=nrow(tax_adj))

  for (mercy in 1:nrow(tax_adj)) {

     Neworder[mercy]<-grep(taxa_names(otu_table(x))[mercy], tax_adj[,1])
  
  }

  tax_adj<-tax_adj[Neworder,]

  #Rename rownames of tax table
  rownames(tax_adj) <- taxa_names(x)

  #Remove sample name column in tax_adj (rows are already named)
  tax_adj<-tax_adj[,-1]

  #Save the object for ease of re-running the script

  saveRDS(tax_adj, file=paste("Output\\BackUps\\Tax_Table_adj", sign_save, "_", cutoff,".rds", sep=""))

} else {
  
  writeLines("\nAre you re-running the script at a new cut-off or just re-running the script b/c there was an error?")
  switch(menu(c("I'm re-running at a new cut-off", "I'm just re-running the script")), redo<-("Yes"), redo<-("No"))
  
    if (redo == "Yes") {
      
      tax<-readRDS(paste("Output\\BackUps\\Tax_Table_", sign_save, "_", cutoff,".rds", sep=""))
      ###Subset Taxonomy data to matching read names 
      #This step is necessary because the otu table now contains a sub-set of OTUs based on whichever OTU cut-off you used. 
      #It will not be in synch with the newly input taxonomy table and this will cause errors.
      OTU_overlap<-intersect(tax[,1], taxa_names(x))
      
      #Create empty matrix for containing intersection
      tax_adj <- matrix(ncol = ncol(tax), nrow = length(OTU_overlap))
      
      for (i in 1:length(OTU_overlap)) {
        
        tax_adj[i,] <- tax[grep(paste("^",OTU_overlap[i],"$",sep=""), tax[,1]),]
        
      }
      
      #Order matrices according to read names in the phyloseq objects (had to do this b/c of issues of losing the phyloseq class with otu_table and phy_tree present and "import_mothur_tree" not working independently.)
      Neworder<- vector(mode="numeric", length=nrow(tax_adj))
      
      for (mercy in 1:nrow(tax_adj)) {
        
        Neworder[mercy]<-grep(taxa_names(otu_table(x))[mercy], tax_adj[,1])
        
      }
      
      tax_adj<-tax_adj[Neworder,]
      
      #Rename rownames of tax table
      rownames(tax_adj) <- taxa_names(x)
      
      #Remove sample name column in tax_adj (rows are already named)
      tax_adj<-tax_adj[,-1]
      
      #Save the object for ease of re-running the script
      
      saveRDS(tax_adj, file=paste("Output\\BackUps\\Tax_Table_adj", sign_save, "_", cutoff,".rds", sep=""))
  
    } else {
    
      tax_adj<-readRDS(paste("Output\\BackUps\\Tax_Table_adj", sign_save, "_", cutoff,".rds", sep=""))
      
  }
}

###Options to read in taxonomic levels
writeLines("\nDo you wish to name the taxonomic hiearachy (i.e. levels) for downstream analysis? \n\nNote: If taxonomy was assigned usin the Silva database, you will not be able to use a standard set of levels b/c the columns are not consistent to any specific taxa level. I recommend using GreenGenes (as of Jan 2014)")
standard_levels<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
print(standard_levels)
switch(menu(c("Yes, I'll input it myself", "Yes, use the standard levels you have listed", "No, thank you")), tax_opt<-("self_input"), tax_opt<-("standard_levels"), tax_opt<-("No"))

if (tax_opt == "self_input") {
 writeLines("\nHow many levels would you like to enter?\nExample: 7")
 level_numb<-scan(n=1, what = numeric())
 writeLines("\nPlease enter your levels one-by-one, one line for each.\nExample: Domain (ENTER) Phylum (ENTER) Sub-phylum (ENTER) ... etc.")
 standard_levels<-scan(n=level_numb, what = character())
 
 ##Must ensure that every column in the taxonomy table has a header, so will artificially create additional species levels if empty columns are present
  if ((ncol(tax_adj)-length(standard_levels)) != 0) {
 
    for (unique in 1:(ncol(tax_adj)-length(standard_levels))) {
   
      standard_levels<-c(standard_levels, paste("Species_level", unique, sep=""))  
   
    }
  } else { 
    
    colnames(tax_adj)<-standard_levels[count]
         
  }
} else {

  if (tax_opt == "standard_levels") {
  
    ##Must ensure that every column in the taxonomy table has a header, so will artificially create additional species levels if empty columns are present
    
    if ((ncol(tax_adj)-length(standard_levels)) != 0) {
      
      for (unique in 1:(ncol(tax_adj)-length(standard_levels))) {
        
        standard_levels<-c(standard_levels, paste("Species_level", unique, sep=""))  
        
      }
    }
    
    colnames(tax_adj)<-standard_levels
      
  } 
}

###Create object class Tax_table using phyloseq command
tax = tax_table(tax_adj)

###Construct a phyloseq-class from OTU table and Taxonomy table
p<-merge_phyloseq(x, tax)

##Remove Samples which are not part of the data set
writeLines("\n\nDo you have samples which you do not wish to analyze present in this set?")
switch(menu(c("Yes","No")), set_rem<-("Yes"), set_rem<-("No"))

if (set_rem == "Yes") {
  writeLines("\n\nWe will use regular expressions to identify which samples you wish to remove. Thus, you can stipulate exact names of samples or portions of a sample name which correspond to a larger set you'd like to remove.\nCan you remove all the samples with a single entry?")
  switch(menu(c("Yes","No")), set_rem2<-("Yes"), set_rem2<-("No"))
  
  if (set_rem == "Yes") {
    writeLines("\n\nPlease enter the exact sample name or pattern you would like to remove from your dataset.\nExample : BL[0-9][0-9][0-9][0-9].27f to remove all samples titled BL0010, BL008, BL102 etc.\nNote: use [0-9] to refer to any number and [a-z] or [A-Z] for any letter. These can be modified to any range.")
    pattern <-scan(n=1, what = character())
    writeLines(paste("Your pattern will remove :", length(grep(pattern, sample_names(p))), "sample(s) from your dataset."))
    p<-prune_samples(colnames(otu_table(p))[-grep(pattern, sample_names(p))], p)

  } else {
    
    writeLines("\n\nHow many separate searches and removals would you like to use?\n")
    searches<-scan(n=1, what = numeric())
    
    for (count in 1:searches) {
      
      writeLines("\n\nPlease enter the exact sample name or pattern you would like to remove from your dataset.\nExample : BL[0-9][0-9][0-9][0-9].27f to remove all samples titled BL0010, BL008, BL102 etc.\nNote: use [0-9] to refer to any number and [a-z] or [A-Z] for any letter. These can be modified to any range.")
      pattern <-scan(n=1, what = character())
      writeLines(paste("Your pattern will remove :", length(grep(pattern, sample_names(p))), "samples from your dataset."))
      p<-prune_samples(colnames(otu_table(p))[-grep(pattern, sample_names(p))], p)
      
    }
  }
}

###Import factor from .csv 
writeLines("\nPlease input a file containing the various factor levels accorded to your samples.\nThe file should be \".csv\" format and the data in a typical form of a design matrix where each row correponds to a sample, and each colum a factor level. Be sure to remove any factors which correspond to samples you removed in previous steps.")
#Display order of samples:
writeLines("\nHere is the current order of your samples. Please list them in your spread sheet in the same order.")
print(data.frame(colnames(otu_table(p))))

#Select file
Factors<-read.csv(file.choose(), header = T)
Factors <- Factors[,colSums(is.na(Factors))<nrow(Factors)] #Remove accidental NA-filled columns

#Transpose sample names to rownames
Factors <- data.frame(Factors, row.names = sample_names(p), stringsAsFactors=F)
Factors<-Factors[,-1]

###Convert into sample_data-class phyloseq object
sample_data <- sample_data(Factors)

#Merge with existing phyloseq object
p<-merge_phyloseq(p, sample_data)


writeLines("\nWould you like to visually inspect your OTU table? (Only first 40 entries showns)")
switch(menu(c("Yes","No")), disp<-("Yes"), disp<-("No"))

if (disp == "Yes") {
  
  print(head(otu_table(p), 15))
  
}

writeLines("\nWould you like to visually inspect the stored sample information table?")
switch(menu(c("Yes","No")), disp<-("Yes"), disp<-("No"))

if (disp == "Yes") {
  
  print(sample_data(p))
  
}

writeLines("\nWould you like to visually inspect your taxonomy table? (Only first 40 entries showns)")
switch(menu(c("Yes","No")), disp<-("Yes"), disp<-("No"))

if (disp == "Yes") {
  
  print(head(tax_table(p), 40))
  
}

###Back-up complete phyloSeq object
saveRDS(p, file=paste("Output\\BackUps\\phyloseq_", sign_save, "_", cutoff,".rds", sep=""))
p<-readRDS(paste("Output\\BackUps\\phyloseq_", sign_save, "_", cutoff,".rds", sep=""))

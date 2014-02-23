library(DESeq2)
library(phyloseq)
library(ggplot2)
library(reshape2)

writeLines("Are you re-running the script starting at Block02?")

switch(menu(c("Yes","No")), rerun<-("Yes"), rerun<-("No"))

if (rerun == "Yes") {
  writeLines(paste("\nPlease re-select your factor list (i.e. design matrix)"))
  Factors<-read.csv(file.choose(), header = T)
  Factors <- Factors[,colSums(is.na(Factors))<nrow(Factors)] #Remove accidental NA-filled columns
  
  #Transpose sample names to rownames
  Factors <- data.frame(Factors, row.names = sample_names(p), stringsAsFactors=F)
  Factors<-Factors[,-1]
  
  writeLines("\nPlease enter the cut-off you have been using (i.e. 0.01)")
  cutoff<-scan(n=1, what = character())
  writeLines("\nPlease enter the distinguishing name you originally used (i.e. \"Bacteria\")")
  sign_save<-scan(n=1, what = character())
  
  p<-readRDS(paste("Output\\BackUps\\phyloseq_", sign_save, "_", cutoff,".rds", sep=""))
  
}

###Filter out rare OTUs and 13C OTUs not present in OTU table to reduce number of OTUs
writeLines("The first step in identifying OTUs will be to reduce the overall number of OTUs, b/c the differential expression package we will use (DESeq2) crashes with OTU tables with > ~ 5000 samples.\n")
writeLines("Your factor list is :")
print(colnames(sample_data(p)))
writeLines("\nPlease supply the factor label you've given to enriched status.\n")
status<-scan(n=1, what = character())
status_list<-grep(status, colnames(sample_data(p)))
unique_status<-which(duplicated(sample_data(p)[,status_list]) == FALSE)
writeLines("Your enrichment status variables are :")
print(sample_data(p)[unique_status,status_list])
writeLines("\nPlease supply the label you used to demarcate 13C enriched samples. Example : C13\n")
label<-scan(n=1, what = character())

#Subset Samples in order to treat C13 sample differently
p_13C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] == paste(label))], p)
p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)

writeLines("\nWe've separated your OTU table based on enrichments status. Please verify this has been done correctly.\n")
print(sample_data(p_13C))
print(sample_data(p_12C))
writeLines("\nIf this is not correct, you'll have to restart the script and re-enter the information properly, OR re-craft your design matrix to successful discriminate between enrichment status.\n")

###Pre-filtering steps
#Determine which OTU are not found in the OTU table (after C12 sample have been removed)
writeLines("\nWe will first remove all OTUs which have zero representation in the C13 group.")
total<-nrow(otu_table(p))
number<-length(which(rowSums(otu_table(p_13C)) == 0))
percent1<-format(round(((number/total)*100), 2), nsmall = 2)
p_13C<-otu_table(p_13C)[-which(rowSums(otu_table(p_13C)) == 0),]

writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))

##Filter based on abundance, maybe rowSum > 5
writeLines("\nNext we will filter by a minimum abundance value across all C13 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C13 samples.\n\n Please select a threshold number. [Recommended: 5] ")
threshold<-scan(n=1, what = numeric())
number<-length(which(rowSums(otu_table(p_13C)) < threshold))
p_13C<-otu_table(p_13C)[-which(rowSums(otu_table(p_13C)) < threshold),]
percent2<-format(round(((number/total)*100), 2), nsmall = 2)
writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))

##Filter based on occurrence in minimum number of samples > 3
writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C13 samples.\n\n Please select a minimum number. [Recommended: 3]")
min_number<-scan(n=1, what = numeric())
present_absent<-otu_table(p_13C)
present_absent[present_absent > 0] <- 1
number<-length(which(rowSums(present_absent) < min_number))
p_13C<-otu_table(p_13C)[-which(rowSums(present_absent) < min_number),]
percent3<-format(round(((number/total)*100), 2), nsmall = 2)
remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C13 dataset."))

#Return to original phyloseq object
p_enr<-intersect(rownames(p_13C), taxa_names(p_12C))
p_enr<-prune_taxa(p_enr, p)

writeLines("\nYour data table now looks like this:\n")
print(p_enr)

### Currently have to normalize counts in DESeq without the phyloseq function phyloseq_to_deseq2
### This is likely to change shortly once they release this functionality

#This step is only made necessary by the fact I have to independently use DESeq, which conflicts with the formatting of a phyloSeq object
rows<-nrow(otu_table(p_enr))
cols<-ncol(otu_table(p_enr))
counts<-as.matrix(otu_table(p_enr)[,])
attributes(counts)<-NULL  #Remove phyloseq parameters
counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
rownames(counts)<-taxa_names(p_enr)   #Return row and column names
colnames(counts)<-sample_names(p_enr)
counts<-data.frame(counts)

#Prepare the condition matrix 
conds<-matrix(nrow=cols, ncol=1) #Rows are equal to columns here b/c were talking about samples not OTUs
conds[which(sample_data(p_enr)[,status] == paste(label))] <- paste(label)
conds[which(sample_data(p_enr)[,status] != paste(label))] <- "C12"
conds<-data.frame(conds, as.vector(Factors[, which(colnames(Factors) != paste(status))]))
rownames(conds) <- colnames(counts)

#Rename column names
colnames(conds)[1]<- paste(status)

#DESeq calculates the geometric mean, which means all the zeros in my columns are problematic
counts<-counts+1

###Create DESeq Object According to Enrichment any other treatment variable
writeLines("\nDifferential expression profiling, (package \"DESeq\") will be performed to identify which OTU are differentially expressed according to the enrichment status you have provided.\nYou may provide multiple factors to DESeq in order to reduce the amount of \'within-group\'-type error and focus on enrichment status.\n Because the DESeq function can not handle calling functions as arguments (i.e. as.symbolic), you will have to manually ammend the DESeq Command in a script window which will open.\n\nHave you completed editing this file?")
file.edit("DESeq Command.R")
switch(menu(c("Yes")), proceed<-("Yes"))

if (proceed == "Yes") {
  
  source("DESeq Command.R")
  
}

###Run Test for Differential Expression
dds <- DESeq(DE)
DE_results<-results(dds)
DE_results<-DE_results[order(DE_results$log2FoldChange),]

writeLines("\nDESeq will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially expressed in favour of C13 enrichment.")

writeLines("\nHere is the count data for your differentially abundant OTU.")
print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])

writeLines("\nIs it showing abundance in your C13 samples?")

switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))

if (order == "No") {
  
  DE_results$log2FoldChange<-DE_results$log2FoldChange*(-1)
  DE_results<-DE_results[order(DE_results$log2FoldChange),]
  writeLines("\nNow it should!")
  print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
  
}

#Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
enr_otu<-setdiff(union(which(DE_results$log2FoldChange < -1.25), which(DE_results$padj < 0.05)), which(DE_results$log2FoldChange > 0))
highlikely<-which(DE_results$log2FoldChange < -1.5)
lowlikely<-intersect(which(DE_results$log2FoldChange < 0 & DE_results$log2FoldChange > -1.5), enr_otu)

enr_otu<-rownames(DE_results)[enr_otu]
highlikely<-rownames(DE_results)[highlikely]
highlikely<-data.frame(highlikely)
colnames(highlikely)<-"OTU"
lowlikely<-rownames(DE_results)[lowlikely]

writeLines(paste("\nIn total you have recovered ", length(enr_otu)," enriched OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be enriched and ", length(lowlikely)," will have to be manually inspected and confirmed."))

###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
d=enr_otu
OTU_location<- vector(mode="numeric", length=length(d))

for (mercy in 1:length(d)) {
  
  OTU_location[mercy]<-grep(enr_otu[mercy], rownames(counts))
  
}

enr_otu<-counts[OTU_location,]

#Return all Count Data to Original status be subtracting One
enr_otu <- enr_otu-1

#Make into data.frame
enr_otu<-data.frame(enr_otu)
enr_otu$likelihood<-NA
OTU=1

#Identify demarcate high and low likelihood OTUs and add in the factor information
for (OTU in 1:nrow(enr_otu)) {
  
  if (length(grep(rownames(enr_otu)[OTU], highlikely$OTU)) == 0) {
    
    enr_otu$likelihood[OTU]<- "Low"  
    
  } else {
    
    enr_otu$likelihood[OTU]<- "High"
    
  }
}    

#Separate high rollers from less clear differentiation
enr_otu_high<-subset(enr_otu, likelihood == "High")
enr_otu_low<-subset(enr_otu, likelihood == "Low")
enr_otu_high$likelihood<- NULL
enr_otu_low$likelihood<- NULL
enr_otu_high<-t(enr_otu_high)
enr_otu_low<-t(enr_otu_low)

#Add in factors
enr_otu_high<-data.frame(enr_otu_high, as.vector(Factors))
enr_otu_low<-data.frame(enr_otu_low, as.vector(Factors))

#Graph (Highly Likely OTUs as boxplots)
melted_low<-melt(enr_otu_low)
melted_high<-melt(enr_otu_high)

writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))

#### Set-up the auto-remove loop
writeLines("\nYou can now visually verify all of the OTUs which are borderline enriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
writeLines(paste("\nWould you like to manually examine your", length(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
switch(menu(c("Visualize and edit","Do it myself later")), pain<-("Yes"), pain<-("No"))

if (pain == "Yes") {
  writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
  writeLines("Your factor list is :")
  print(colnames(sample_data(p)))
  factor_2<-scan(n=1, what = character())

  lowlikely_remove<-data.frame(lowlikely, keep=NA)
  colnames(lowlikely_remove)[1]<-"OTU"
  writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
  for (T in 1:length(lowlikely)) {
  
    graphics.off()
    print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
    
    writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
    print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
    writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
    switch(menu(c("Keep","Discard")), keep<-("Yes"), keep<-("No"))
  
    if (keep == "Yes") {
    
      lowlikely_remove$keep[T]<- "1"
    
    } else {
    
      lowlikely_remove$keep[T]<- "0"
  
    }
  }
}

writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
switch(menu(c("Visualize and edit","Do it myself later")), pain<-("Yes"), pain<-("No"))

if (pain == "Yes") {
  writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
  writeLines("Your factor list is :")
  print(colnames(sample_data(p)))
  factor_2<-scan(n=1, what = character())
  
  highlikely_remove<-data.frame(highlikely, keep=NA)
  colnames(highlikely_remove)[1]<-"OTU"
  writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
  
  for (T in 1:nrow(highlikely)) {
    
    graphics.off()
    print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
    
    writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
    print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
    writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
    switch(menu(c("Keep","Discard")), keep<-("Yes"), keep<-("No"))
    
    if (keep == "Yes") {
      
      highlikely_remove$keep[T]<- "1"
      
    } else {
      
      highlikely_remove$keep[T]<- "0"
      
    }
  }
}


##Keep the keepers alongside the highlikely OTUs
manual<-rbind(highlikely_remove, lowlikely_remove)
manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
write.csv(manual, file=paste("Output\\C13 OTUs - List of Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
enr_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep == 1]), as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep == 1])))

#Prune phyloseq to contain enriched OTUs
p_enr<-prune_taxa(enr_otu, p)

saveRDS(p_enr, file=paste("Output\\BackUps\\p_All_Enriched_OTU_", sign_save, "_", cutoff,".rds", sep=""))
write.csv(taxa_names(p_enr), file= paste("Output\\p_All_Enriched_OTU_", sign_save, "_", cutoff,".csv", sep=""))

writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))

if (proceed == "Yes") {
  
  if (order == "No") {
    
    DE_results$log2FoldChange<-DE_results$log2FoldChange*(-1)
    DE_results<-DE_results[order(DE_results$log2FoldChange),]
    
  }
  
  #Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
  C12_otu<-setdiff(union(which(DE_results$log2FoldChange < -1.25), which(DE_results$padj < 0.05)), which(DE_results$log2FoldChange > 0))
  highlikely<-which(DE_results$log2FoldChange < -1.5)
  lowlikely<-intersect(which(DE_results$log2FoldChange < 0 & DE_results$log2FoldChange > -1.5), C12_otu)
  
  C12_otu<-rownames(DE_results)[C12_otu]
  highlikely<-rownames(DE_results)[highlikely]
  highlikely<-data.frame(highlikely)
  colnames(highlikely)<-"OTU"
  lowlikely<-rownames(DE_results)[lowlikely]
  
  writeLines(paste("\nIn total you have recovered ", length(C12_otu)," C12-only OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be C12-only and ", length(lowlikely)," are marginal."))
  
  ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
  d=C12_otu
  OTU_location<- vector(mode="numeric", length=length(d))
  
  for (mercy in 1:length(d)) {
    
    OTU_location[mercy]<-grep(C12_otu[mercy], rownames(counts))
    
  }
  
  C12_otu<-counts[OTU_location,]
  
  #Return all Count Data to Original status be subtracting One
  C12_otu <- C12_otu-1
  
  #Make into data.frame
  C12_otu<-data.frame(C12_otu)
  C12_otu$likelihood<-NA
  OTU=1
  
  #Identify demarcate high and low likelihood OTUs and add in the factor information
  for (OTU in 1:nrow(C12_otu)) {
    
    if (length(grep(rownames(C12_otu)[OTU], highlikely$OTU)) == 0) {
      
      C12_otu$likelihood[OTU]<- "Low"  
      
    } else {
      
      C12_otu$likelihood[OTU]<- "High"
      
    }
  }    
  
  #Separate high rollers from less clear differentiation
  C12_otu_high<-subset(C12_otu, likelihood == "High")
  C12_otu_low<-subset(C12_otu, likelihood == "Low")
  C12_otu_high$likelihood<- NULL
  C12_otu_low$likelihood<- NULL
  C12_otu_high<-t(C12_otu_high)
  C12_otu_low<-t(C12_otu_low)
  
  #Add in factors
  C12_otu_high<-data.frame(C12_otu_high, as.vector(Factors))
  C12_otu_low<-data.frame(C12_otu_low, as.vector(Factors))
  
  #Graph (Highly Likely OTUs as boxplots)
  melted_low<-melt(C12_otu_low)
  melted_high<-melt(C12_otu_high)
  
  writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as C12-only.")
  print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
  writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be C12-only? There are ", nrow(highlikely),"of these OTUs?"))
  switch(menu(c("Visualize and edit","Do it myself later")), pain<-("Yes"), pain<-("No"))
  
  if (pain == "Yes") {
    writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
    writeLines("Your factor list is :")
    print(colnames(sample_data(p)))
    factor_2<-scan(n=1, what = character())
    
    highlikely_remove<-data.frame(highlikely, keep=NA)
    colnames(highlikely_remove)[1]<-"OTU"
    writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
    
    for (T in 1:nrow(highlikely)) {
      
      graphics.off()
      print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
      
      writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
      print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
      writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
      switch(menu(c("Keep","Discard")), keep<-("Yes"), keep<-("No"))
      
      if (keep == "Yes") {
        
        highlikely_remove$keep[T]<- "1"
        
      } else {
        
        highlikely_remove$keep[T]<- "0"
        
      }
    }
  }
  
  
  ##Keep the keepers 
  manual<-highlikely_remove
  manual<-cbind(manual, Likelihood_C12=c(rep("high", nrow(highlikely_remove))))
  write.csv(manual, file=paste("Output\\C12 OTUs - List of Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
  C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep == 1]))
  
  #Prune phyloseq to contain enriched OTUs
  p_C12<-prune_taxa(C12_otu, p)
  
  print(p_C12)
  
  saveRDS(p_C12, file=paste("Output\\BackUps\\p_All_c12_OTU_", sign_save, "_", cutoff,".rds", sep=""))
  write.csv(taxa_names(p_C12), file= paste("Output\\p_All_C12_OTU_", sign_save, "_", cutoff,".csv", sep=""))
  
}
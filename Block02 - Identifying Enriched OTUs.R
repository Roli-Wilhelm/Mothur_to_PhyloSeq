library(DESeq2)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(limma)
library(edgeR)

writeLines("Are you re-running the script starting at Block02?")

switch(menu(c("Yes","No")), rerun<-("Yes"), rerun<-("No"))

if (rerun == "Yes") {
  writeLines(paste("\nPlease ENTER your factor list (i.e. design matrix)"))
  Factors<-read.csv(file.choose(), header = T)
  Factors <- Factors[,colSums(is.na(Factors))<nrow(Factors)] #Remove accidental NA-filled columns
  
  #Input saved phyloseq object
  writeLines(paste("\nPlease ENTER your phyloseq object (It is likely located in /Output/BackUps/ and is your original imported data, unadultered)"))
  p<-readRDS(file.choose())
  
  #Transpose sample names to rownames
  Factors <- data.frame(Factors, row.names = sample_names(p), stringsAsFactors=F)
  Factors<-Factors[,-1]
  
  writeLines("\nPlease enter the cut-off you have been using (i.e. 0.01)")
  cutoff<-scan(n=1, what = character())
  writeLines("\nPlease enter the distinguishing name you originally used (i.e. \"Bacteria\")")
  sign_save<-scan(n=1, what = character())
  
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
writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C13 group.")
total<-nrow(otu_table(p))
number<-length(which(rowSums(otu_table(p_13C)) == 0))
percent1<-format(round(((number/total)*100), 2), nsmall = 2)
p_13C<-otu_table(p_13C)[-which(rowSums(otu_table(p_13C)) == 0),]

writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))

##Filter based on abundance, maybe rowSum > 5
writeLines("\nNext, we will filter by a minimum abundance value across all C13 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C13 samples.\n\n Please select a threshold number. [Recommended: 5] ")
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

#####
#Fork between Limma and DESeq and Relative Abundance Calculations
#####
writeLines("\nWould you like to use DESeq, Limma-voom or Relative Abundance to 12C Controls for differential expression profiling? [To do each one, you must re-run the script]")
switch(menu(c("DESeq", "Limma-voom", "Relative Abundance")), Method<-("DESeq"), Method<-("Limma"), Method<-("Relative Abundance"))

if (Method == "DESeq") {

#############################
###DESeq Differential Abundance Profiling
#############################
### Currently have to normalize counts in DESeq without the phyloseq function phyloseq_to_deseq2
### This is likely to change shortly once they release this functionality
#############################

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
writeLines("\nDifferential expression profiling, (package \"DESeq\") will be performed to identify which OTU are differentially expressed according to the enrichment status you have provided.\n Because the DESeq function can not handle calling functions as arguments in R, you will have to manually ENTER the column name in your factor list attributing enrichment status [NOTE: Please SAVE the script before answering following question].\n\n Further, you are able to modify the linear model [not recommended]. Have you added additional terms to linear model?")
file.edit("DESeq Command.R")
switch(menu(c("Yes", "No")), proceed<-("Yes"), proceed<-("No"))

if (proceed == "Yes") {
  
  writeLines("\nIf you entered a second factor in addition to \"enrichment\" please ENTER it now: [This is simply used for naming, so if multiple factors were used, enter the first and make a note of it.]")
  DESeq_factor<-scan(n=1, what = character())
  
  count<-counts
  cond<-conds
  
  ##Run DESeq
  source("DESeq Command.R")
  
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
  
  if (length(enr_otu) != 0) {
    
    ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
    ##Loop to parse the location of each designated OTU and pull that info into a new matrix
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
    
    #Prepare for Graph 
    melted_low<-melt(enr_otu_low)
    melted_high<-melt(enr_otu_high)
    
    writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
    print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
    
    #### Set-up the auto-remove loop
    writeLines("\nYou can now visually verify all of the OTUs which are borderline enriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
    writeLines(paste("\nWould you like to manually examine your", length(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
    switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
    
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
        test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
        
        if (logb(max(test$value), base=3) != 0) { 
          print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
          
          writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
          print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
          writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
          switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
          
          if (keep == "High") {
            
            lowlikely_remove$keep[T]<- "2"
            
          } 
          if (keep == "Yes") {
            
            lowlikely_remove$keep[T]<- "1"
            
          }
          if (keep == "No") {
            
            lowlikely_remove$keep[T]<- "0"
            
          }
        } else {
          
          writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
          
        }
      }
    }
    
    writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
    switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
    
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
        test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
        
        if (logb(max(test$value), base=3) != 0) { 
          print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
          
          writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
          print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
          writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
          switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
          
          if (keep == "High") {
            
            highlikely_remove$keep[T]<- "2"
            
          } 
          
          if (keep == "Yes") {
            
            highlikely_remove$keep[T]<- "1"
            
          }
          
          if (keep == "No") {
            
            highlikely_remove$keep[T]<- "0"
            
          }
        } else {
          
          writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
        }
      }
    }
    
    if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE) {
      ##Keep the keepers alongside the highlikely OTUs
      manual<-rbind(highlikely_remove, lowlikely_remove)
      manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
      write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
      enr_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]), as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1])))
      
      #Prune phyloseq to contain enriched OTUs
      p_enr_DESeq<-prune_taxa(enr_otu, p)
      
      saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
      
    } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
      ##Keep the keepers alongside the highlikely OTUs
      manual<-lowlikely_remove
      manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
      write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
      enr_otu<-as.vector(as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1]))
      
      #Prune phyloseq to contain enriched OTUs
      p_enr_DESeq<-prune_taxa(enr_otu, p)
      
      saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
      
    } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
      ##Keep the keepers alongside the highlikely OTUs
      manual<-highlikely_remove
      manual<-cbind(manual, Likelihood_Enriched=rep("high", nrow(highlikely_remove)))
      write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
      enr_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1])))
      
      #Prune phyloseq to contain enriched OTUs
      p_enr_DESeq<-prune_taxa(enr_otu, p)
      
      saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
      
    }
        
  } else {
    
    writeLines("\n\nWARNING : There were no significant OTUs attributed to your C13 Enrichment status. Try a different differential abundance profiling method.\n\n")
    
  }
    
  #####
  #C12 OTUs
  #####
  
  writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
  switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))
  
  if (proceed == "Yes") {
      p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)
      
      ###Pre-filtering steps
      #Determine which OTU are not found in the OTU table (after C13 sample have been removed)
      writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C12 group.")
      total<-nrow(otu_table(p))
      number<-length(which(rowSums(otu_table(p_12C)) == 0))
      percent1<-format(round(((number/total)*100), 2), nsmall = 2)
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) == 0),]
      
      writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))
      
      ##Filter based on abundance, maybe rowSum > 5
      writeLines("\nNext, we will filter by a minimum abundance value across all C12 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C12 samples.\n\n Please select a threshold number. [Recommended: 5] ")
      threshold<-scan(n=1, what = numeric())
      number<-length(which(rowSums(otu_table(p_12C)) < threshold))
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) < threshold),]
      percent2<-format(round(((number/total)*100), 2), nsmall = 2)
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))
      
      ##Filter based on occurrence in minimum number of samples > 3
      writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C12 samples.\n\n Please select a minimum number. [Recommended: 3]")
      min_number<-scan(n=1, what = numeric())
      present_absent<-otu_table(p_12C)
      present_absent[present_absent > 0] <- 1
      number<-length(which(rowSums(present_absent) < min_number))
      p_12C<-otu_table(p_12C)[-which(rowSums(present_absent) < min_number),]
      percent3<-format(round(((number/total)*100), 2), nsmall = 2)
      remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C12 dataset."))
      
      #Return to original phyloseq object
      p_12C<-intersect(rownames(p_12C), taxa_names(p_12C))
      p_12C<-prune_taxa(p_12C, p)
      
      writeLines("\nYour data table now looks like this:\n")
      print(p_12C)
      
      #This step is only made necessary by the fact I have to independently use DESeq, which conflicts with the formatting of a phyloSeq object
      rows<-nrow(otu_table(p_12C))
      cols<-ncol(otu_table(p_12C))
      counts<-as.matrix(otu_table(p_12C)[,])
      attributes(counts)<-NULL  #Remove phyloseq parameters
      counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
      rownames(counts)<-taxa_names(p_12C)   #Return row and column names
      colnames(counts)<-sample_names(p_12C)
      counts<-data.frame(counts)
      
      #Prepare the condition matrix 
      conds<-matrix(nrow=cols, ncol=1) #Rows are equal to columns here b/c were talking about samples not OTUs
      conds[which(sample_data(p_12C)[,status] == paste(label))] <- paste(label)
      conds[which(sample_data(p_12C)[,status] != paste(label))] <- "C12"
      conds<-data.frame(conds, as.vector(Factors[, which(colnames(Factors) != paste(status))]))
      rownames(conds) <- colnames(counts)
      
      #Rename column names
      colnames(conds)[1]<- paste(status)
      
      #DESeq calculates the geometric mean, which means all the zeros in my columns are problematic
      counts<-counts+1
      
      ##Run DESeq
      source("DESeq Command.R")
      
      ###Run Test for Differential Expression
      dds <- DESeq(DE)
      DE_results<-results(dds)
      DE_results<-DE_results[order(DE_results$log2FoldChange),]
      
      writeLines("\nDESeq will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C12.")
      
      writeLines("\nHere is the count data for your differentially abundant OTU.")
      print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
      
      writeLines("\nIs it showing abundance in your C12 samples?")
      
      switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
      
      if (order == "No") {
        
        DE_results$log2FoldChange<-DE_results$log2FoldChange*(-1)
        DE_results<-DE_results[order(DE_results$log2FoldChange),]
        writeLines("\nNow it should!")
        print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
        
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
      
      if (length(C12_otu) != 0) {
        
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
        switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
        
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
            test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
            
            if (logb(max(test$value), base=3) != 0) { 
              print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
              
              writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
              print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
              writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
              switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
              
              if (keep == "High") {
                
                highlikely_remove$keep[T]<- "2"
                
              } 
              if (keep == "Yes") {
                
                highlikely_remove$keep[T]<- "1"
                
              }
              if (keep == "No") {
                
                highlikely_remove$keep[T]<- "0"
                
              }
            } else {
              
              writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
            }
          }
        }
        
        #DO SAME FOR MARGINAL
        writeLines(paste("\nWould you like to manually examine the OTUs which were deemed marginally C12-only? There are ", length(lowlikely),"of these OTUs?"))
        switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
        
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
            test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
            
            if (logb(max(test$value), base=3) != 0) { 
              print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
              
              writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
              print(tax_table(p)[taxa_names(p) == lowlikely_remove$OTU[T],])
              writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
              switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
              
              if (keep == "High") {
                
                lowlikely_remove$keep[T]<- "2"
                
              } 
              if (keep == "Yes") {
                
                lowlikely_remove$keep[T]<- "1"
                
              }
              if (keep == "No") {
                
                lowlikely_remove$keep[T]<- "0"
                
              }
            } else {
              
              writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
            }
          }
        }
        
        if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE) {
          
          ##Keep the keepers 
          manual<-rbind(highlikely_remove, lowlikely_remove)
          manual<-cbind(manual, Likelihood_C12=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
          write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
          C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1] & lowlikely[lowlikely_remove$keep >= 1]))
          
          #Prune phyloseq to contain enriched OTUs
          p_C12_DESeq<-prune_taxa(C12_otu, p)
          
          print(p_C12_DESeq)
          
          saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
          write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
                    
        } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
          ##Keep the keepers 
          manual<-lowlikely_remove
          manual<-cbind(manual, Likelihood_C12=rep("low", nrow(lowlikely_remove)))
          write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
          C12_otu<-as.vector(as.matrix(lowlikely[lowlikely_remove$keep >= 1]))
          
          #Prune phyloseq to contain enriched OTUs
          p_C12_DESeq<-prune_taxa(C12_otu, p)
          
          print(p_C12_DESeq)
          
          saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
          write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
                    
        } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
          ##Keep the keepers 
          manual<-highlikely_remove
          manual<-cbind(manual, Likelihood_C12=rep("high", nrow(highlikely_remove)))
          write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", DESeq_factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
          C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]))
          
          #Prune phyloseq to contain enriched OTUs
          p_C12_DESeq<-prune_taxa(C12_otu, p)
          
          print(p_C12_DESeq)
          
          saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
          write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", DESeq_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
                    
        }
                
      } else {
        
        writeLines("\n\nWARNING : There were no significant OTUs attributed to your C12 unenriched status. Try a different differential abundance profiling method if you are expecting for results.\n\n")
        
      }
    }
  }

if (proceed == "No") {
  
  #Prompt for specific factor to subset according to separate factors
  writeLines("\nWould you like to specific factor to subset your data by and apply DESeq separately to? If so, please enter it now. If you do not wish to, just hit ENTER [recommended for first pass].")
  writeLines("Your factor list is :")
  print(colnames(sample_data(p)))
  DESeq_divide<-scan(n=1, what = character())
  
  if (length(DESeq_divide) != 0) {
    
    writeLines("\nAre you returning to re-analyze the C12 (unenriched) community? If you would you like to skip the 13C OTU selection step, do so now.")
    switch(menu(c("Do 13C OTU selection", "Skip to 12C")), skip<-("No"), skip<-("Yes"))
    
    if (skip == "No") {
    
      divide_list<-grep(DESeq_divide, colnames(conds))
      unique_factor<-as.vector(conds[which(duplicated(conds[,divide_list]) == FALSE),divide_list])
      
      for (go in 1:length(unique_factor)) {
        
        ##Prepare Count Data
        #Append factor data to count data
        count<-t(counts)
        count<-cbind(count, conds)
        
        #Split count data into single factor
        divide<-unique_factor[go]
        count<-subset(count, get(DESeq_divide) == paste(divide))
        
        #Transpose back into OTU x Sample table
        count<-data.frame(count[,1:(ncol(count)-ncol(conds))])
        count<-data.frame(t(count))
        
        ##Prepare Factor Data
        cond<-subset(conds, get(DESeq_divide) == paste(divide))
        
        ##Run DESeq
        source("DESeq Command.R")
        
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
        
        if (length(enr_otu) != 0) {
          
          ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
          
          ##Loop to parse the location of each designated OTU and pull that info into a new matrix
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
          
          #Prepare for Graph 
          melted_low<-melt(enr_otu_low)
          melted_high<-melt(enr_otu_high)
          
          writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
          print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
          
          #### Set-up the auto-remove loop
          writeLines("\nYou can now visually verify all of the OTUs which are borderline enriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
          writeLines(paste("\nWould you like to manually examine your", length(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  lowlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  lowlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  lowlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
                           
              }
            }
          }
          
          writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
              
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  highlikely_remove$keep[T]<- "2"
                  
                } 
                
                if (keep == "Yes") {
                  
                  highlikely_remove$keep[T]<- "1"
                  
                }
                
                if (keep == "No") {
                  
                  highlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
          
          if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE){
            ##Keep the keepers alongside the highlikely OTUs
            manual<-rbind(highlikely_remove, lowlikely_remove)
            manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
            manual$Subsetted_by<- paste(divide)
          
          } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
            ##Keep the keepers alongside the highlikely OTUs
            manual<-highlikely_remove
            manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
            manual$Subsetted_by<- paste(divide)
            
          } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
            ##Keep the keepers alongside the highlikely OTUs
            manual<-highlikely_remove
            manual<-cbind(manual, Likelihood_Enriched=rep("high", nrow(highlikely_remove)))
            manual$Subsetted_by<- paste(divide)
            
          }
                    
          assign(divide, manual)
          
        } else { 
          
          manual <- data.frame(OTU="None Selected", keep = NA, Likelihood_Enriched = NA, Subsetted_by = paste(divide))
          assign(divide, manual)
        
        }
      
      }
         
      manual<-do.call(rbind, lapply(unique_factor, as.symbol))
      
      if (any(duplicated(manual$OTU) == TRUE)) {
        
        manual<-manual[-which(duplicated(manual$OTU)),]
        
      }
           
      ##Keep the keepers alongside the highlikely OTUs
      write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", DESeq_divide, "_", sign_save, "_", cutoff, ".csv", sep=""))
      enr_otu<-as.vector(manual[which(manual$keep >= 1),1])
      
      #Prune phyloseq to contain enriched OTUs
      p_enr_DESeq<-prune_taxa(enr_otu, p)
      
      saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", DESeq_divide, "_",sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", DESeq_divide, "_",sign_save, "_", cutoff,".csv", sep=""))
      
    }
    
    ####
    ##Do C12 OTUs
    ####
        
    writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
    switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))
    
    if (proceed == "Yes") {
      
      #Carrying forward original info in case it has been modified
      p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)
      divide_list<-grep(DESeq_divide, colnames(conds))
      unique_factor<-as.vector(conds[which(duplicated(conds[,divide_list]) == FALSE),divide_list])
      
      ###Pre-filtering steps
      #Determine which OTU are not found in the OTU table (after C13 sample have been removed)
      writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C12 group.")
      total<-nrow(otu_table(p))
      number<-length(which(rowSums(otu_table(p_12C)) == 0))
      percent1<-format(round(((number/total)*100), 2), nsmall = 2)
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) == 0),]
      
      writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))
      
      ##Filter based on abundance, maybe rowSum > 5
      writeLines("\nNext, we will filter by a minimum abundance value across all C12 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C12 samples.\n\n Please select a threshold number. [Recommended: 5] ")
      threshold<-scan(n=1, what = numeric())
      number<-length(which(rowSums(otu_table(p_12C)) < threshold))
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) < threshold),]
      percent2<-format(round(((number/total)*100), 2), nsmall = 2)
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))
            
      ##Filter based on occurrence in minimum number of samples > 3
      writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C12 samples.\n\n Please select a minimum number. [Recommended: 3]")
      min_number<-scan(n=1, what = numeric())
      present_absent<-otu_table(p_12C)
      present_absent[present_absent > 0] <- 1
      number<-length(which(rowSums(present_absent) < min_number))
      p_12C<-otu_table(p_12C)[-which(rowSums(present_absent) < min_number),]
      percent3<-format(round(((number/total)*100), 2), nsmall = 2)
      remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C12 dataset."))
      
      #Return to original phyloseq object
      p_12C<-intersect(rownames(p_12C), taxa_names(p_12C))
      p_12C<-prune_taxa(p_12C, p)
      
      writeLines("\nYour data table now looks like this:\n")
      print(p_12C)
      
      #This step is only made necessary by the fact I have to independently use DESeq, which conflicts with the formatting of a phyloSeq object
      rows<-nrow(otu_table(p_12C))
      cols<-ncol(otu_table(p_12C))
      counts<-as.matrix(otu_table(p_12C)[,])
      attributes(counts)<-NULL  #Remove phyloseq parameters
      counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
      rownames(counts)<-taxa_names(p_12C)   #Return row and column names
      colnames(counts)<-sample_names(p_12C)
      counts<-data.frame(counts)
      
      #Prepare the condition matrix 
      conds<-matrix(nrow=cols, ncol=1) #Rows are equal to columns here b/c were talking about samples not OTUs
      conds[which(sample_data(p_12C)[,status] == paste(label))] <- paste(label)
      conds[which(sample_data(p_12C)[,status] != paste(label))] <- "C12"
      conds<-data.frame(conds, as.vector(Factors[, which(colnames(Factors) != paste(status))]))
      rownames(conds) <- colnames(counts)
      
      #Rename column names
      colnames(conds)[1]<- paste(status)
      
      #DESeq calculates the geometric mean, which means all the zeros in my columns are problematic
      counts<-counts+1
      
      for (go in 1:length(unique_factor)) {
        
        ##Prepare Count Data
        #Append factor data to count data
        count<-t(count)
        count<-cbind(count, conds)
        
        #Split count data into single factor
        divide<-unique_factor[go]
        count<-subset(count, get(DESeq_divide) == paste(divide))
        
        #Transpose back into OTU x Sample table
        count<-data.frame(counts_foo[,1:(ncol(count)-ncol(conds))])
        count<-data.frame(t(count))
        
        ##Prepare Factor Data
        cond<-subset(conds, get(DESeq_divide) == paste(divide))
        
        ##Run DESeq
        source("DESeq Command.R")
        
        ###Run Test for Differential Expression
        dds <- DESeq(DE)
        DE_results<-results(dds)
        DE_results<-DE_results[order(DE_results$log2FoldChange),]
      
        writeLines("\nDESeq will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C12.")
      
        writeLines("\nHere is the count data for your differentially abundant OTU.")
        print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
      
        writeLines("\nIs it showing abundance in your C12 samples?")
      
        switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
      
        if (order == "No") {
        
          DE_results$log2FoldChange<-DE_results$log2FoldChange*(-1)
          DE_results<-DE_results[order(DE_results$log2FoldChange),]
          writeLines("\nNow it should!")
          print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
        
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
      
        if (length(C12_otu) !=0) {
          
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
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
              
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  highlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  highlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  highlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
          
          #DO SAME FOR MARGINAL
          writeLines(paste("\nWould you like to manually examine the OTUs which were deemed marginally C12-only? There are ", length(lowlikely),"of these OTUs?"))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p)[taxa_names(p) == lowlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  lowlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  lowlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  lowlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
          
          if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE){
            ##Make a dataframe containing the likelihood status and the "keep" status
            manual<-rbind(highlikely_remove, lowlikely_remove)
            manual<-cbind(manual, Likelihood_C12=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
            manual$Subsetted_by<- paste(divide)
            
          } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
            ##Make a dataframe containing the likelihood status and the "keep" status
            manual<-lowlikely_remove
            manual<-cbind(manual, Likelihood_C12=rep("low", nrow(lowlikely_remove)))
            manual$Subsetted_by<- paste(divide)
            
          } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
            ##Make a dataframe containing the likelihood status and the "keep" status
            manual<-highlikely_remove
            manual<-cbind(manual, Likelihood_C12=rep("high", nrow(highlikely_remove)))
            manual$Subsetted_by<- paste(divide)
            
          }
          
              
          assign(divide, manual)
          
        } else {
    
          manual <- data.frame(OTU="None Selected", keep = NA, Likelihood_C12 = NA, Subsetted_by = paste(divide))
          assign(divide, manual)
          
        }
      }
      
      manual<-do.call(rbind, lapply(unique_factor, as.symbol))
      
      if (any(duplicated(manual$OTU)) == TRUE) {
      
        manual<-manual[-which(duplicated(manual$OTU)),]
        
      }
      
      ##Keep the keepers 
      write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", DESeq_divide, "_", sign_save, "_", cutoff, ".csv", sep=""))
      C12_otu<-as.vector(manual[which(manual$keep >= 1),1]) 
      
      #Prune phyloseq to contain enriched OTUs
      p_C12_DESeq<-prune_taxa(C12_otu, p)
      
      print(p_C12_DESeq)
      
      saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", DESeq_divide, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", DESeq_divide, "_", sign_save, "_", cutoff,".csv", sep=""))
          
      }
    }

  if (length(DESeq_divide) == 0) {
    
    count<-counts
    cond<-conds
    
    ##Run DESeq
    source("DESeq Command.R")
    
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
    
    if (length(enr_otu) != 0) {
      ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
      ##Loop to parse the location of each designated OTU and pull that info into a new matrix
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
      
      #Prepare for Graph 
      melted_low<-melt(enr_otu_low)
      melted_high<-melt(enr_otu_high)
      
      writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
      print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
      
      #### Set-up the auto-remove loop
      writeLines("\nYou can now visually verify all of the OTUs which are borderline enriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
      writeLines(paste("\nWould you like to manually examine your", length(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
      switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
      
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
          test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
          
          if (logb(max(test$value), base=3) != 0) { 
            print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
          
            writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
            print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
            writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
            switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
            
            if (keep == "High") {
              
              lowlikely_remove$keep[T]<- "2"
              
            } 
            if (keep == "Yes") {
              
              lowlikely_remove$keep[T]<- "1"
              
            }
            if (keep == "No") {
              
              lowlikely_remove$keep[T]<- "0"
              
            }
          } else {
            writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
            
          }
        }
      }
      
      writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
      switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
      
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
          test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
          
          if (logb(max(test$value), base=3) != 0) { 
            print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
            
            writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
            print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
            writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
            switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
            
            if (keep == "High") {
              
              highlikely_remove$keep[T]<- "2"
              
            } 
            
            if (keep == "Yes") {
              
              highlikely_remove$keep[T]<- "1"
              
            }
            
            if (keep == "No") {
              
              highlikely_remove$keep[T]<- "0"
              
            }
          } else {
            
            writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
          }
        }
      }
      
      if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE){
        ##Keep the keepers alongside the highlikely OTUs
        manual<-rbind(highlikely_remove, lowlikely_remove)
        manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
        write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        enr_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]), as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1])))
        
        #Prune phyloseq to contain enriched OTUs
        p_enr_DESeq<-prune_taxa(enr_otu, p)
        
        saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
        write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
      
      } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
        ##Keep the keepers alongside the highlikely OTUs
        manual<-lowlikely_remove
        manual<-cbind(manual, rep("low", nrow(lowlikely_remove)))
        write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        enr_otu<-as.vector(as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1]))
        
        #Prune phyloseq to contain enriched OTUs
        p_enr_DESeq<-prune_taxa(enr_otu, p)
        
        saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
        write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
                
      } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
        ##Keep the keepers alongside the highlikely OTUs
        manual<-highlikely_remove
        manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove))))
        write.csv(manual, file=paste("Output\\C13 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        enr_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]))
        
        #Prune phyloseq to contain enriched OTUs
        p_enr_DESeq<-prune_taxa(enr_otu, p)
        
        saveRDS(p_enr_DESeq, file=paste("Output\\BackUps\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
        write.csv(taxa_names(p_enr_DESeq), file= paste("Output\\p_Enriched_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
                
      }
            
    } else {
      
      writeLines("\n\nWARNING : There were no significant OTUs attributed to your C13 Enrichment status. Try a different differential abundance profiling method.\n\n")
      
    }
    
    writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
    switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))
    
    if (proceed == "Yes") {
      p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)
      
      ###Pre-filtering steps
      #Determine which OTU are not found in the OTU table (after C13 sample have been removed)
      writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C12 group.")
      total<-nrow(otu_table(p))
      number<-length(which(rowSums(otu_table(p_12C)) == 0))
      percent1<-format(round(((number/total)*100), 2), nsmall = 2)
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) == 0),]
      
      writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))
      
      ##Filter based on abundance, maybe rowSum > 5
      writeLines("\nNext, we will filter by a minimum abundance value across all C12 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C12 samples.\n\n Please select a threshold number. [Recommended: 5] ")
      threshold<-scan(n=1, what = numeric())
      number<-length(which(rowSums(otu_table(p_12C)) < threshold))
      p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) < threshold),]
      percent2<-format(round(((number/total)*100), 2), nsmall = 2)
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))
      
      ##Filter based on occurrence in minimum number of samples > 3
      writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C12 samples.\n\n Please select a minimum number. [Recommended: 3]")
      min_number<-scan(n=1, what = numeric())
      present_absent<-otu_table(p_12C)
      present_absent[present_absent > 0] <- 1
      number<-length(which(rowSums(present_absent) < min_number))
      p_12C<-otu_table(p_12C)[-which(rowSums(present_absent) < min_number),]
      percent3<-format(round(((number/total)*100), 2), nsmall = 2)
      remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
      writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C12 dataset."))
      
      #Return to original phyloseq object
      p_12C<-intersect(rownames(p_12C), taxa_names(p_12C))
      p_12C<-prune_taxa(p_12C, p)
      
      writeLines("\nYour data table now looks like this:\n")
      print(p_12C)
      
      #This step is only made necessary by the fact I have to independently use DESeq, which conflicts with the formatting of a phyloSeq object
      rows<-nrow(otu_table(p_12C))
      cols<-ncol(otu_table(p_12C))
      counts<-as.matrix(otu_table(p_12C)[,])
      attributes(counts)<-NULL  #Remove phyloseq parameters
      counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
      rownames(counts)<-taxa_names(p_12C)   #Return row and column names
      colnames(counts)<-sample_names(p_12C)
      counts<-data.frame(counts)
      
      #Prepare the condition matrix 
      conds<-matrix(nrow=cols, ncol=1) #Rows are equal to columns here b/c were talking about samples not OTUs
      conds[which(sample_data(p_12C)[,status] == paste(label))] <- paste(label)
      conds[which(sample_data(p_12C)[,status] != paste(label))] <- "C12"
      conds<-data.frame(conds, as.vector(Factors[, which(colnames(Factors) != paste(status))]))
      rownames(conds) <- colnames(counts)
      
      #Rename column names
      colnames(conds)[1]<- paste(status)
      
      #DESeq calculates the geometric mean, which means all the zeros in my columns are problematic
      counts<-counts+1
      
      ##Run DESeq
      source("DESeq Command.R")
      
      ###Run Test for Differential Expression
      dds <- DESeq(DE)
      DE_results<-results(dds)
      DE_results<-DE_results[order(DE_results$log2FoldChange),]
      
      writeLines("\nDESeq will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C12.")
      
      writeLines("\nHere is the count data for your differentially abundant OTU.")
      print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
      
      writeLines("\nIs it showing abundance in your C12 samples?")
      
      switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
      
      if (order == "No") {
        
        DE_results$log2FoldChange<-DE_results$log2FoldChange*(-1)
        DE_results<-DE_results[order(DE_results$log2FoldChange),]
        writeLines("\nNow it should!")
        print(otu_table(p)[rownames(otu_table(p)) == rownames(DE_results)[1]])
        
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
      
        if (length(C12_otu) != 0) {
          
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
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  highlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  highlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  highlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
                
              }
            }
          }
          
          #DO SAME FOR MARGINAL
          writeLines(paste("\nWould you like to manually examine the OTUs which were deemed marginally C12-only? There are ", length(lowlikely),"of these OTUs?"))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p)[taxa_names(p) == lowlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(lowlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  lowlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  lowlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  lowlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
          
          if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE){
            
            ##Keep the keepers 
            manual<-rbind(highlikely_remove, lowlikely_remove)
            manual<-cbind(manual, Likelihood_C12=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
            
            write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
            C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1] & lowlikely[lowlikely_remove$keep >= 1]))
            
            #Prune phyloseq to contain enriched OTUs
            p_C12_DESeq<-prune_taxa(C12_otu, p)
            
            print(p_C12_DESeq)
            
            saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
            write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
            
          } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
            
            ##Keep the keepers 
            manual<-lowlikely_remove
            manual<-cbind(manual, Likelihood_C12=rep("low", nrow(lowlikely_remove)))
            
            write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
            C12_otu<-as.vector(as.matrix(lowlikely[lowlikely_remove$keep >= 1]))
            
            #Prune phyloseq to contain enriched OTUs
            p_C12_DESeq<-prune_taxa(C12_otu, p)
            
            print(p_C12_DESeq)
            
            saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
            write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
            
          } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
            
            ##Keep the keepers 
            manual<-rbind(highlikely_remove)
            manual<-cbind(manual, Likelihood_C12=rep("high", nrow(highlikely_remove)))
            
            write.csv(manual, file=paste("Output\\C12 OTUs - DESeq - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
            C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1] & lowlikely[lowlikely_remove$keep >= 1]))
            
            #Prune phyloseq to contain enriched OTUs
            p_C12_DESeq<-prune_taxa(C12_otu, p)
            
            print(p_C12_DESeq)
            
            saveRDS(p_C12_DESeq, file=paste("Output\\BackUps\\p_c12_OTU_DESeq_", sign_save, "_", cutoff,".rds", sep=""))
            write.csv(taxa_names(p_C12_DESeq), file= paste("Output\\p_C12_OTU_DESeq_", sign_save, "_", cutoff,".csv", sep=""))
          }
          
          
        } else { 
        
          writeLines("\n\nWARNING : There were no significant OTUs attributed to your C12 unenriched status. Try a different differential abundance profiling method if you are expecting for results.\n\n")
                    
        }
      }
    }
  
  }
}

if (Method == "Limma") {

#############################
###Limma + Voom Differential Abundance Profiling
#############################
### Voom is newly released (2014) and supposedly equally as good at limiting the # of false-positive as DESeq, but less 
### conservative and more able to identify differentially expressed genes.
#############################

#Get count data in correct format
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

#Prompt for specific factor to subset according to separate factors
writeLines("\nWould you like to specific factor to subset your data by and apply Limma-Voom separately to? If so, please enter it now. If you do not wish to, simply press ENTER [recommended for first pass].")
writeLines("Your factor list is :")
print(colnames(sample_data(p)))
VOOM_divide<-scan(n=1, what = character())

if (length(VOOM_divide) != 0) {
  
  writeLines("\nAre you returning to re-analyze the C12 (unenriched) community? If you would you like to skip the 13C OTU selection step, do so now.")
  switch(menu(c("Do 13C OTU selection", "Skip to 12C")), skip<-("No"), skip<-("Yes"))
  
  if (skip == "No") {
    
    divide_list<-grep(VOOM_divide, colnames(conds))
    unique_factor<-as.vector(conds[which(duplicated(conds[,divide_list]) == FALSE),divide_list])
    
    for (go in 1:length(unique_factor)) {
      ##Prepare Count Data
      #Append factor data to count data
      counts_foo<-t(counts)
      counts_foo<-cbind(counts_foo, conds)
      
      #Split count data into single factor
      divide<-unique_factor[go]
      counts_foo<-subset(counts_foo, get(VOOM_divide) == paste(divide))
      
      #Transpose back into OTU x Sample table
      counts_foo<-data.frame(counts_foo[,1:(ncol(counts_foo)-ncol(conds))])
      counts_foo<-data.frame(t(counts_foo))
      
      ##Prepare Factor Data
      conds_foo<-subset(conds, get(VOOM_divide) == paste(divide))
      
      #Make design matrix
      Enrichment <-conds_foo[,colnames(conds_foo) == status]
      design <- model.matrix(~Enrichment)
      
      #the voom-limma process
      norm.factor <- calcNormFactors(counts_foo)
      dat.voomed <- voom(counts_foo, design, plot = TRUE, lib.size = colSums(counts_foo) * norm.factor)
      fit <- lmFit(dat.voomed, design)
      fit <- eBayes(fit)
      
      TopHits<-topTable(fit, coef = paste("Enrichment", label, sep=""), sort.by = "logFC", adjust.method="BH", number=Inf)
      TopHits<-TopHits[with(TopHits, order(-logFC)), ]
      
      writeLines("\nLimma-voom will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C13.")
      writeLines("\nHere is the count data for your differentially abundant OTU.")
      print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
      
      writeLines("\nIs it showing abundance in your C13 samples?")
      switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
      
      if (order == "No") {
        TopHits$logFC<-TopHits$logFC*(-1)
        TopHits<-TopHits[with(TopHits, order(-logFC)), ]
        writeLines("\nNow it should!")
        print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
        
      }
      
      if (any(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0) != FALSE) {
        
        #Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
        enr_otu<-TopHits[which(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0),]
        highlikely<-which(enr_otu$logFC > 1.5)
        lowlikely<-which(enr_otu$logFC < 1.5)
        
        highlikely<-rownames(enr_otu)[highlikely]
        highlikely<-data.frame(highlikely)
        colnames(highlikely)<-"OTU"
        lowlikely<-rownames(enr_otu)[lowlikely]
        lowlikely<-data.frame(lowlikely)
        colnames(lowlikely)<-"OTU"
        
        enr_otu<-rownames(enr_otu)
        
        writeLines(paste("\nIn total you have recovered ", nrow(enr_otu)," enriched OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be enriched and ", nrow(lowlikely)," will have to be manually inspected and confirmed."))
        
        ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
        d=enr_otu
        OTU_location<- vector(mode="numeric", length=length(d))
        
        for (mercy in 1:length(d)) {
          
          OTU_location[mercy]<-grep(enr_otu[mercy], rownames(counts))
          
        }
        
        enr_otu<-counts[OTU_location,]
        
        #Make into data.frame
        enr_otu<-data.frame(enr_otu)
        enr_otu$likelihood<-NA
        
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
        
        if (ncol(enr_otu_low) != 0) {
          
          #Add in factors
          enr_otu_low<-data.frame(enr_otu_low, as.vector(Factors))
          
          #Graph (Highly Likely OTUs as boxplots)
          melted_low<-melt(enr_otu_low)
          
          writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
          print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
          
          #### Set-up the auto-remove loop
          writeLines("\nYou can now visually verify all of the OTUs which are borderline enriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
          writeLines(paste("\nWould you like to manually examine your", nrow(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
          if (pain == "Yes") {
            writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
            writeLines("Your factor list is :")
            print(colnames(sample_data(p)))
            factor_2<-scan(n=1, what = character())
            
            lowlikely_remove<-data.frame(lowlikely, keep=NA)
            colnames(lowlikely_remove)[1]<-"OTU"
            writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
            for (T in 1:nrow(lowlikely)) {
              
              graphics.off()
              test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(lowlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  lowlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  lowlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  lowlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
          
        }
        
        if (ncol(enr_otu_high) != 0) {
         
          #Add in factors
          enr_otu_high<-data.frame(enr_otu_high, as.vector(Factors))
          
          #Visualize
          melted_high<-melt(enr_otu_high)
          
          writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
          print(ddply(melted_high, ~ Enrichment, summarize, value=sum(value)))
          
          writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
          switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
          
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
              test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
              
              if (logb(max(test$value), base=3) != 0) { 
                print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
                
                writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
                print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
                writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
                switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
                
                if (keep == "High") {
                  
                  highlikely_remove$keep[T]<- "2"
                  
                } 
                if (keep == "Yes") {
                  
                  highlikely_remove$keep[T]<- "1"
                  
                }
                if (keep == "No") {
                  
                  highlikely_remove$keep[T]<- "0"
                  
                }
              } else {
                
                writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
              }
            }
          }
        }
        
        if (ncol(enr_otu_high) != 0 & (ncol(enr_otu_low) != 0)) {
          
          ##Keep the keepers alongside the highlikely OTUs
          manual<-rbind(highlikely_remove, lowlikely_remove)
          manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
          manual$Subsetted_by<- paste(divide)
          
        }
        
        if (ncol(enr_otu_high) != 0 & (ncol(enr_otu_low) == 0)) {
          
          ##Keep the keepers alongside the highlikely OTUs
          manual<-rbind(highlikely_remove)
          manual<-cbind(manual, Likelihood_Enriched=rep("high", nrow(highlikely_remove)))
          manual$Subsetted_by<- paste(divide)
          
        }
        
        if (ncol(enr_otu_high) == 0 & (ncol(enr_otu_low) != 0)) {
          
          ##Keep the keepers alongside the highlikely OTUs
          manual<-rbind(lowlikely_remove)
          manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
          manual$Subsetted_by<- paste(divide)
          
        }
                
      } else { 
        
        manual <- data.frame(OTU=0,keep=0,Likelihood_Enriched=0,Subsetted_by=0)
        
      }
      
      assign(divide, manual)

    }
      
    manual<-do.call(rbind, lapply(unique_factor, as.symbol))
    
    if (any(duplicated(manual$OTU)) == TRUE) {
      manual<-manual[-which(duplicated(manual$OTU)),]
    }
    
    if (any(manual$OTU == 0) == TRUE) {
    manual<-manual[-which(manual$OTU == 0),]
    
    }
    
    ##Keep the keepers alongside the highlikely OTUs
    write.csv(manual, file=paste("Output\\C13 OTUs - voom - Manual Edits_", VOOM_divide, "_", sign_save, "_", cutoff, ".csv", sep=""))
    enr_otu<-as.vector(manual[which(manual$keep >= 1),1])
      
    #Prune phyloseq to contain enriched OTUs
    p_enr_voom<-prune_taxa(enr_otu, p)
      
    saveRDS(p_enr_voom, file=paste("Output\\BackUps\\p_Enriched_OTU_voom_", VOOM_divide, "_", sign_save, "_", cutoff,".rds", sep=""))
    write.csv(taxa_names(p_enr_voom), file= paste("Output\\p_Enriched_OTU_voom_", VOOM_divide, "_", sign_save, "_", cutoff,".csv", sep=""))
      
    }
    
  ####
  #Repeat for C12 OTUs
  ####
  
  writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
  switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))
  
  if (proceed == "Yes") {
    p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)
    divide_list<-grep(VOOM_divide, colnames(conds))
    unique_factor<-as.vector(conds[which(duplicated(conds[,divide_list]) == FALSE),divide_list])
    
    ###Pre-filtering steps
    #Determine which OTU are not found in the OTU table (after C12 sample have been removed)
    writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C12 group.")
    total<-nrow(otu_table(p))
    number<-length(which(rowSums(otu_table(p_12C)) == 0))
    percent1<-format(round(((number/total)*100), 2), nsmall = 2)
    p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) == 0),]
    
    writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))
    
    ##Filter based on abundance, maybe rowSum > 5
    writeLines("\nNext, we will filter by a minimum abundance value across all C12 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C12 samples.\n\n Please select a threshold number. [Recommended: 5] ")
    threshold<-scan(n=1, what = numeric())
    number<-length(which(rowSums(otu_table(p_12C)) < threshold))
    p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) < threshold),]
    percent2<-format(round(((number/total)*100), 2), nsmall = 2)
    writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))
    
    ##Filter based on occurrence in minimum number of samples > 3
    writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C12 samples.\n\n Please select a minimum number. [Recommended: 3]")
    min_number<-scan(n=1, what = numeric())
    present_absent<-otu_table(p_12C)
    present_absent[present_absent > 0] <- 1
    number<-length(which(rowSums(present_absent) < min_number))
    p_12C<-otu_table(p_12C)[-which(rowSums(present_absent) < min_number),]
    percent3<-format(round(((number/total)*100), 2), nsmall = 2)
    remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
    writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C12 dataset."))
    
    #Return to original phyloseq object
    p_12C<-intersect(rownames(p_12C), taxa_names(p_12C))
    p_12C<-prune_taxa(p_12C, p)
    
    writeLines("\nYour data table now looks like this:\n")
    print(p_12C)
    
    #Only necessary due to conflicts with the formatting of a phyloSeq object
    rows<-nrow(otu_table(p_12C))
    cols<-ncol(otu_table(p_12C))
    counts<-as.matrix(otu_table(p_12C)[,])
    attributes(counts)<-NULL  #Remove phyloseq parameters
    counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
    rownames(counts)<-taxa_names(p_12C)   #Return row and column names
    colnames(counts)<-sample_names(p_12C)
    counts<-data.frame(counts)
    
    for (go in 1:length(unique_factor)) {
      
      ##Prepare Count Data
      #Append factor data to count data
      counts_foo<-counts
      
      #Remove any columns and rows which are equal or less than one
      if (any(rowSums(counts_foo) == 0) == TRUE) {
        counts_foo<- counts_foo[-which(rowSums(counts_foo) == 0),]
      }
      
      if (any(colSums(counts_foo) <= 1) == TRUE) {
        factor_erase<- colnames(counts_foo)[which(colSums(counts_foo) <= 1)]
        counts_foo<- counts_foo[,-which(colSums(counts_foo) <= 1)]
      }
      
      counts_foo<-t(counts_foo)
      
      if (exists("factor_erase") == TRUE) {
        counts_foo<-cbind(counts_foo, conds[-which(rownames(conds) == factor_erase),])
      } else {
        counts_foo<-cbind(counts_foo, conds)
      }
      
      #Split count data into single factor
      divide<-unique_factor[go]
      counts_foo<-subset(counts_foo, get(VOOM_divide) == paste(divide))
      
      #Transpose back into OTU x Sample table
      counts_foo<-data.frame(counts_foo[,1:(ncol(counts_foo)-ncol(conds))])
      counts_foo<-data.frame(t(counts_foo))
      
      ##Prepare Factor Data
      if (exists("factor_erase") == TRUE) {
        conds_foo<-subset(conds[-which(rownames(conds) == factor_erase),], get(VOOM_divide) == paste(divide))
      } else {
        conds_foo<-subset(conds, get(VOOM_divide) == paste(divide))
      }
                  
      #Make design matrix
      Enrichment <-conds_foo[,colnames(conds_foo) == status]
      design <- model.matrix(~Enrichment)
      
      #the voom-limma process
      norm.factor <- calcNormFactors(counts_foo)
      dat.voomed <- voom(counts_foo, design, plot = TRUE, lib.size = colSums(counts_foo) * norm.factor)
      fit <- lmFit(dat.voomed, design)
      fit <- eBayes(fit)
      
      TopHits<-topTable(fit, coef = paste("Enrichment", label, sep=""), sort.by = "logFC", adjust.method="BH", number=Inf)
      TopHits<-TopHits[with(TopHits, order(-logFC)), ]
      
      writeLines("\nLimma-voom will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C13.")
      writeLines("\nHere is the count data for your differentially abundant OTU.")
      print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
      
      writeLines("\nIs it showing abundance in your C12 samples?")
      switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
      
      if (order == "No") {
        TopHits$logFC<-TopHits$logFC*(-1)
        TopHits<-TopHits[with(TopHits, order(-logFC)), ]
        writeLines("\nNow it should!")
        print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
        
      }
      
      if (any(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0) != FALSE) {
        
      #Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
      C12_otu<-TopHits[which(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0),]
      highlikely<-which(C12_otu$logFC > 1.5)
      lowlikely<-which(C12_otu$logFC < 1.5)
      
      highlikely<-rownames(C12_otu)[highlikely]
      highlikely<-data.frame(highlikely)
      colnames(highlikely)<-"OTU"
      lowlikely<-rownames(C12_otu)[lowlikely]
      lowlikely<-data.frame(lowlikely)
      colnames(lowlikely)<-"OTU"
      
      C12_otu<-rownames(C12_otu)
      
      writeLines(paste("\nIn total you have recovered ", nrow(C12_otu)," enriched OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be enriched and ", nrow(lowlikely)," is marginal and (unlike 13C selection) will be ignored."))
      
      ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
      d=C12_otu
      OTU_location<- vector(mode="numeric", length=length(d))
      
      for (mercy in 1:length(d)) {
        
        OTU_location[mercy]<-grep(C12_otu[mercy], rownames(counts))
        
      }
      
      C12_otu<-counts[OTU_location,]
      
      #Make into data.frame
      C12_otu<-data.frame(C12_otu)
      C12_otu$likelihood<-NA
      
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
    
      if (ncol(C12_otu_low) != 0) {
        
        #Add in factors
        C12_otu_low<-data.frame(C12_otu_low, as.vector(Factors))
        
        #Graph 
        melted_low<-melt(C12_otu_low)
        
        writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
        print(ddply(melted_low, ~ Enrichment, summarize, value=sum(value)))
        
        #### Set-up the auto-remove loop
        writeLines("\nYou can now visually verify all of the OTUs which are borderline C12 unenriched. Whatever method you employ here, ensure that it is consistent. This is a time-saveing step so that you do not have to employ a variety of logical filters to your data. Depending on the number of OTUs you have at this point, it may make sense to encode some filters.\n\nYou will be able to discard or keep each OTU on the spot.") 
        writeLines(paste("\nWould you like to manually examine your", nrow(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
        switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
        
        if (pain == "Yes") {
          writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
          writeLines("Your factor list is :")
          print(colnames(sample_data(p)))
          factor_2<-scan(n=1, what = character())
          
          lowlikely_remove<-data.frame(lowlikely, keep=NA)
          colnames(lowlikely_remove)[1]<-"OTU"
          writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
          for (T in 1:nrow(lowlikely)) {
            
            graphics.off()
            test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
            
            if (logb(max(test$value), base=3) != 0) { 
              print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
            
              writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
              print(tax_table(p)[taxa_names(p) == lowlikely_remove$OTU[T],])
              writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(lowlikely), sep=" "))
              switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
              
              if (keep == "High") {
                
                lowlikely_remove$keep[T]<- "2"
                
              } 
              if (keep == "Yes") {
                
                lowlikely_remove$keep[T]<- "1"
                
              }
              if (keep == "No") {
                
                lowlikely_remove$keep[T]<- "0"
                
              }
            } else {
              
              writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
            }
          }
        }
        
      }
      
      if (ncol(C12_otu_high) != 0) {
        
        #Add in factors
        C12_otu_high<-data.frame(C12_otu_high, as.vector(Factors))
        
        #Graph (Highly Likely OTUs as boxplots)
        melted_high<-melt(C12_otu_high)
        
        writeLines("Here is a quick look at how the total counts breakdown with the OTUs that have been selected as enriched")
        print(ddply(melted_high, ~ Enrichment, summarize, value=sum(value)))
        
        writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be C12 unenriched? There are ", nrow(highlikely),"of these OTUs?"))
        switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
        
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
            test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
            
            if (logb(max(test$value), base=3) != 0) { 
              print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
            
              writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
              print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
              writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
              switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
              
              if (keep == "High") {
                
                highlikely_remove$keep[T]<- "2"
                
              } 
              if (keep == "Yes") {
                
                highlikely_remove$keep[T]<- "1"
                
              }
              if (keep == "No") {
                
                highlikely_remove$keep[T]<- "0"
                
              }
            } else {
                            
              writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
                         
            }
          }
        }
      }
      
      if (ncol(C12_otu_high) != 0 & (ncol(C12_otu_low) != 0)) {
        
        ##Keep the keepers alongside the highlikely OTUs
        manual<-rbind(highlikely_remove, lowlikely_remove)
        manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
        manual$Subsetted_by<- paste(divide)
        
      }
      
      if (ncol(C12_otu_high) != 0 & (ncol(C12_otu_low) == 0)) {
        
        ##Keep the keepers alongside the highlikely OTUs
        manual<-rbind(highlikely_remove)
        manual<-cbind(manual, Likelihood_Enriched=rep("high", nrow(highlikely_remove)))
        manual$Subsetted_by<- paste(divide)
        
      }
      
      if (ncol(C12_otu_high) == 0 & (ncol(C12_otu_low) != 0)) {
        
        ##Keep the keepers alongside the highlikely OTUs
        manual<-rbind(lowlikely_remove)
        manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
        manual$Subsetted_by<- paste(divide)
        
      }
      
      } else { 
      
      manual <- data.frame(OTU=0,keep=0,Likelihood_Enriched=0,Subsetted_by=0)
      
    }
    
    assign(divide, manual)
      
    }
    
    manual<-do.call(rbind, lapply(unique_factor, as.symbol))
    
    if (any(duplicated(manual$OTU)) == TRUE) {
      manual<-manual[-which(duplicated(manual$OTU)),]
    }
    
    if (any(manual$OTU == 0) == TRUE) {
    manual<-manual[-which(manual$OTU == 0),]
    }
    
    ##Keep the keepers 
    write.csv(manual, file=paste("Output\\C12 OTUs - voom - Manual Edits_", VOOM_divide, "_", sign_save, "_", cutoff, ".csv", sep=""))
    C12_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]))
    enr_otu<-as.vector(manual[which(manual$keep >= 1),1])
      
    #Prune phyloseq to contain enriched OTUs
    p_C12_voom<-prune_taxa(C12_otu, p)
      
    print(p_C12_voom)
      
    saveRDS(p_C12_voom, file=paste("Output\\BackUps\\p_c12_OTU_voom_", VOOM_divide, "_", sign_save, "_", cutoff,".rds", sep=""))
    write.csv(taxa_names(p_C12_voom), file= paste("Output\\p_C12_OTU_voom_", VOOM_divide, "_", sign_save, "_", cutoff,".csv", sep=""))
    
  }

} 

if (length(VOOM_divide) == 0) {
  
  #Make design matrix
  Enrichment <-Factors[,colnames(Factors) == status]
  design <- model.matrix(~Enrichment)
  
  #Get count data in correct format
  rows<-nrow(otu_table(p_enr))
  cols<-ncol(otu_table(p_enr))
  counts<-as.matrix(otu_table(p_enr)[,])
  attributes(counts)<-NULL  #Remove phyloseq parameters
  counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
  rownames(counts)<-taxa_names(p_enr)   #Return row and column names
  colnames(counts)<-sample_names(p_enr)
  counts<-data.frame(counts)
  
  #the voom-limma process
  norm.factor <- calcNormFactors(counts)
  dat.voomed <- voom(counts, design, plot = TRUE, lib.size = colSums(counts) * norm.factor)
  fit <- lmFit(dat.voomed, design)
  fit <- eBayes(fit)
  
  TopHits<-topTable(fit, coef = paste("Enrichment", label, sep=""), sort.by = "logFC", adjust.method="BH", number=Inf)
  TopHits<-TopHits[with(TopHits, order(-logFC)), ]
  
  writeLines("\nLimma-voom will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C13.")
  writeLines("\nHere is the count data for your differentially abundant OTU.")
  print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
  
  writeLines("\nIs it showing abundance in your C13 samples?")
  switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
  
  if (order == "No") {
    TopHits$logFC<-TopHits$logFC*(-1)
    TopHits<-TopHits[with(TopHits, order(-logFC)), ]
    writeLines("\nNow it should!")
    print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
    
  }
  
  #Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
  enr_otu<-TopHits[which(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0),]
  highlikely<-which(enr_otu$logFC > 1.5)
  lowlikely<-which(enr_otu$logFC < 1.5)
  
  highlikely<-rownames(enr_otu)[highlikely]
  highlikely<-data.frame(highlikely)
  colnames(highlikely)<-"OTU"
  lowlikely<-rownames(enr_otu)[lowlikely]
  lowlikely<-data.frame(lowlikely)
  colnames(lowlikely)<-"OTU"
  
  enr_otu<-rownames(enr_otu)
  
  writeLines(paste("\nIn total you have recovered ", length(enr_otu)," enriched OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be enriched and ", nrow(lowlikely)," will have to be manually inspected and confirmed."))
  
  ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
  d=enr_otu
  OTU_location<- vector(mode="numeric", length=length(d))
  
  for (mercy in 1:length(d)) {
    
    OTU_location[mercy]<-grep(enr_otu[mercy], rownames(counts))
    
  }
  
  enr_otu<-counts[OTU_location,]
  
  #Make into data.frame
  enr_otu<-data.frame(enr_otu)
  enr_otu$likelihood<-NA
  
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
  writeLines(paste("\nWould you like to manually examine your", nrow(lowlikely),"marginal OTUs? If not, do consider applying your own filter to the data."))
  switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
  
  if (pain == "Yes") {
    writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
    writeLines("Your factor list is :")
    print(colnames(sample_data(p)))
    factor_2<-scan(n=1, what = character())
    
    lowlikely_remove<-data.frame(lowlikely, keep=NA)
    colnames(lowlikely_remove)[1]<-"OTU"
    writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
    
    for (T in 1:nrow(lowlikely)) {
      
      graphics.off()
      test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
      
      if (logb(max(test$value), base=3) != 0) { 
        print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value,base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(counts)"))
        writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
        print(tax_table(p_enr)[taxa_names(p_enr) == lowlikely_remove$OTU[T],])
        writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(lowlikely), sep=" "))
        switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
        
        if (keep == "High") {
          
          lowlikely_remove$keep[T]<- "2"
          
        } 
        if (keep == "Yes") {
          
          lowlikely_remove$keep[T]<- "1"
          
        }
        if (keep == "No") {
          
          lowlikely_remove$keep[T]<- "0"
          
        }
      } else {
        
        writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
        
      }
    }
  }
  
  writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be enriched? There are ", nrow(highlikely),"of these OTUs?"))
  switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
  
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
      test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
      
      if (logb(max(test$value), base=3) != 0) { 
        print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
        
        writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
        print(tax_table(p_enr)[taxa_names(p_enr) == highlikely_remove$OTU[T],])
        writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
        switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
        
        if (keep == "High") {
          
          highlikely_remove$keep[T]<- "2"
          
        } 
        if (keep == "Yes") {
          
          highlikely_remove$keep[T]<- "1"
          
        }
        if (keep == "No") {
          
          highlikely_remove$keep[T]<- "0"
          
        }
      } else {
        
        writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
                   
      }
    } 
  }
  
  if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE) {
    ##Keep the keepers alongside the highlikely OTUs
    manual<-rbind(highlikely_remove, lowlikely_remove)
    manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
    write.csv(manual, file=paste("Output\\C13 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
    enr_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]), as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1])))
    
  } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
    ##Keep the keepers alongside the highlikely OTUs
    manual<-lowlikely_remove
    manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
    write.csv(manual, file=paste("Output\\C13 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
    enr_otu<-as.vector(as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1]))
    
  } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
    ##Keep the keepers alongside the highlikely OTUs
    manual<-highlikely_remove
    manual<-cbind(manual, Likelihood_Enriched=rep("high", nrow(highlikely_remove)))
    write.csv(manual, file=paste("Output\\C13 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
    enr_otu<-as.vector(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]))    
  }
  
  #Prune phyloseq to contain enriched OTUs
  p_enr_voom<-prune_taxa(enr_otu, p)
  
  saveRDS(p_enr_voom, file=paste("Output\\BackUps\\p_Enriched_OTU_voom_", sign_save, "_", cutoff,".rds", sep=""))
  write.csv(taxa_names(p_enr_voom), file= paste("Output\\p_Enriched_OTU_voom_", sign_save, "_", cutoff,".csv", sep=""))
  
  ####
  #Repeat for C12 OTUs
  ####
  
  writeLines("Would you like to perform a similar analyses to verify which OTUs are more commonly present in the C12 Unenriched Samples? This may be useful in contrasting your findings.")
  switch(menu(c("Yes","No")), proceed<-("Yes"), proceed<-("No"))
  
  if (proceed == "Yes") {
    p_12C<-prune_samples(colnames(otu_table(p))[which(sample_data(p)[,status] != paste(label))], p)
    
    ###Pre-filtering steps
    #Determine which OTU are not found in the OTU table (after C12 sample have been removed)
    writeLines("\nFirst, we will first remove all OTUs which have zero representation in the C12 group.")
    total<-nrow(otu_table(p))
    number<-length(which(rowSums(otu_table(p_12C)) == 0))
    percent1<-format(round(((number/total)*100), 2), nsmall = 2)
    p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) == 0),]
    
    writeLines(paste("\nWe removed", number," OTUs at this step, which corresponds to", percent1,"percent of OTUs"))
    
    ##Filter based on abundance, maybe rowSum > 5
    writeLines("\nNext, we will filter by a minimum abundance value across all C12 samples, namely, the OTU must be found at an abundance > than the threshold in TOTAL for your C12 samples.\n\n Please select a threshold number. [Recommended: 5] ")
    threshold<-scan(n=1, what = numeric())
    number<-length(which(rowSums(otu_table(p_12C)) < threshold))
    p_12C<-otu_table(p_12C)[-which(rowSums(otu_table(p_12C)) < threshold),]
    percent2<-format(round(((number/total)*100), 2), nsmall = 2)
    writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent2,"percent of the original (pre-filtered) OTUs."))
    
    ##Filter based on occurrence in minimum number of samples > 3
    writeLines("\nNext we will filter by the minimum number of times an OTU must occur in C12 samples.\n\n Please select a minimum number. [Recommended: 3]")
    min_number<-scan(n=1, what = numeric())
    present_absent<-otu_table(p_12C)
    present_absent[present_absent > 0] <- 1
    number<-length(which(rowSums(present_absent) < min_number))
    p_12C<-otu_table(p_12C)[-which(rowSums(present_absent) < min_number),]
    percent3<-format(round(((number/total)*100), 2), nsmall = 2)
    remaining<-format(round(total*(1-((as.numeric(percent1)/100)+(as.numeric(percent2)/100)+(as.numeric(percent3)/100))),0))
    writeLines(paste("\nWe removed", number,"OTUs at this step, which corresponds to", percent3,"percent of the original (pre-filtered) OTUs. Combined, you have removed ", as.numeric(percent1)+as.numeric(percent2)+as.numeric(percent3)," percent of your OTUs, which leaves a total of", remaining,"OTUs in your C12 dataset."))
    
    #Return to original phyloseq object
    p_12C<-intersect(rownames(p_12C), taxa_names(p_12C))
    p_12C<-prune_taxa(p_12C, p)
    
    writeLines("\nYour data table now looks like this:\n")
    print(p_12C)
    
    #Make design matrix
    Enrichment <-Factors[,colnames(Factors) == status]
    design <- model.matrix(~Enrichment)
    
    #Only necessary due to conflicts with the formatting of a phyloSeq object
    rows<-nrow(otu_table(p_12C))
    cols<-ncol(otu_table(p_12C))
    counts<-as.matrix(otu_table(p_12C)[,])
    attributes(counts)<-NULL  #Remove phyloseq parameters
    counts<-matrix(counts, nrow=rows, ncol=cols)  #Create matrix
    rownames(counts)<-taxa_names(p_12C)   #Return row and column names
    colnames(counts)<-sample_names(p_12C)
    counts<-data.frame(counts)
    
    #the voom-limma process
    norm.factor <- calcNormFactors(counts)
    dat.voomed <- voom(counts, design, plot = TRUE, lib.size = colSums(counts) * norm.factor)
    fit <- lmFit(dat.voomed, design)
    fit <- eBayes(fit)
    
    TopHits<-topTable(fit, coef = paste("Enrichment", label, sep=""), sort.by = "logFC", adjust.method="BH", number=Inf)
    TopHits<-TopHits[with(TopHits, order(-logFC)), ]
    
    writeLines("\nLimma-voom will assign negative or positive values based on the initial ordering of your factor level, so, we will need to manually verify which OTUs are differentially abundant in favour of C13.")
    writeLines("\nHere is the count data for your differentially abundant OTU.")
    print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
    
    writeLines("\nIs it showing abundance in your C12 samples?")
    switch(menu(c("Yes","No")), order<-("Yes"), order<-("No"))
    
    if (order == "No") {
      TopHits$logFC<-TopHits$logFC*(-1)
      TopHits<-TopHits[with(TopHits, order(-logFC)), ]
      writeLines("\nNow it should!")
      print(otu_table(p)[rownames(otu_table(p)) == rownames(TopHits)[1]])
      
    }
    
    #Extract taxa names for all log change <-1.5 | padj <0.05 (it has to be done in this strange way because p-values are missing for many of the interesting OTU)
    C12_otu<-TopHits[which(TopHits$adj.P.Val < 0.05 & TopHits$logFC > 0),]
    highlikely<-which(C12_otu$logFC > 1.5)
    lowlikely<-which(C12_otu$logFC < 1.5)
    
    highlikely<-rownames(C12_otu)[highlikely]
    highlikely<-data.frame(highlikely)
    colnames(highlikely)<-"OTU"
    lowlikely<-rownames(C12_otu)[lowlikely]
    lowlikely<-data.frame(lowlikely)
    colnames(lowlikely)<-"OTU"
    
    C12_otu<-rownames(C12_otu)
    
    writeLines(paste("\nIn total you have recovered ", length(C12_otu)," enriched OTUs. Of these OTUs, ", nrow(highlikely), " are highly likely to be enriched and ", nrow(lowlikely)," is marginal and (unlike 13C selection) will be ignored."))
    
    ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
    d=C12_otu
    OTU_location<- vector(mode="numeric", length=length(d))
    
    for (mercy in 1:length(d)) {
      
      OTU_location[mercy]<-grep(C12_otu[mercy], rownames(counts))
      
    }
    
    C12_otu<-counts[OTU_location,]
    
    #Make into data.frame
    C12_otu<-data.frame(C12_otu)
    C12_otu$likelihood<-NA
    
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
    switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
    
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
        test<-melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),]
        
        if (logb(max(test$value), base=3) != 0) { 
          print(ggplot(melted_high[grep(highlikely_remove$OTU[T],melted_high$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
          
          writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
          print(tax_table(p)[taxa_names(p) == highlikely_remove$OTU[T],])
          writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(highlikely), sep=" "))
          switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
          
          if (keep == "High") {
            
            highlikely_remove$keep[T]<- "2"
            
          } 
          if (keep == "Yes") {
            
            highlikely_remove$keep[T]<- "1"
            
          }
          if (keep == "No") {
            
            highlikely_remove$keep[T]<- "0"
            
          }
        } else {
          
          writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
                     
        }
      }
    }
    
    writeLines(paste("\nWould you like to manually examine the OTUs which were deemed highly likely to be C12-only? There are ", nrow(lowlikely),"of these OTUs?"))
    switch(menu(c("Visualize and edit","Do it myself later (sometimes causes bugs)")), pain<-("Yes"), pain<-("No"))
    
    if (pain == "Yes") {
      writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
      writeLines("Your factor list is :")
      print(colnames(sample_data(p)))
      factor_2<-scan(n=1, what = character())
      
      lowlikely_remove<-data.frame(lowlikely, keep=NA)
      colnames(lowlikely_remove)[1]<-"OTU"
      writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
      
      for (T in 1:nrow(lowlikely)) {
        
        graphics.off()
        test<-melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),]
        
        if (logb(max(test$value), base=3) != 0) { 
          print(ggplot(melted_low[grep(lowlikely_remove$OTU[T],melted_low$variable),], aes(variable, logb(value, base=3), colour=get(status), shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02)) + ylab("Log3(Counts)"))
          
          writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
          print(tax_table(p)[taxa_names(p) == lowlikely_remove$OTU[T],])
          writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(lowlikely), sep=" "))
          switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
          
          if (keep == "High") {
            
            lowlikely_remove$keep[T]<- "2"
            
          } 
          if (keep == "Yes") {
            
            lowlikely_remove$keep[T]<- "1"
            
          }
          if (keep == "No") {
            
            lowlikely_remove$keep[T]<- "0"
            
          }
        } else {
          
          writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
          
        }
      }
    }
    
    if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == TRUE) {
        ##Keep the keepers 
        manual<-rbind(highlikely_remove, lowlikely_remove)
        manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove)), rep("low", nrow(lowlikely_remove))))
        write.csv(manual, file=paste("Output\\C12 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        C12_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1]), as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1])))
        
      } else if (exists("highlikely_remove") == FALSE & exists("lowlikely_remove") == TRUE) {
        ##Keep the keepers 
        manual<-lowlikely_remove
        manual<-cbind(manual, Likelihood_Enriched=rep("low", nrow(lowlikely_remove)))
        write.csv(manual, file=paste("Output\\C12 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        C12_otu<-as.vector(as.matrix(lowlikely_remove$OTU[lowlikely_remove$keep >= 1]))
        
      } else if (exists("highlikely_remove") == TRUE & exists("lowlikely_remove") == FALSE) {
        ##Keep the keepers 
        manual<-highlikely_remove
        manual<-cbind(manual, Likelihood_Enriched=c(rep("high", nrow(highlikely_remove))))
        write.csv(manual, file=paste("Output\\C12 OTUs - voom - Manual Edits_", sign_save, "_", cutoff, ".csv", sep=""))
        C12_otu<-as.vector(rbind(as.matrix(highlikely$OTU[highlikely_remove$keep >= 1])))
                
      }
      
      #Prune phyloseq to contain enriched OTUs
      p_C12_voom<-prune_taxa(C12_otu, p)
      
      print(p_C12_voom)
      
      saveRDS(p_C12_voom, file=paste("Output\\BackUps\\p_c12_OTU_voom_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(taxa_names(p_C12_voom), file= paste("Output\\p_C12_OTU_voom_", sign_save, "_", cutoff,".csv", sep=""))
    
  }
  
}

}

if (Method == "Relative Abundance") {
  
  #########
  #Using Relative Abundance to Identify Enriched OTUs
  #########
  
  ##Use proportion of total counts for each
  p_rel<-p  
  
  for (n in 1:ncol(otu_table(p_rel))) {
    
    otu_table(p_rel)[,n] <- otu_table(p_rel)[,n] / sum(otu_table(p_rel)[,n])
    
  }
  
  ##Calculate relative abundance in between paired C13 and C12 samples 
  
  #Simplifying terms for loop
  #Input which phyloseq object has your paired C13 and C12 samples
  d<-p_rel  ##The one piece of required information!!
  
  writeLines("Your sample labels are:")
  print(rownames(sample_data(d)))
  writeLines("\nPlease supply the label you've used in your sample names to distinguish enriched status.For example, \"LogLake01_13C\" would be \"13C\"\n")
  label_C13<-scan(n=1, what = character())
  writeLines("\nPlease supply the label you used to demarcate 12C enriched samples. Example : C12\n")
  label_C12<-scan(n=1, what = character())
  
  names<-colnames(otu_table(d))  #note that "d" is used in this term
  C13_samples<-colnames(otu_table(d))[grep(label_C13, names)] 
  C12_samples<-colnames(otu_table(d))[grep(label_C12, names)] 
  
  C12_sample_holder<- data.frame(C12_samples)
  C13_sample_holder<- data.frame(C13_samples)
  
  C13_samples<-gsub(label_C13, "", C13_samples)
  C12_samples<-gsub(label_C12, "", C12_samples) 
  
  C12_sample_holder$subbed<- C12_samples
  C13_sample_holder$subbed<- C13_samples
  
  p_rel<-matrix(nrow= nrow(otu_table(d)))
  
  #Must identify orphaned pairs before running main loop
  #Remove Sample with only C13 representatives
  C13_only<-setdiff(C13_samples, C12_samples) 
  
  #For orphaned 13C, we'll try to use an average of the the most related C12 samples as a reference
  
  if (length(C13_only) != 0) {
    
    for (adopt in C13_only) {
      
      writeLines(paste("\nThe following sample is orphaned:",adopt))
      writeLines("\nPlease input a string which may uniquely identify a group of related samples to use as reference for orphaned 13C samples.\n\nFeel free to use regular expression fit for R, such as [A-Z] or [0-9], to manage variation in sample titles.\n\nFor instance, \"LH0[0-9][1,3,5,7,9]\" will use the following fictitious samples \"LH007, LH009, LH011\" to reference \"BL005_C12\" against.")
      neighbour<-scan(n=1, what = character())
      avg_neighbour<-C12_samples[grep(neighbour, C12_samples)]
      
      for (nextdoor in avg_neighbour) {
        
        nextdoor_full <- C12_sample_holder[which(C12_sample_holder$subbed == nextdoor), 1]  #Return to original sample name
        foo<-otu_table(d)[,colnames(otu_table(d)) == nextdoor_full]
        assign(paste(nextdoor),foo)
        
      }
      
      avg_neighbour<-data.frame(do.call(cbind, lapply(avg_neighbour, as.symbol)))
      avg_neighbour<-data.frame(mean=rowMeans(avg_neighbour))
      
      ###NOTE: if your samples weren't named exactly SAMPLEID_13C, then you may get an error in the calls below which use "adopt". alter paste(adopt, "_13C") to match your actual sample names
      adopt_C13 <- C13_sample_holder[which(C13_sample_holder$subbed == adopt),1]
      
      ##Calculate relative abundance in C13 samples
      C13_rel<-otu_table(d)[,colnames(otu_table(d)) == adopt_C13] / avg_neighbour[,1]
      
      #Treat Inf as actual count data (i.e. count / one, rather than zero) 
      C13_rel[is.infinite(C13_rel)] <- otu_table(d)[which(is.infinite(C13_rel)), grep(adopt_C13, colnames(otu_table(d)))]
      
      #Set NaN to zero
      C13_rel[is.nan(C13_rel)] <- 0 
      
      #Set all fractions to zero (B/c these will be replaced with negative rel abundances)
      C13_rel[C13_rel[,1] < 1 & C13_rel[,1] > 0] <- 0 
      
      ##Calculate relative abundance in C12 samples
      C12_rel<-avg_neighbour[,1] / otu_table(d)[,colnames(otu_table(d)) == adopt_C13]
      
      #Repeat above transformations
      C12_rel[is.infinite(C12_rel)] <- avg_neighbour[which(is.infinite(C12_rel)),1]
      C12_rel[is.nan(C12_rel)] <- 0 
      C12_rel[C12_rel[,1] < 1 & C12_rel[,1] > 0] <- 0
      
      #Make C12 negative
      C12_rel[is.finite(C12_rel)] <- C12_rel[is.finite(C12_rel)]*-1
      
      #Merge tables
      list <- as.matrix(C13_rel + C12_rel)
      p_rel <- cbind(p_rel, list)
    }
  }
  
  #Remove Sample with only C12 representatives
  C12_only<-setdiff(C12_samples, C13_samples)
  
  if (length(C12_only) != 0) {
    
    for (kill in C12_only) {  
      
      d<-prune_samples(colnames(otu_table(d))[-grep(kill, colnames(otu_table(d)))], d)
      
    } 
  }
  
  ##NOW The main loop for all paired samples
  #Subset only paired samples
  paired_C13<-intersect(C13_sample_holder$subbed, C12_sample_holder$subbed) 
  
  #Run it
  for(C13 in paired_C13) {
    
    C12 <- C12_sample_holder[which(C12_sample_holder$subbed == C13),1]
    C13 <- C13_sample_holder[which(C13_sample_holder$subbed == C13),1]
    
    ##Calculate relative abundance in C13 samples
    C13_rel<-otu_table(d)[,grep(C13, names)] / otu_table(d)[,grep(C12, names)]
    
    #Treat Inf as actual count data (i.e. count / one, rather than zero) 
    C13_rel[is.infinite(C13_rel)] <- otu_table(d)[which(is.infinite(C13_rel)), grep(C13, colnames(otu_table(d)))]
    
    #Set NaN to zero
    C13_rel[is.nan(C13_rel)] <- 0 
    
    #Set all fractions to zero (B/c these will be replaced with negative rel abundances)
    C13_rel[C13_rel[,1] < 1 & C13_rel[,1] > 0] <- 0 
    
    ##Calculate relative abundance in C12 samples
    C12_rel<-otu_table(d)[,grep(C12, names)] / otu_table(d)[,grep(C13, names)]
    
    #Repeat above transformations
    if (any(is.infinite(C12_rel)) == TRUE) {
      
      C12_rel[is.infinite(C12_rel)] <- otu_table(d)[which(is.infinite(C12_rel)), grep(C12, colnames(otu_table(d)))]
      C12_rel[is.nan(C12_rel)] <- 0 
      C12_rel[C12_rel[,1] < 1 & C12_rel[,1] > 0] <- 0
      
    } 
    
    #Make C12 negative
    C12_rel[is.finite(C12_rel)] <- C12_rel[is.finite(C12_rel)]*-1
    
    #Merge tables
    list <- as.matrix(C13_rel + C12_rel)
    p_rel <- cbind(p_rel, list)
    
  }
    
  p_rel <- p_rel[,-1]
  
  saveRDS(p_rel, file=paste("Output\\BackUps\\p_rel_abund_OTU_", sign_save, "_", cutoff,".rds", sep=""))
  write.csv(taxa_names(p_rel), file= paste("Output\\p_rel_abund_OTU_", sign_save, "_", cutoff,".csv", sep=""))  
  
  ######
  # Select the enrOTU and c12 OTU in a similar visualization process
  ######
  
  #Subset Sample Variables for only 13C enriched (i.e. 12C controls are now built into 13C)
  design <- as.matrix(sample_data(prune_samples(colnames(p_rel), p)))
  
  #Transpose and create mother to then spring off
  main <- data.frame(cbind(t(p_rel), design))
  
  #Prompt for specific factor to break down rel. abund with
  writeLines("\nIn this pipeline, relative abundance is summed across all samples for each OTU.\n\nIf you would rather sum within a specific factor, please enter it now. If there is none, just hit ENTER.")
  writeLines("Your factor list is :")
  print(colnames(sample_data(p)))
  factor<-scan(n=1, what = character())
  
  writeLines("\nAre you returning to re-analyze the C12 (unenriched) community? If you would you like to skip the 13C OTU selection step, do so now.")
  switch(menu(c("Do 13C OTU selection", "Skip to 12C")), skip<-("No"), skip<-("Yes"))
  
  if (skip == "No") {
    
    #####
    #enrOTU selection
    #####
    
    if (length(factor) != 0) {
      
      writeLines("All OTUs with a sum of relative abundance greater than 3 are being carried forward as putatively enriched degraders.")
      factor_list<-grep(factor, colnames(main))
      unique_factor<-as.vector(main[which(duplicated(main[,factor_list]) == FALSE),factor_list])
      
      writeLines("Because we will subset the data for each factor, there may be cases were an OTU is only represented once within all samples. This can be undesirable, so if you would like to set a threshold for the number of samples an OTU must appear, ENTER it now (or enter 0, for no threshold).")
      threshold<-scan(n=1, what = numeric())
      
      for (go in 1:length(unique_factor)) {
        
        #Split main into single factor
        divide<-unique_factor[go]
        main_foo<-subset(main, get(factor) == paste(divide))
        
        #Transpose back into OTU x Sample table and make data numeric
        main_foo<-t(main_foo)
        foo_rownames<-rownames(main_foo)
        foo_colnames<-colnames(main_foo)
        main_foo<-apply(main_foo, 2, function(x) as.numeric(x))
        main_foo<-data.frame(main_foo)  
        rownames(main_foo)<-foo_rownames
        
        #Apply threshold criteria
        main_foo<-main_foo[which(apply(main_foo, 1, function(x) length(which(x != 0)) > threshold)),]
        
        #Construct a Filter based on standard deviation?
        
        #Rank Total based on SUM of relative abundance of specific factor
        Total_rank<-apply(main_foo, 1, function(x) sum(x))
        Total_rank<-Total_rank[rev(order(Total_rank))]
        if (length(which(is.na(Total_rank))) != 0) {
          
          Total_rank<-Total_rank[-which(is.na(Total_rank))]
        }  
        
        if (length(which(Total_rank == 0)) != 0) {
          
          Total_rank<-Total_rank[-which(Total_rank == 0)]
        }
        
        Total_rank<-Total_rank[which(Total_rank > 3)]
        assign(divide, Total_rank)
        
      }
      
      Total_rank<-do.call(c, lapply(unique_factor, as.symbol))
      
      if (any(duplicated(names(Total_rank))) == TRUE) {
        
        Total_rank<-Total_rank[-which(duplicated(names(Total_rank)))]
        
      }
      
    } else {
      
      writeLines("All OTUs with a sum of relative abundance greater than 3 are being carried forward as putatively enriched degraders.")
      
      #Rank Total based on SUM of relative abundance
      Total_rank<-apply(p_rel, 1, function(x) sum(x))
      Total_rank<-Total_rank[rev(order(Total_rank))]
      
      if (length(which(is.na(Total_rank))) != 0) {
        
        Total_rank<-Total_rank[-which(is.na(Total_rank))]
      }  
      
      if (length(which(Total_rank == 0)) != 0) {
        
        Total_rank<-Total_rank[-which(Total_rank == 0)]
      }
      Total_rank<-Total_rank[which(Total_rank > 3)]
      
    }
    
    ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
    ##Loop to parse the location of each designated OTU and pull that info into a new matrix
    
    d=Total_rank
    OTU_location<- vector(mode="numeric", length=length(d))
    
    for (mercy in 1:length(d)) {
      
      OTU_location[mercy]<-grep(names(d)[mercy], rownames(p_rel))
      
    }
    
    enr_otu<-otu_table(p)[OTU_location,]
    
    #Make into data.frame
    enr_otu<-data.frame(enr_otu)
    
    #Transpose and Add Factors
    enr_otu<-t(enr_otu)
    
    #Add in factors
    enr_otu<-data.frame(enr_otu, as.matrix(sample_data(p)))
    
    #Melt and Add together Relative abundance and Total Counts
    melted<-melt(enr_otu)
    
    #### Let User Manage the Number of OTUs
    writeLines("\nYou can now visually verify all of the OTUs that have been selected based on relative abundance.\n\nYou will be able to discard or keep each OTU on the spot.") 
    writeLines(paste("\nYou have recovered", length(Total_rank),"OTUs. Is this number manageable for visual confirmation, or would you like to set a different threshold for identifying putative C13 enriched OTU? [Note: if you don't plan on visually verifying, just proceed]."))
    switch(menu(c("Enter New Threshold","Proceed")), new_thresh<-("Yes"), new_thresh<-("No"))
    
    while (new_thresh == "Yes") {
      
      if (length(factor) != 0) {
        
        writeLines("Please select a new threshold for which the total sum of relative abundances across a sample must be below [note: original was set at 3]")
        new_threshold<-scan(n=1, what = numeric())
        
        factor_list<-grep(factor, colnames(main))
        unique_factor<-as.vector(main[which(duplicated(main[,factor_list]) == FALSE),factor_list])
        
        writeLines("Because we will subset the data for each factor, there may be cases were an OTU is only represented once within all samples. This can be undesirable, so if you would like to set a threshold for the number of samples an OTU must appear, ENTER it now (or enter 0, for no threshold).")
        threshold<-scan(n=1, what = numeric())
        
        for (go in 1:length(unique_factor)) {
          
          #Split main into single factor
          divide<-unique_factor[go]
          main_foo<-subset(main, get(factor) == paste(divide))
          
          #Transpose back into OTU x Sample table and make data numeric
          main_foo<-t(main_foo)
          foo_rownames<-rownames(main_foo)
          foo_colnames<-colnames(main_foo)
          main_foo<-apply(main_foo, 2, function(x) as.numeric(x))
          main_foo<-data.frame(main_foo)  
          rownames(main_foo)<-foo_rownames
          
          #Apply threshold criteria
          main_foo<-main_foo[which(apply(main_foo, 1, function(x) length(which(x != 0)) > threshold)),]
          
          #Construct a Filter based on standard deviation?
          
          #Rank Total based on SUM of relative abundance of specific factor
          Total_rank<-apply(main_foo, 1, function(x) sum(x))
          Total_rank<-Total_rank[rev(order(Total_rank))]
          if (length(which(is.na(Total_rank))) != 0) {
            
            Total_rank<-Total_rank[-which(is.na(Total_rank))]
          }  
          
          if (length(which(Total_rank == 0)) != 0) {
            
            Total_rank<-Total_rank[-which(Total_rank == 0)]
          }
          
          Total_rank<-Total_rank[which(Total_rank > new_threshold)]
          assign(divide, Total_rank)
          
        }
        
        Total_rank<-do.call(c, lapply(unique_factor, as.symbol))
        
        if (any(duplicated(names(Total_rank))) == TRUE) {
          
          Total_rank<-Total_rank[-which(duplicated(names(Total_rank)))]
          
        }
        
      } else {
        
        writeLines("Please select a new threshold for which the total sum of relative abundances across a sample must be below [note: original was set at 3]")
        new_threshold<-scan(n=1, what = numeric())
        
        #Rank Total based on SUM of relative abundance
        Total_rank<-apply(p_rel, 1, function(x) sum(x))
        Total_rank<-Total_rank[rev(order(Total_rank))]
        
        if (length(which(is.na(Total_rank))) != 0) {
          
          Total_rank<-Total_rank[-which(is.na(Total_rank))]
        }  
        
        if (length(which(Total_rank == 0)) != 0) {
          
          Total_rank<-Total_rank[-which(Total_rank == 0)]
        }
        Total_rank<-Total_rank[which(Total_rank > new_threshold)]
        
      }
      
      ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
      ##Loop to parse the location of each designated OTU and pull that info into a new matrix
      
      d=Total_rank
      OTU_location<- vector(mode="numeric", length=length(d))
      
      for (mercy in 1:length(d)) {
        
        OTU_location[mercy]<-grep(names(d)[mercy], rownames(p_rel))
        
      }
      
      enr_otu<-otu_table(p)[OTU_location,]
      
      #Make into data.frame
      enr_otu<-data.frame(enr_otu)
      
      #Transpose and Add Factors
      enr_otu<-t(enr_otu)
      
      #Add in factors
      enr_otu<-data.frame(enr_otu, as.matrix(sample_data(p)))
      
      #Melt and Add together Relative abundance and Total Counts
      melted<-melt(enr_otu)
      
      #### Change Status to Exit Loop
      writeLines(paste("\nYou have recovered", length(Total_rank),"OTUs. Is this number manageable for visual confirmation, or would you like to set a different threshold for identifying putative C12 OTU? [Note: if you don't plan on visually verifying, just proceed]."))
      switch(menu(c("Enter New Threshold","Proceed")), new_thresh<-("Yes"), new_thresh<-("No"))
      
    }
    
    #### Visualization loop
    writeLines("\nYou can now visually verify all of the OTUs that have been selected based on relative abundance.\n\nYou will be able to discard or keep each OTU on the spot.") 
    writeLines(paste("\nWould you like to manually examine your", length(Total_rank),"OTUs? If not, do consider applying your own filter to the data."))
    switch(menu(c("Visualize and edit","Do it myself later")), pain<-("Yes"), pain<-("No"))
    
    if (pain == "Yes") {
      writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
      writeLines("Your factor list is :")
      print(colnames(sample_data(p)))
      factor_2<-scan(n=1, what = character())
      
      Total_rank_remove<-data.frame(names(Total_rank), keep=NA)
      colnames(Total_rank_remove)[1]<-"OTU"
      writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
      
      for (T in 1:length(Total_rank)) {
        
        graphics.off()
        test<-melted[grep(Total_rank_remove$OTU[T],melted$variable),]
        
        if (logb(max(test$value), base=3) != 0) { 
          print(ggplot(melted[grep(Total_rank_remove$OTU[T], melted$variable),], aes(variable, logb(value,base=3), colour=Enrichment, shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02))) #+ facet_grid(. ~ Data_Type))
          
          writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
          print(tax_table(p)[taxa_names(p) == Total_rank_remove$OTU[T],])
          writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(Total_rank), sep=" "))
          switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
          
          if (keep == "High") {
            
            Total_rank_remove$keep[T]<- "2"
            
          }
          if (keep == "Yes") {
            
            Total_rank_remove$keep[T]<- "1"
            
          }
          if (keep == "No") {
            
            Total_rank_remove$keep[T]<- "0"
            
          }
        } else {
          writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
        }
      }
    }
    
    Total_rank_remove_CSV <- as.vector(Total_rank_remove)
    Total_rank_remove <- as.vector(Total_rank_remove[which(Total_rank_remove$keep >= 1),1])
    
    p_enr_rel<-prune_taxa(Total_rank_remove, p)
    
    if (exists("new_threshold") == TRUE) {
      
      saveRDS(p_enr_rel, file=paste("Output\\BackUps\\p_Enriched_OTU_relabund_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(Total_rank_remove_CSV, file=paste("Output\\C13 OTUs - RelAbund - Manual Edits_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
      write.csv(taxa_names(p_enr_rel), file= paste("Output\\p_Enriched_OTU_relabund_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff,".csv", sep=""))
      
      writeLines("Three files have been created in your Output and Output\\BackUps folder. One is a .rds file, and two are .csv. One of the .csv contains the full information (which can be used in Block02c), the other is designed to be used in Block05.\n\nYou HAVE Changed the THRESHOLD and this is included in the name.")
      
    } else {
      
      saveRDS(p_enr_rel, file=paste("Output\\BackUps\\p_Enriched_OTU_relabund_", factor, "_", sign_save, "_", cutoff,".rds", sep=""))
      write.csv(Total_rank_remove_CSV, file=paste("Output\\C13 OTUs - RelAbund - Manual Edits_", factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
      write.csv(taxa_names(p_enr_rel), file= paste("Output\\p_Enriched_OTU_relabund_", factor, "_", sign_save, "_", cutoff,".csv", sep=""))
      
      writeLines("Three files have been created in your Output and Output\\BackUps folder. One is a .rds file, and two are .csv. One of the .csv contains the full information (which can be used in Block02c), the other is designed to be used in Block05.")
      
    }
    
  }
  
  #####
  #C12 Unenriched OTU Selection
  #####
  
  if (length(factor) != 0) {
    
    writeLines("All OTUs with a sum of relative abundance less than -3 are being carried forward as putative non-cellulose using members.")
    factor_list<-grep(factor, colnames(main))
    unique_factor<-as.vector(main[which(duplicated(main[,factor_list]) == FALSE),factor_list])
    
    writeLines("Because we will subset the data for each factor, there may be cases were an OTU is only represented once within all samples. This can be undesirable, so if you would like to set a threshold for the number of samples an OTU must appear, ENTER it now (or enter 0, for no threshold).")
    threshold<-scan(n=1, what = numeric())
    
    for (go in 1:length(unique_factor)) {
      
      #Split main into single factor
      divide<-unique_factor[go]
      main_foo<-subset(main, get(factor) == paste(divide))
      
      #Transpose back into OTU x Sample table and make data numeric
      main_foo<-t(main_foo)
      foo_rownames<-rownames(main_foo)
      foo_colnames<-colnames(main_foo)
      main_foo<-apply(main_foo, 2, function(x) as.numeric(x))
      main_foo<-data.frame(main_foo)  
      rownames(main_foo)<-foo_rownames
      
      #Apply threshold criteria
      main_foo<-main_foo[which(apply(main_foo, 1, function(x) length(which(x != 0)) > threshold)),]
      
      #Construct a Filter based on standard deviation?
      
      #Rank Total based on SUM of relative abundance of specific factor
      Total_rank<-apply(main_foo, 1, function(x) sum(x))
      Total_rank<-Total_rank[order(Total_rank)]
      
      if (length(which(is.na(Total_rank))) != 0) {
        
        Total_rank<-Total_rank[-which(is.na(Total_rank))]
      }  
      
      if (length(which(Total_rank == 0)) != 0) {
        
        Total_rank<-Total_rank[-which(Total_rank == 0)]
      }
      
      Total_rank<-Total_rank[which(Total_rank < -3)]
      assign(divide, Total_rank)
      
    }
    
    Total_rank<-do.call(c, lapply(unique_factor, as.symbol))
    
    if (any(duplicated(names(Total_rank))) == TRUE) {
      
      Total_rank<-Total_rank[-which(duplicated(names(Total_rank)))]
      
    }
    
  } else {
    
    writeLines("All OTUs with a sum of relative abundance less than -3 are being carried forward as putative non-cellulose using members.")
    
    #Rank Total based on SUM of relative abundance
    Total_rank<-apply(p_rel, 1, function(x) sum(x))
    Total_rank<-Total_rank[order(Total_rank)]
    
    if (length(which(is.na(Total_rank))) != 0) {
      
      Total_rank<-Total_rank[-which(is.na(Total_rank))]
    }  
    
    if (length(which(Total_rank == 0)) != 0) {
      
      Total_rank<-Total_rank[-which(Total_rank == 0)]
    }
    
    #Apply criteria (ie. sum must be less than) [Note: I tried variuos approaches based on mean and median, and due to the variability across samples, none were very effective. Still, it may be worthe exploring additional methods.]
    Total_rank<-Total_rank[which(Total_rank < -3)]
    
  }
  
  ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
  ##Loop to parse the location of each designated OTU and pull that info into a new matrix
  
  d=Total_rank
  OTU_location<- vector(mode="numeric", length=length(d))
  
  for (mercy in 1:length(d)) {
    
    OTU_location[mercy]<-grep(names(d)[mercy], rownames(p_rel))
    
  }
  
  C12_otu<-otu_table(p)[OTU_location,]
  
  #Make into data.frame
  C12_otu<-data.frame(C12_otu)
  
  #Transpose and Add Factors
  C12_otu<-t(C12_otu)
  
  #Add in factors
  C12_otu<-data.frame(C12_otu, as.matrix(sample_data(p)))
  
  #Melt and Add together Relative abundance and Total Counts
  melted<-melt(C12_otu)
  
  #### Let User Manage the Number of OTUs
  writeLines("\nYou can now visually verify all of the OTUs that have been selected based on relative abundance.\n\nYou will be able to discard or keep each OTU on the spot.") 
  writeLines(paste("\nYou have recovered", length(Total_rank),"OTUs. Is this number manageable for visual confirmation, or would you like to set a different threshold for identifying putative C12 OTU? [Note: if you don't plan on visually verifying, just proceed]."))
  switch(menu(c("Enter New Threshold","Proceed")), new_thresh<-("Yes"), new_thresh<-("No"))
  
  while (new_thresh == "Yes") {
    
    if (length(factor) != 0) {
      
      writeLines("Please select a new threshold for which the total sum of relative abundances across a sample must be below. [Note: original was set at -3 (i.e. negative), since C12 are considered depleted relative to C13]")
      new_threshold<-scan(n=1, what = numeric())
      
      factor_list<-grep(factor, colnames(main))
      unique_factor<-as.vector(main[which(duplicated(main[,factor_list]) == FALSE),factor_list])
      
      writeLines("Because we will subset the data for each factor, there may be cases were an OTU is only represented once within all samples. This can be undesirable, so if you would like to set a threshold for the number of samples an OTU must appear, ENTER it now (or enter 0, for no threshold).")
      threshold<-scan(n=1, what = numeric())
      
      for (go in 1:length(unique_factor)) {
        
        #Split main into single factor
        divide<-unique_factor[go]
        main_foo<-subset(main, get(factor) == paste(divide))
        
        #Transpose back into OTU x Sample table and make data numeric
        main_foo<-t(main_foo)
        foo_rownames<-rownames(main_foo)
        foo_colnames<-colnames(main_foo)
        main_foo<-apply(main_foo, 2, function(x) as.numeric(x))
        main_foo<-data.frame(main_foo)  
        rownames(main_foo)<-foo_rownames
        
        #Apply threshold criteria
        main_foo<-main_foo[which(apply(main_foo, 1, function(x) length(which(x != 0)) > threshold)),]
        
        #Construct a Filter based on standard deviation?
        
        #Rank Total based on SUM of relative abundance of specific factor
        Total_rank<-apply(main_foo, 1, function(x) sum(x))
        Total_rank<-Total_rank[order(Total_rank)]
        
        if (length(which(is.na(Total_rank))) != 0) {
          
          Total_rank<-Total_rank[-which(is.na(Total_rank))]
        }  
        
        if (length(which(Total_rank == 0)) != 0) {
          
          Total_rank<-Total_rank[-which(Total_rank == 0)]
        }
        
        Total_rank<-Total_rank[which(Total_rank < new_threshold)]
        assign(divide, Total_rank)
        
      }
      
      Total_rank<-do.call(c, lapply(unique_factor, as.symbol))
      
      if (any(duplicated(names(Total_rank))) == TRUE) {
        
        Total_rank<-Total_rank[-which(duplicated(names(Total_rank)))]
        
      }
      
      
    } else {
      
      writeLines("Please select a new threshold for which the total sum of relative abundances across a sample must be below. [Note: original was set at -3 (i.e. negative), since C12 are considered depleted relative to C13]")
      new_threshold<-scan(n=1, what = numeric())
      
      #Rank Total based on SUM of relative abundance
      Total_rank<-apply(p_rel, 1, function(x) sum(x))
      Total_rank<-Total_rank[order(Total_rank)]
      
      if (length(which(is.na(Total_rank))) != 0) {
        
        Total_rank<-Total_rank[-which(is.na(Total_rank))]
      }  
      
      if (length(which(Total_rank == 0)) != 0) {
        
        Total_rank<-Total_rank[-which(Total_rank == 0)]
      }
      
      #Apply criteria (ie. sum must be less than) [Note: I tried variuos approaches based on mean and median, and due to the variability across samples, none were very effective. Still, it may be worthe exploring additional methods.]
      Total_rank<-Total_rank[which(Total_rank < new_threshold)]
      
    }
    
    ###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
    ##Loop to parse the location of each designated OTU and pull that info into a new matrix
    
    d=Total_rank
    OTU_location<- vector(mode="numeric", length=length(d))
    
    for (mercy in 1:length(d)) {
      
      OTU_location[mercy]<-grep(names(d)[mercy], rownames(p_rel))
      
    }
    
    C12_otu<-otu_table(p)[OTU_location,]
    
    #Make into data.frame
    C12_otu<-data.frame(C12_otu)
    
    #Transpose and Add Factors
    C12_otu<-t(C12_otu)
    
    #Add in factors
    C12_otu<-data.frame(C12_otu, as.matrix(sample_data(p)))
    
    #Melt and Add together Relative abundance and Total Counts
    melted<-melt(C12_otu)
    
    #### Change Status
    writeLines(paste("\nYou have recovered", length(Total_rank),"OTUs. Is this number manageable for visual confirmation, or would you like to set a different threshold for identifying putative C12 OTU? [Note: if you don't plan on visually verifying, just proceed]."))
    switch(menu(c("Enter New Threshold","Proceed")), new_thresh<-("Yes"), new_thresh<-("No"))
    
  }
  
  ##Visualization Loop
  writeLines(paste("\nWould you like to manually examine your", length(Total_rank),"OTUs? If not, do consider applying your own filter to the data."))
  switch(menu(c("Visualize and edit","Do it myself later")), pain<-("Yes"), pain<-("No"))
  
  if (pain == "Yes") {
    writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
    writeLines("Your factor list is :")
    print(colnames(sample_data(p)))
    factor_2<-scan(n=1, what = character())
    
    Total_rank_remove<-data.frame(names(Total_rank), keep=NA)
    colnames(Total_rank_remove)[1]<-"OTU"
    writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
    
    for (T in 1:length(Total_rank)) {
      
      graphics.off()
      test<-melted[grep(Total_rank_remove$OTU[T],melted$variable),]
      
      if (logb(max(test$value), base=3) != 0) { 
        print(ggplot(melted[grep(Total_rank_remove$OTU[T], melted$variable),], aes(variable, logb(value,base=3), colour=Enrichment, shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02))) #+ facet_grid(. ~ Data_Type))
        
        writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
        print(tax_table(p)[taxa_names(p) == Total_rank_remove$OTU[T],])
        writeLines(paste("\nYou have curated", T, "OTUs out of a total of", length(Total_rank), sep=" "))
        switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
        
        if (keep == "High") {
          
          Total_rank_remove$keep[T]<- "2"
          
        }
        if (keep == "Yes") {
          
          Total_rank_remove$keep[T]<- "1"
          
        }
        if (keep == "No") {
          
          Total_rank_remove$keep[T]<- "0"
          
        }
      } else {
        
        writeLines(paste("\nThis OTU has a maximum count of 1, and would break the graph if you visualized it")) 
      }
    }
  }
  
  Total_rank_remove_CSV <- as.vector(Total_rank_remove)
  Total_rank_remove <- as.vector(Total_rank_remove[which(Total_rank_remove$keep >= 1),1])
  
  p_C12<-prune_taxa(Total_rank_remove, p)
  
  if (exists("new_threshold") == TRUE) {
    
    saveRDS(p_C12, file=paste("Output\\BackUps\\p_C12_OTU_relabund_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff,".rds", sep=""))
    write.csv(Total_rank_remove_CSV, file=paste("Output\\C12 OTUs - RelAbund - Manual Edits_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
    write.csv(taxa_names(p_C12), file= paste("Output\\p_C12_OTU_relabund_thresh_", new_threshold, "_", factor, "_", sign_save, "_", cutoff,".csv", sep=""))
    
    writeLines("Three files have been created in your Output and Output\\BackUps folder. One is a .rds file, and two are .csv. One of the .csv contains the full information (which can be used in Block02c), the other is designed to be used in Block05.\n\nYou HAVE Changed the THRESHOLD and this is included in the name.")
      
  } else {
    
    saveRDS(p_C12, file=paste("Output\\BackUps\\p_C12_OTU_relabund_", factor, "_", sign_save, "_", cutoff,".rds", sep=""))
    write.csv(Total_rank_remove_CSV, file=paste("Output\\C12 OTUs - RelAbund - Manual Edits_", factor, "_", sign_save, "_", cutoff, ".csv", sep=""))
    write.csv(taxa_names(p_C12), file= paste("Output\\p_C12_OTU_relabund_", factor, "_", sign_save, "_", cutoff,".csv", sep=""))
    
    writeLines("Three files have been created in your Output and Output\\BackUps folder. One is a .rds file, and two are .csv. One of the .csv contains the full information (which can be used in Block02c), the other is designed to be used in Block05.")
      
  }
}
library(limma)
library(phyloseq)

#Set Names
cutoff <- "0.01"
sign_save <- "Bacteria"
added_factor <- 

#Import Files
DESeq_13C <- read.csv(file= file.choose())
DESeq_12C <- read.csv(file= file.choose())

Voom_13C <- read.csv(file= file.choose())
Voom_12C <- read.csv(file= file.choose())

Rel_13C <- read.csv(file= file.choose())
Rel_12C <- read.csv(file= file.choose())

C13_List <- c("DESeq_13C", "Voom_13C", "Rel_13C")
C12_List <- c("DESeq_12C", "Voom_12C", "Rel_12C")

C13_List_no_Rel <- c("DESeq_13C", "Voom_13C")
C12_List_no_Rel <- c("DESeq_12C", "Voom_12C")

####
#Compare Efficacy of Each Method
####

#Compare Number of "High" to "Low" Assigned Based on Threshold (not really comparable since 1.25 was used for DESeq and 1.5 for Limma)

for (i in C13_List_no_Rel) {
  
  foo<-data.frame(get(i))
  total<- length(foo$Likelihood_Enriched)
  total_high <- length(which(foo$Likelihood_Enriched == "high"))
  percent_high <- total_high / total
  stat <- data.frame(method=paste(i), total=total, total_high=total_high, percent_high=percent_high)
  
  assign(paste(i,"_temp", sep=""), stat)
  
}

High_Lo_stat<-do.call(rbind, lapply(paste(C13_List_no_Rel,"_temp", sep=""), as.symbol))

#Compare Number of Manual Rejections

for (i in C13_List) {
  
  foo<-data.frame(get(i))
  zeros<-length(which(foo$keep == 0))
  ones<- length(which(foo$keep == 1))
  twos<- length(which(foo$keep == 2))
  total <- length(foo$keep)
    
  stat <- data.frame(method=paste(i), percent_removed=round((zeros/total)*100, 1), percent_moderate_certainty=round((ones/total)*100, 1), percent_high_certainty=round((twos/total)*100, 1))
  
  assign(paste(i,"_temp", sep=""), stat)
  
}

Visualized_keepers_stat<-do.call(rbind, lapply(paste(C13_List,"_temp", sep=""), as.symbol))

####
#Compare Amount of Overlap
####

#By OTU name
#List of all non-redundant OTUs 
OTU_Total <- rbind(DESeq_13C[,1:3], Voom_13C[,1:3], Rel_13C[,1:3])
OTU_Total <- OTU_Total[OTU_Total$keep != 0, ]

if (any(duplicated(OTU_Total$OTU))) {
  OTU_Total <- OTU_Total[-which(duplicated(OTU_Total$OTU)),]
}

#Identify which OTUs are identified by which method
DESeq<- matrix(c(rep(NA, length(OTU_Total$OTU))))
Voom<- matrix(c(rep(NA, length(OTU_Total$OTU))))
Rel<- matrix(c(rep(NA, length(OTU_Total$OTU))))

for (x in 1:length(OTU_Total$OTU)) {
  
  DESeq[x]<- length(intersect(OTU_Total$OTU[x], DESeq_13C$OTU)) == 1
  Voom[x]<- length(intersect(OTU_Total$OTU[x], Voom_13C$OTU)) == 1
  Rel[x]<- length(intersect(OTU_Total$OTU[x], Rel_13C$OTU)) == 1
  
}

##Return to Logical Vectors
DESeq <- as.logical(DESeq)
Voom <- as.logical(Voom)
Rel <- as.logical(Rel)

#Make FIRST version of Venn Diagram (MAKE SECOND ONE FOLLOWING REFINEMENT)
venn <- cbind(DESeq, Voom, Rel)
venn_counts <- vennCounts(venn)

vennDiagram(venn_counts, include = "both", 
            names = c("DESeq", "Limma-Voom", "Relative Abundance"), 
            cex = 1, counts.col = "red")

####
#Refining Process to Reach Final enrOTU list
####

#Re-Visualize all OTU identified by only a single method 
re_vis <- data.frame(OTU= rep("OTU", nrow(venn)), Source= rep("Method", nrow(venn)), stringsAsFactors=FALSE)

for (oi in 1:nrow(venn)) {
  
  if (length(which(venn[oi,] == TRUE)) == 1) {
    
    re_vis$OTU[oi] <- as.vector(OTU_Total$OTU[oi])
    re_vis$Source[oi] <- colnames(venn)[which(venn[oi,] == TRUE)]
      
  }
}

re_vis<-re_vis[-which(re_vis$OTU == "OTU" | re_vis$Source == "Method"),]

###The following (convoluted script) will prepare the data to visualize the the differnece between enriched and unenriched OTUs
##Loop to parse the location of each designated OTU and pull that info into a new matrix

#ENTER original phyloseq object
p<-readRDS(file.choose())

d=re_vis
OTU_location<- vector(mode="numeric", length=nrow(d))

for (mercy in 1:nrow(d)) {
  
  OTU_location[mercy]<-grep(d$OTU[mercy], rownames(otu_table(p)))
  
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

#### Visualization loop
writeLines("\nIs there any other factor you'd like to visualize on your graph? The symbol shape can be made to correspond to a second factor.\nPlease enter it now, if there is none, just hit ENTER.")
writeLines("Your factor list is :")
print(colnames(sample_data(p)))
factor_2<-scan(n=1, what = character())
  
re_vis_remove<-data.frame(OTU=re_vis$OTU, keep=NA)
writeLines("\nPlease note that the location of the points has arbitrarily been slightly staggered on the x-axis, so symbols do not overlap each other")
  
for (T in 1:nrow(re_vis)) {
    
  graphics.off()
  print(ggplot(melted[grep(re_vis_remove$OTU[T], melted$variable),], aes(variable, logb(value,base=3), colour=Enrichment, shape=get(factor_2))) + geom_jitter(size = 5, position = position_jitter(width = 0.02))) #+ facet_grid(. ~ Data_Type))
  
  writeLines("Please select whether you would like to keep or discard this OTU.\n\nIt's corresponding taxonomy is:")
  print(tax_table(p)[taxa_names(p) == re_vis_remove$OTU[T],])
  writeLines(paste("\nYou have curated", T, "OTUs out of a total of", nrow(re_vis), sep=" "))
  writeLines(paste("\nIt was selected as enriched by ONLY the", re_vis$Source[T], "method", sep=" "))
  
  switch(menu(c("Highly Enriched", "Moderately Enriched","Discard")), keep<-("High"), keep<-("Yes"), keep<-("No"))
    
  if (keep == "High") {
    
    re_vis_remove$keep[T]<- "2"
      
  }
  if (keep == "Yes") {
      
    re_vis_remove$keep[T]<- "1"
      
  }
  if (keep == "No") {
      
    re_vis_remove$keep[T]<- "0"
      
  }    
}

####
##MAKE FINAL enrOTU phyloseq object and other stats
####

#STAT Showing Number of Removed Upon 2nd Review
Review_remove_stat<-data.frame(Percent_Total_Removed = round((length(which(re_vis_remove$keep == 0)) / nrow(venn))*100, 1), Remaining_Total = nrow(venn)-length(which(re_vis_remove$keep == 0)))

#Remove OTUs rejected upon second review
for (takedown in re_vis_remove$OTU[which(re_vis_remove$keep == 0)]) {
  
  OTU_Total <- OTU_Total[-grep(takedown, OTU_Total$OTU),]
  
}

#Prune phyloseq object TO FINAL VERSION OF enrOTU
p_enr_gold <- prune_taxa(as.vector(OTU_Total$OTU), p)

#STAT - Proportion of enrOTU composing the total diversity?
#How many counts required to justify presence of OTU?
Threshold<-10

Total_enrOTU<- nrow(otu_table(p_enr_gold))
Total_Taxa_in_Libraries <- length(which(rowSums(otu_table(p)) > Threshold))

Proportion_Cellulose_Degrading_Taxa_stat <- data.frame(Threshold_Counts_To_Qualify_OTU = Threshold, Percent_of_Total = round((Total_enrOTU/Total_Taxa_in_Libraries)*100, 2))


#####
#Make Venn Diagram with REFINED DATA
#####

#Identify which OTUs are identified by which method
DESeq<- matrix(c(rep(NA, length(OTU_Total$OTU))))
Voom<- matrix(c(rep(NA, length(OTU_Total$OTU))))
Rel<- matrix(c(rep(NA, length(OTU_Total$OTU))))

for (x in 1:length(OTU_Total$OTU)) {
  
  DESeq[x]<- length(intersect(OTU_Total$OTU[x], DESeq_13C$OTU)) == 1
  Voom[x]<- length(intersect(OTU_Total$OTU[x], Voom_13C$OTU)) == 1
  Rel[x]<- length(intersect(OTU_Total$OTU[x], Rel_13C$OTU)) == 1
  
}

##Return to Logical Vectors
DESeq <- as.logical(DESeq)
Voom <- as.logical(Voom)
Rel <- as.logical(Rel)

#Make FIRST version of Venn Diagram (MAKE SECOND ONE FOLLOWING REFINEMENT)
venn <- cbind(DESeq, Voom, Rel)
venn_counts <- vennCounts(venn)

vennDiagram(venn_counts, include = "both", 
            names = c("DESeq", "Limma-Voom", "Relative Abundance"), 
            cex = 1, counts.col = "red")
####
#Compare Taxonomic Profile According to Each Method
####
Taxa <- matrix(nrow=nrow(OTU_Total), ncol=7)
colnames(Taxa)<-colnames(tax_table(p))
rownames(Taxa)<-OTU_Total$OTU

for (q in 1:nrow(Taxa)) {
  
  Taxa[q,] <- tax_table(p)[grep(OTU_Total$OTU[q], rownames(tax_table(p))),]
  
}


#Remove Bootstrap information and taxonomic level place holder ("g__")
Taxa[]<-gsub("g__", "", Taxa[])
Taxa[]<-gsub("([0-9][0-9][0-9])", "", Taxa[])
Taxa[]<-gsub("([0-9][0-9])", "", Taxa[])
Taxa[]<-gsub("[(]", "", Taxa[])
Taxa[]<-gsub("[)]", "", Taxa[])

#Remove Species (too short a sequence to trust a species classification)
Taxa<-data.frame(Taxa)
Taxa$Species<-NULL

#Remove Redundancy in Taxonomic Information 
Taxa_collapsed<-Taxa[-which(duplicated(Taxa)),]
venn_collapsed<-venn[-which(duplicated(Taxa)),] 

#Identify which are found in only one data set
One_Method_Only_Taxa<-Taxa_collapsed[which(rowSums(venn_collapsed) == 1),]

for (q in 1:nrow(One_Method_Only_Taxa)) {
  
  One_Method_Only_Taxa$Method[q] <- re_vis[grep(rownames(One_Method_Only_Taxa)[q], re_vis$OTU), 2]

}

#Use "Order" as Taxonomic Level for Comparison (some phylum don't have characterized at the family level)
write.csv(One_Method_Only_Taxa, file = "Output\\Final\\Taxa_Exclusive_to_Single_Selection_Method.csv")

#Identify which are found in more than one data set
Multi_Method_Taxa<-Taxa_collapsed[which(rowSums(venn_collapsed) > 1),]

for (q in 1:nrow(Multi_Method_Taxa)) {
  
  Multi_Method_Taxa$Method[q] <- re_vis[grep(rownames(Multi_Method_Taxa)[q], re_vis$OTU), 2]
  
}

#Use "Order" as Taxonomic Level for Comparison (some phylum don't have characterized at the family level)
write.csv(Multi_Method_Taxa, file = "Output\\Final\\Taxa_Multiple_Selection_Methods.csv")
write.csv(Taxa_collapsed, file = "Output\\Final\\Taxa_ALL_Selection_Methods.csv")


############################## NOT FINISHED
############################## NOT FINISHED NOT FINISHED
############################## NOT FINISHED NOT FINISHED NOT FINISHED
############################## NOT FINISHED NOT FINISHED
############################## NOT FINISHED
####
#Calculate False-positive Rate (somehow) OR Compare Occurrence of same Taxonomy Between C12 and C13 (use C12 as )
####

#Theory: Random taxa have the same probability of being detected as C12 enriched, so the overlap can serve as a liberal control (i.e. we don't have direct method for attributing false-positive) 
#Noticeably, many of the taxa unique to a specific selection method are also highly common in the C12 libraries.

#Mark OTUs with taxonomy shared by C12 and C13 lists that were identified only once or twice as putative false-positives

####
#Write Files
####

#Stat Files
dir.create("Output\\Final")
write.csv(Visualized_keepers_stat, file="Output\\Final\\Stats_Quality_of_Selection_Method_Manual.csv")
write.csv(High_Lo_stat, file="Output\\Final\\Stats_Quality_of_Selection_Method_Automated.csv")
write.csv(Review_remove_stat, file="Output\\Final\\Stats_Percent_OTU_Removed_by_Review.csv")
write.csv(Proportion_Cellulose_Degrading_Taxa_stat, file="Output\\Final\\Stats_Proportion_of_OTU_Identified_as_Enriched.csv")

#enrOTU List File
if (exists("added_factor") == TRUE) {
  
  saveRDS(p_enr_gold, file=paste("Output\\BackUps\\p_Enriched_OTU_Gold_", added_factor, "_", sign_save, "_", cutoff,".rds", sep=""))
  write.csv(taxa_names(p_enr_gold), file= paste("Output\\Final\\p_Enriched_OTU_Gold_", added_factor, "_", sign_save, "_", cutoff,".csv", sep=""))
  
} else {
  
  saveRDS(p_enr_gold, file=paste("Output\\BackUps\\p_Enriched_OTU_Gold_", sign_save, "_", cutoff,".rds", sep=""))
  write.csv(taxa_names(p_enr_gold), file= paste("Output\\Final\\p_Enriched_OTU_Gold_", sign_save, "_", cutoff,".csv", sep=""))
    
}





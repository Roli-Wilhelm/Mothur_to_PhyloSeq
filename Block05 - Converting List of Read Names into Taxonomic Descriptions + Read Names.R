##Import the results from greping the representative enriched OTU namess (from p_enr) agains the original names files
##This script will then label all of the reads with "SIP_Cellulose" enriched
writeLines("Please select the \".csv\" file containing all names of reads associated with your enriched OTUs. Note: if you get an error saying problem with matching, try re-saving your .csv using excel as a .csv.")
x<-read.csv(file=file.choose(), header=F, stringsAsFactors = F)
x[,1]<-NULL

writeLines("Please input the name of your enrichment status. Example, \"SIP_Cellulose\".")
name<-scan(n=1, what = character())

writeLines("Please select the phyloseq file which contains your enriched OTU data. Example, ~\"\\Output\\BackUps\\p_All_Enriched_OTU_Bacteria_0.01.rds\".")
p<-readRDS(file=file.choose())

#Grep the taxonomy information from tax_table and condense to important info

for (count in 1:nrow(x)) {
  
  #Input taxa data
  tax_info<-tax_table(p)[grep(x[count,1], rownames(tax_table(p))),]
  
  #Remove bootstrap values and other extraneous stuff
  tax_info<-gsub("k__|p__|c__|o__|f__|g__|s__", "", tax_info)
  tax_info<-gsub("[0-9][0-9][0-9]", "", tax_info)
  tax_info<-gsub("[0-9][0-9]", "", tax_info)
  tax_info<-gsub("[(]", "", tax_info)
  tax_info<-gsub("[)]", "", tax_info)
  tax_info<-gsub("[0-9][0-9]", "", tax_info)
  tax_info[,][which(tax_info[,] == "")]<-"unclassified"
    
  if (length(grep("unclassified", tax_info)) != 0) {
    
    for (join in 1:ncol(x)) {
      
      if (x[count, join] != "") {
        
        x[count, join]<-paste(">", x[count, join], " ", name, "_Enriched", count, " Unclassified ", tax_info[,min(grep("unclassified", tax_info))-1]," partial 16S rRNA gene", sep="")
      
      }
    }
  } else {
    
    for (join in 1:ncol(x)) {
      
      if (x[count, join] != "") {
        
        x[count, join]<-paste(">", x[count, join], " ", name, "_Enriched", count," ", tax_info[,ncol(tax_info)]," partial 16S rRNA gene", sep="")  
      
      }
    }
  }
}

writeLines(paste("A .csv with your re-named reads has been created in the output directory entitled: \"Fasta_enr_", name, ".csv", sep=""))
o<-as.vector(as.matrix(x))
o<-o[o != ""]
write.csv(o, file=paste("Output\\Fasta_enr_", name, ".csv", sep=""), quote=F,  row.names=F)


     

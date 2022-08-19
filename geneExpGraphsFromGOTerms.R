#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


#Script for making gene expression graphs from gene lists
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)


#setwd("/Users/joannachustecki/Documents/RNAseqAnalysis/")
#CWTvsVWT <- read.csv("VRN2_RNA_seq_data_Untreated_vs_Treated/CWTvsVWT_0FC.csv")
#Cprt6vsVprt6 <- read.csv("VRN2_RNA_seq_data_Untreated_vs_Treated/Cprt6vsVprt6_0FC.csv")

#comp1.name <- "VRN2_RNA_seq_data_Untreated_vs_Treated/cPRT6vsvPRT6_0FC.csv"
#comp2.name <- "VRN2_RNA_seq_data_Untreated_vs_Treated/CDMvsVDM_0FC.csv"

comp1.name <- args[1]
comp2.name <- args[2]

comp1 <- read.csv(paste(getwd(),comp1.name,sep="/"))
comp2 <- read.csv(paste(getwd(),comp2.name,sep="/"))

#with our particular data set, the WT comparison data frame included two extra FRI.C1/V1 samples, affecting the normalisation. All 
# of our other data frames have the same values across the dataframe, so it won't matter, just as long as we're not using the WT frame for gathering read count data.
#So I've introduced a specific catch for that data frame
if(comp1.name == "VRN2_RNA_seq_data_AllTogether/cWT_vs_vWT_0FC.csv"){
  comp1 <- comp2
  print("CWTvsVWT in comparison")
}
if(comp2.name == "VRN2_RNA_seq_data_AllTogether/cWT_vs_vWT_0FC.csv"){
  comp2 <- comp1
  print("CWTvsVWT in comparison")
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#For making gene expression graphs from GO terms

overallComparison<- args[3]
specificComparison <- args[8]
#overallComparison<- "cPRT6_vs_vPRT6_cDM_vs_vDM" 
#specificComparison <- "cDM_vs_vDMintersectcPRT6_vs_vPRT6"

FCweFilteredBy<- as.numeric(args[9])

GOTermsQuantity <- read.delim(paste("GOTermAnalysis/",overallComparison,"/",specificComparison,"/result_GO-geneList-",specificComparison,"_",FCweFilteredBy,"FC.csv_BP_Fisher+KSclassic_test.txt",sep=""))
GOTermsNames <-  read.delim(paste("GOTermAnalysis/",overallComparison,"/",specificComparison,"/result_GO-Individual_Genes-geneList-",specificComparison,"_",FCweFilteredBy,"FC.csv_BP_Fisher+KSclassic_test.txt",sep=""))

#overallComparison<- "cWT_vs_vWT_cPRT6_vs_vPRT6"
#specificComparison <- "cWT_vs_vWTOnly"
#GOTermsQuantity <- read.delim(paste("GOTermAnalysis/",overallComparison,"/",specificComparison,"/result_GO-geneList-cWT_vs_vWTONLY_2FC.csv_BP_Fisher+KSclassic_test.txt",sep=""))
#GOTermsNames <-  read.delim(paste("GOTermAnalysis/",overallComparison,"/",specificComparison,"/result_GO-Individual_Genes-geneList-cWT_vs_vWTONLY_2FC.csv_BP_Fisher+KSclassic_test.txt",sep=""))
GOTerms <- data.frame(GOTermsQuantity,GOTermsNames$V3 )

#filter for only those GO terms with >10 genes associated with them
geneListAll <- GOTerms %>% filter(Significant > 10)


readCountGraphMaker <- function(geneList,GOterm, width){
  
  comp1.1genes.tmp <- comp1[ match(geneList, comp1[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[4],1,sep=""),paste(args[4],2,sep=""),paste(args[4],3,sep=""))]
  comp1.2genes.tmp <- comp1[ match(geneList, comp1[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[5],1,sep=""),paste(args[5],2,sep=""),paste(args[5],3,sep=""))]
  
  comp2.1genes.tmp <- comp2[ match(geneList, comp2[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[6],1,sep=""),paste(args[6],2,sep=""),paste(args[6],3,sep=""))]
  comp2.2genes.tmp <- comp2[ match(geneList, comp2[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[7],1,sep=""),paste(args[7],2,sep=""),paste(args[7],3,sep=""))]
  
  rowmeans<-function(df) {
    for ( i in 1:nrow(df) ) {
      df$mean[i] <- mean( as.numeric(df[i,c(3,4,5)]) ) 
      df$sd[i] <- sd( as.numeric(df[i,c(3,4,5)]) ) 
    } 
    return(df)
  }
  
  comp1.1genes <- rowmeans(comp1.1genes.tmp)
  comp1.2genes <- rowmeans(comp1.2genes.tmp)
  
  comp2.1genes<- rowmeans(comp2.1genes.tmp)
  comp2.2genes<- rowmeans(comp2.2genes.tmp)
  
  
  #if the tair_symbol is empty for whatever reason, we'll use it's ensembl_gene_id to fill in
  comp1.1genes$tair_symbol[which(comp1.1genes$tair_symbol %in% "")] <- comp1.1genes$ensembl_gene_id[which(comp1.1genes$tair_symbol %in% "")]
  comp1.2genes$tair_symbol[which(comp1.2genes$tair_symbol %in% "")] <- comp1.2genes$ensembl_gene_id[which(comp1.2genes$tair_symbol %in% "")]
  
  comp2.1genes$tair_symbol[which(comp2.1genes$tair_symbol %in% "")] <- comp2.1genes$ensembl_gene_id[which(comp2.1genes$tair_symbol %in% "")]
  comp2.2genes$tair_symbol[which(comp2.2genes$tair_symbol %in% "")] <- comp2.2genes$ensembl_gene_id[which(comp2.2genes$tair_symbol %in% "")]
  
  plotDF <- data.frame(comp1.1genes$mean,comp1.2genes$mean, comp2.1genes$mean, comp2.2genes$mean, comp1.1genes$sd,comp1.2genes$sd, comp2.1genes$sd, comp2.2genes$sd, comp2.2genes$tair_symbol)
  colnames(plotDF) <- c(args[4],args[5],args[6],args[7],args[4],args[5],args[6],args[7],"geneName")
  df1 <- melt(plotDF[,c(1:4,9)], id.vars='geneName')
  df2 <- melt(plotDF[,c(5:8,9)], id.vars='geneName')
  df1[,"sd"] <- df2$value
  
  #bar width and error bar position dodge must be the same
  g <- ggplot(data = df1, mapping = aes(x = geneName, y = value,fill = variable ,width=0.8)) + 
    geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
    ggtitle(paste(GOterm,"- genes mean Norm Expression data, +/- sd, from genes at ",FCweFilteredBy,"FC filter",sep="")) +
    theme(axis.text.x=element_text(angle = -70, hjust = 0, size=10)) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd,x = geneName), col="gray48", width=0.6,
                  position=position_dodge(0.8)) +
    labs(y= "mean normalised read count")
 plot(g)
 #need to check if the GOterm has a "/" in it, otherwise it won;t print into file correctly, it will think there's a directory coming up.
 if(grepl("/", GOterm, fixed = TRUE)){
   GOterm <- str_remove(GOterm, "/")
 }
 ggsave(paste("/Users/joannachustecki/Documents/RNAseqAnalysis/GOTermAnalysis/",overallComparison,"/",specificComparison,"/",GOterm,"byReadCount_",specificComparison,".pdf",sep=""),g, width=width1,height=5)
}

for( g in 1:nrow(geneListAll)){
  geneList<- unlist(str_split(geneListAll[g,9],", "))
  GOterm <- geneListAll[g,2]
  
  #This next loop is just for plotting purporses.
  if (length(geneList) > 100 ){
    width1 = 45
  } else {
    width1 = 15
  }
  
  graph <- readCountGraphMaker(geneList, GOterm, width1)
  
}

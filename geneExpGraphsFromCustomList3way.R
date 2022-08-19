#The same script as geneExpGraphsFromGoTerms, but adapted for three comparisons, instead of the regular 2.

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


#Script for making gene expression graphs from gene lists
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)

# currently coming in as input

# args[1] <- "VRN2_RNA_seq_data_Untreated_vs_Treated/cWT_vs_cPRT6_0FC.csv" 
# args[2] <- "VRN2_RNA_seq_data_Untreated_vs_Treated/cWT_vs_vWT_0FC.csv" 
# args[3] <- "VRN2_RNA_seq_data_Untreated_vs_Treated/cDM_vs_vDM_0FC.csv" 
# args[4] <- "cWT_vs_vWT_cPRT6_vs_vPRT6_cDM_vs_vDM"
# args[5] <-"WT.C"
# args[6] <-"WT.V"
# args[7] <-"prt6.C"
# args[8] <-"prt6.V"
# args[9] <-"DM.C"
# args[10] <-"DM.V"
# args[11] <- "cWT_vs_vWTunique_vs_cPRT6_vs_vPRT6andcDM_vs_vDM"


#setwd("/Users/joannachustecki/Documents/RNAseqAnalysis/")
#CWTvsVWT <- read.csv("VRN2_RNA_seq_data_Untreated_vs_Treated/CWTvsVWT_0FC.csv")
#Cprt6vsVprt6 <- read.csv("VRN2_RNA_seq_data_Untreated_vs_Treated/Cprt6vsVprt6_0FC.csv")

#comp1.name <- "VRN2_RNA_seq_data_Untreated_vs_Treated/cPRT6vsvPRT6_0FC.csv"
#comp2.name <- "VRN2_RNA_seq_data_Untreated_vs_Treated/CDMvsVDM_0FC.csv"

comp1.name <- args[1]
comp2.name <- args[2]
comp3.name <- args[3]

comp1 <- read.csv(paste(getwd(),comp1.name,sep="/"))
comp2 <- read.csv(paste(getwd(),comp2.name,sep="/"))
comp3 <- read.csv(paste(getwd(),comp3.name,sep="/"))


#with our particular data set, the WT comparison data frame included two extra FRI.C1/V1 samples, affecting the normalisation. All 
# of our other data frames have the same values across the dataframe, so it won't matter, just as long as we're not using the WT frame for gathering read count data.
#So I've introduced a specific catch for that data frame
if(comp1.name == "VRN2_RNA_seq_data_Untreated_vs_Treated/cWTvsvWT_0FC.csv"){
  comp1 <- comp2
  print("CWTvsVWT in comparison")
}
if(comp2.name == "VRN2_RNA_seq_data_Untreated_vs_Treated/cWTvsvWT_0FC.csv"){
  comp2 <- comp1
  print("CWTvsVWT in comparison")
}
if(comp3.name == "VRN2_RNA_seq_data_Untreated_vs_Treated/cWTvsvWT_0FC.csv"){
  comp3 <- comp1
  print("CWTvsVWT in comparison")
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#For making gene expression graphs from GO terms

overallComparison<- args[4]
specificComparison <- args[11]

FCweFilteredBy<- as.numeric(args[12])

#overallComparison<- "cWT_vs_vWT_cPRT6_vs_vPRT6_cDM_vs_vDM" 
#specificComparison <- "intersectonlycWT_vs_vWT_cPRT6_vs_vPRT6"

#geneList-intersectonlycWT_vs_vWT_cPRT6_vs_vPRT6_2FC

geneList_ <- read.csv(paste("GOTermAnalysis/",overallComparison,"/",specificComparison,"/geneList-",specificComparison,"_",FCweFilteredBy,"FC.csv",sep=""))
geneListAll <- unlist(list(geneList_$x))

readCountGraphMaker <- function(geneList, width){
  
  comp1.1genes.tmp <- comp1[ match(geneList, comp1[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[5],1,sep=""),paste(args[5],2,sep=""),paste(args[5],3,sep=""))]
  comp1.2genes.tmp <- comp1[ match(geneList, comp1[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[6],1,sep=""),paste(args[6],2,sep=""),paste(args[6],3,sep=""))]
  
  comp2.1genes.tmp <- comp2[ match(geneList, comp2[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[7],1,sep=""),paste(args[7],2,sep=""),paste(args[7],3,sep=""))]
  comp2.2genes.tmp <- comp2[ match(geneList, comp2[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[8],1,sep=""),paste(args[8],2,sep=""),paste(args[8],3,sep=""))]
  
  comp3.1genes.tmp <- comp3[ match(geneList, comp3[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[9],1,sep=""),paste(args[9],2,sep=""),paste(args[9],3,sep=""))]
  comp3.2genes.tmp <- comp3[ match(geneList, comp3[,2]) , c("ensembl_gene_id","tair_symbol",paste(args[10],1,sep=""),paste(args[10],2,sep=""),paste(args[10],3,sep=""))]
  
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
  
  comp3.1genes<- rowmeans(comp3.1genes.tmp)
  comp3.2genes<- rowmeans(comp3.2genes.tmp)
  
  
  #if the tair_symbol is empty for whatever reason, we'll use it's ensembl_gene_id to fill in
  comp1.1genes$tair_symbol[which(comp1.1genes$tair_symbol %in% "")] <- comp1.1genes$ensembl_gene_id[which(comp1.1genes$tair_symbol %in% "")]
  comp1.2genes$tair_symbol[which(comp1.2genes$tair_symbol %in% "")] <- comp1.2genes$ensembl_gene_id[which(comp1.2genes$tair_symbol %in% "")]
  
  comp2.1genes$tair_symbol[which(comp2.1genes$tair_symbol %in% "")] <- comp2.1genes$ensembl_gene_id[which(comp2.1genes$tair_symbol %in% "")]
  comp2.2genes$tair_symbol[which(comp2.2genes$tair_symbol %in% "")] <- comp2.2genes$ensembl_gene_id[which(comp2.2genes$tair_symbol %in% "")]
  
  comp3.1genes$tair_symbol[which(comp3.1genes$tair_symbol %in% "")] <- comp3.1genes$ensembl_gene_id[which(comp3.1genes$tair_symbol %in% "")]
  comp3.2genes$tair_symbol[which(comp3.2genes$tair_symbol %in% "")] <- comp3.2genes$ensembl_gene_id[which(comp3.2genes$tair_symbol %in% "")]
  
  plotDF <- data.frame(comp1.1genes$mean, comp1.2genes$mean, 
                       comp2.1genes$mean, comp2.2genes$mean, 
                       comp3.1genes$mean, comp3.2genes$mean, 
                       comp1.1genes$sd, comp1.2genes$sd, 
                       comp2.1genes$sd, comp2.2genes$sd, 
                       comp3.1genes$sd, comp3.2genes$sd, 
                       comp2.2genes$tair_symbol)
  colnames(plotDF) <- c(args[5],args[6],args[7],args[8],args[9],args[10],args[5],args[6],args[7],args[8],args[9],args[10],"geneName")
  df1 <- melt(plotDF[,c(1:6,13)], id.vars='geneName')
  df2 <- melt(plotDF[,c(7:12,13)], id.vars='geneName')
  df1[,"sd"] <- df2$value
  
  #bar width and error bar position dodge must be the same
  g <- ggplot(data = df1, mapping = aes(x = geneName, y = value,fill = variable ,width=0.8)) + 
    geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
    ggtitle(paste(length(geneListAll)," ",specificComparison," genes mean Norm Expression data, +/- sd, from genes at ",FCweFilteredBy,"FC filter",sep="")) +
    theme(axis.text.x=element_text(angle = -70, hjust = 0, size=10)) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd,x = geneName), col="gray48", width=0.6,
                  position=position_dodge(0.8)) +
    labs(y= "mean normalised read count")
  plot(g)

  ggsave(paste("/Users/joannachustecki/Documents/RNAseqAnalysis/GOTermAnalysis/",overallComparison,"/",specificComparison,"/ReadCount_",specificComparison,"_2FC.pdf",sep=""),g, width=width1,height=5)
}


  #This next loop is just for plotting purporses.
  if (length(geneListAll) > 100 ){
    width1 = 45
  } else {
    width1 = 15
  }
  
  graph <- readCountGraphMaker(geneListAll, width1)

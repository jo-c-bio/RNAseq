#!/usr/bin/env Rscript

#Load input file

#fileType <- "CWTvsVWT"
#DEgenesList <- read.csv("/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_Untreated_vs_Treated/CWTvsVWT_0FC.csv")
args = commandArgs(trailingOnly=TRUE)

fileType <- args[1]
DEgenesList <- read.csv(args[2])

# we can't use gsub to Remove '#' and '/'  to give nas, as it changes whole column to a character field, so it can't find the appropriate values anymore,

#What fold change are we looking for?
#If a gene has 1.5x more reads in the treated WT, it will have a fold change of -0.5849625
#If a gene has 2x more reads in the treated WT, it will have a fold change of -1
#If a gene has 1.5x less reads in the treated WT, it will have a fold change of 0.5849625
#If a gene has 2x less reads in the treated WT, it will have a fold change of 1

#For each set of genes, we want only those that pass (are lower than) an adjusted significance threshold of 0.5
#Let's pick out genes where the treated sample has more than 1.5x as many reads than control
          
#all of these filters go in 2 steps as the adj values need to be ran as.numeric(), but without the NA values
UPby1_5initial <- DEgenesList[DEgenesList$log2FoldChange <= -0.5849625  & DEgenesList$padj != "#N/A" , ]
UPby1_5 <- UPby1_5initial[as.numeric(UPby1_5initial$padj) <= 0.05 , ]
write.csv(UPby1_5,paste(fileType,"UPby1_5.csv",sep="_"))

#now let's pick out genes that have at least 2x more reads in the treated compared to the control

UPby2initial <- DEgenesList[DEgenesList$log2FoldChange < -1 & DEgenesList$padj != "#N/A" , ]
UPby2 <- UPby2initial[as.numeric(UPby2initial$padj) <= 0.05 , ]
write.csv(UPby2,paste(fileType,"UPby2.csv",sep="_"))

#and 3x 

UPby3initial <- DEgenesList[DEgenesList$log2FoldChange < -1.584963 & DEgenesList$padj != "#N/A" , ]
UPby3 <- UPby3initial[as.numeric(UPby3initial$padj) <= 0.05 , ]
write.csv(UPby3,paste(fileType,"UPby3.csv",sep="_"))

#and 4x

UPby4initial <- DEgenesList[DEgenesList$log2FoldChange < -2  & DEgenesList$padj != "#N/A" , ]
UPby4 <- UPby4initial[as.numeric(UPby4initial$padj) <= 0.05 , ]
write.csv(UPby4,paste(fileType,"UPby4.csv",sep="_"))

#Now, which genes are downregulated
#down by 1.5x

DOWNby1_5initial <- DEgenesList[DEgenesList$log2FoldChange >= 0.58 & DEgenesList$padj != "#N/A" , ]
DOWNby1_5 <- DOWNby1_5initial[as.numeric(DOWNby1_5initial$padj) <= 0.05 , ]
write.csv(DOWNby1_5,paste(fileType,"DOWNby1_5.csv",sep="_"))

#down by 2x

DOWNby2initial <- DEgenesList[DEgenesList$log2FoldChange >= 1 & DEgenesList$padj != "#N/A"  , ]
DOWNby2 <- DOWNby2initial[as.numeric(DOWNby2initial$padj) <= 0.05 , ]
write.csv(DOWNby2,paste(fileType,"DOWNby2.csv",sep="_"))

#down by 3x

DOWNby3initial <- DEgenesList[DEgenesList$log2FoldChange >= 1.58 & DEgenesList$padj != "#N/A" , ]
DOWNby3 <- DOWNby3initial[as.numeric(DOWNby3initial$padj) <= 0.05 , ]
write.csv(DOWNby3,paste(fileType,"DOWNby3.csv",sep="_"))

#down by 4x

DOWNby4initial <- DEgenesList[DEgenesList$log2FoldChange >=  2 & DEgenesList$padj != "#N/A" , ]
DOWNby4 <- DOWNby4initial[as.numeric(DOWNby4initial$padj) <= 0.05 , ]
write.csv(DOWNby4,paste(fileType,"DOWNby4.csv",sep="_"))





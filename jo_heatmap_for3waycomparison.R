#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggpubr)
library(ggvenn)
library(pheatmap)


dataPrep<-function(comp1, comp2, comp3){
  #checking differences-the cWTvsvWT has different amount of genes to the other comparisons.
  setdiff(comp1$ensembl_gene_id,comp2$ensembl_gene_id)
  setdiff(comp2$ensembl_gene_id,comp1$ensembl_gene_id)
  setdiff(comp2$ensembl_gene_id,comp3$ensembl_gene_id)
  
  #so, we'll add all those genes into the cWTvsvWT list with a fold change of 0.
  #Also, cWTvsvWT has one gene, AT5G06930, that all the other comparisons don't have, so we'll add that too to them, with a fold change of 0 
  
  fillerData <- function(comparisonDataFrame,inputDataFrame){
    fillerDataAdd<-data.frame(setdiff(comparisonDataFrame$ensembl_gene_id,inputDataFrame$ensembl_gene_id),
                              rep(0,length(setdiff(comparisonDataFrame$ensembl_gene_id,inputDataFrame$ensembl_gene_id))),
                              rep(0,length(setdiff(comparisonDataFrame$ensembl_gene_id,inputDataFrame$ensembl_gene_id))))
    colnames(fillerDataAdd) <-  c("ensembl_gene_id","padj","log2FoldChange")
    inputDataFrame <- rbind(inputDataFrame[ , c("ensembl_gene_id","padj","log2FoldChange")],fillerDataAdd)
    return(inputDataFrame)
  }
  
  comp1<-fillerData(comp2,comp1)
  comp2<-fillerData(comp1,comp2) 
  comp3<-fillerData(comp2,comp3) 
  
  #order by gene number
  comp1 <- comp1[ order(comp1$ensembl_gene_id), ]
  comp2 <- comp2[ order(comp2$ensembl_gene_id), ]
  comp3 <- comp3[ order(comp3$ensembl_gene_id), ]
  
  #change column names 
  colnames(comp1) <-  lapply(c("ensembl_gene_id","padj","log2"), function(x) paste("comp1", x, sep="_"))
  colnames(comp2) <-  lapply(c("ensembl_gene_id","padj","log2"), function(x) paste("comp2", x, sep="_"))
  colnames(comp3) <-  lapply(c("ensembl_gene_id","padj","log2"), function(x) paste("comp3", x, sep="_"))
  
  #make into one matrix
  matrix1<-bind_cols(comp1,comp2,comp3)
  rownames(matrix1) <- matrix1$comp1_ensembl_gene_id
  matrix <- within(matrix1, rm("comp1_ensembl_gene_id","comp2_ensembl_gene_id","comp3_ensembl_gene_id"))
  
  #This does not get rid of gene repeats, but rather genes that have the exact same expression set as another. 
  #These can be splice variants, or weirdly, there's a mito ribo gene that had the same values as a nuclear ribo gene.
  #matrix3 <- distinct(matrix2)
  
  return(matrix)
}

### data input

#What fold change do you want? ie log2(4) = 2
f <- as.numeric(args[1])
#f <- 2

comparison1<-args[2]
comparison2<-args[3]
comparison3<-args[4]
#comparison1<-"cWT_vs_vWT"
#comparison2<-"cPRT6_vs_vPRT6"
#comparison3<-"cDM_vs_vDM"

comp1 <- read.csv(args[5])
comp2 <- read.csv(args[6])
comp3 <- read.csv(args[7])
#comp1 <- read.csv('/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_Untreated_vs_Treated/cWT_vs_vWT_0FC.csv')
#comp2 <- read.csv("/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_Untreated_vs_Treated/cPRT6_vs_vPRT6_0FC.csv")
#comp3 <- read.csv("/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_Untreated_vs_Treated/cDM_vs_vDM_0FC.csv")


currentWD <- getwd()


#####venn diagram

y <- list(comp1x = comp1 %>%
            filter(padj != "#N/A" ) %>% 
            filter(as.numeric(padj) < 0.05) %>% 
            filter(log2FoldChange > f |  log2FoldChange < -f),
          comp2x = comp2 %>%
            filter( padj != "#N/A" ) %>% 
            filter(as.numeric(padj) < 0.05) %>% 
            filter(log2FoldChange > f |  log2FoldChange < -f),
          comp3x = comp3 %>%
            filter( padj != "#N/A" ) %>% 
            filter(as.numeric(padj) < 0.05) %>% 
            filter(log2FoldChange > f |  log2FoldChange < -f))
# comp4x = comp4x %>%
#   filter(padj < 0.05, padj != "#N/A") %>% 
#   filter(log2FoldChange > f |  log2FoldChange < -f)


y2 <- list(unlist(y$comp1x[,1]), unlist(y$comp2x[,1]), unlist(y$comp3x[,1]))
#,
#cPRT6_vs_vPRT6 = unlist(y$comp3x[,1]),cDM_vs_vDM = unlist(y$comp4x[,1]))
names(y2) <- c(comparison1, comparison2, comparison3)


venn <- ggvenn(
  y2, 
  fill_color = c("plum1", "#EFC000FF","Turquoise"),
  stroke_size = 0.5, set_name_size = 4, text_size = 4) 
#ggtitle(paste("Microarray (hypoxia) and RNAseq (vernalisation) data,","\n", "filtered for log2 foldchange +/- ",f,"",sep=""))
#ggtitle(paste("RNAseq (vernalisation) data,","\n", "filtered for log2 foldchange +/- ",f,"",sep=""))


venn
ggsave(paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison1,"_",comparison2,"_",comparison3,"_FC",f,".pdf",sep=""),venn)

#We want to get these gene lists of overlap between the different comparisons, so use setdiff and intersect. 
#All DEGs from comp1, regardless of comp2 or comp3.(i.e. all WT vernalised RNAseq genes as a control.)
comp1only <- y2[[1]]
#All DEGs from comp2, regardless of comp1 or comp3.(i.e. all prt6 genes changed, just incase I need it.)
comp2only <- y2[[2]]
#All DEGs from comp3, regardless of comp1 or comp2.(i.e. all prt6 hypoxia microarray genes changed, just incase I need it.)
comp3only <- y2[[3]]
#unique to comparison 1 i.e. not intersecting anything else
comp1unique <- setdiff(setdiff(y2[[1]],y2[[2]]),y2[[3]])
#unique to comparison 2 i.e. not intersecting anything else
comp2unique <- setdiff(setdiff(y2[[2]],y2[[1]]),y2[[3]])
#unique to comparison 3 i.e. not intersecting anything else
comp3unique <- setdiff(setdiff(y2[[3]],y2[[2]]),y2[[1]])
#intersection
threeways_intersect <- intersect(intersect(y2[[1]],y2[[2]]),y2[[3]])
#Because these intersections all have the three-way intersect in common, we will remove it to be left with those that only overlap between the specific pairs.
c1c2_intersect<-intersect(y2[[1]],y2[[2]])
c1c2_intersect = c1c2_intersect[ - which(c1c2_intersect %in% threeways_intersect)] 
c2c3_intersect<-intersect(y2[[2]],y2[[3]])
c2c3_intersect = c2c3_intersect[ - which(c2c3_intersect %in% threeways_intersect)] 
c3c1_intersect<-intersect(y2[[3]],y2[[1]])
c3c1_intersect = c3c1_intersect[ - which(c3c1_intersect %in% threeways_intersect)] 

#need to be careful if i'm ever doing heatmaps with this data, that the fold change 
#is currently hypoxia vs normal. So if wanted normal vs hypoxia, would just need to switch all
#the fold change values over to -ve or +ve


#Take these forward for GO term analysis
write.csv(comp1only,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison1,"Only/geneList-",comparison1,"Only_",f,"FC.csv",sep=""))
write.csv(comp2only,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison2,"Only/geneList-",comparison2,"Only_",f,"FC.csv",sep=""))
write.csv(comp3only,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison3,"Only/geneList-",comparison3,"Only_",f,"FC.csv",sep=""))

write.csv(comp1unique,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison1,"unique_vs_",comparison2,"and",comparison3,"/geneList-",comparison1,"unique_vs_",comparison2,"and",comparison3,"_",f,"FC.csv",sep=""))
write.csv(comp2unique,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison2,"unique_vs_",comparison1,"and",comparison3,"/geneList-",comparison2,"unique_vs_",comparison1,"and",comparison3,"_",f,"FC.csv",sep=""))
write.csv(comp3unique,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",comparison3,"unique_vs_",comparison1,"and",comparison2,"/geneList-",comparison3,"unique_vs_",comparison1,"and",comparison2,"_",f,"FC.csv",sep=""))

write.csv(threeways_intersect,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/ALLintersect",comparison1,"_",comparison2,"_",comparison3,"/geneList-ALLintersect",comparison1,"_",comparison2,"_",comparison3,"_",f,"FC.csv",sep=""))

write.csv(c1c2_intersect,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/intersectonly",comparison1,"_",comparison2,"/geneList-intersectonly",comparison1,"_",comparison2,"_",f,"FC.csv",sep=""))
write.csv(c2c3_intersect,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/intersectonly",comparison2,"_",comparison3,"/geneList-intersectonly",comparison2,"_",comparison3,"_",f,"FC.csv",sep=""))
write.csv(c3c1_intersect,paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/intersectonly",comparison3,"_",comparison1,"/geneList-intersectonly",comparison3,"_",comparison1,"_",f,"FC.csv",sep=""))

#heatmap

matrix<-dataPrep(comp1,comp2,comp3)

#we can generate the heatmap to be only intersecting genes in the lists, or we can make them so all the genes in the lists 
#are in the heatmap, and then if they're significant we can see the fold change. 

matrixIntersect <- matrix %>%
  #%>% # filter for only significant DEGs, and ones where the pvalue isnt n/a.
  #This little script give the intersecting list between the two comparisons, as you're asking it to filter comp1, then use that to filter for comp2. 
  filter(comp2_padj != "#N/A", comp1_padj != "#N/A" , comp3_padj != "#N/A" ) %>%
  filter(as.numeric(comp1_padj) < 0.05, as.numeric(comp2_padj) < 0.05, as.numeric(comp3_padj) < 0.05) %>%
  filter(comp1_log2 > f |  comp1_log2 < -f)%>%
  filter(comp2_log2 > f |  comp2_log2 < -f)%>%
  filter(comp3_log2 > f |  comp3_log2 < -f)

conditions<- c("comp1_log2","comp2_log2","comp3_log2")
matrixIntersect <- as.matrix(matrixIntersect[conditions])

#OR we can generate it to be full lists of the genes in the untreated/treated and see what changes. (we'll still filter for significance)
matrixFull  <- matrix %>%
  #%>% # filter for only significant DEGs, and ones where the pvalue isnt n/a.
  #This little script give the intersecting list between the two comparisons, as you're asking it to filter comp1, then use that to filter for comp2. 
  filter(comp2_padj != "#N/A", comp1_padj != "#N/A", comp3_padj != "#N/A") %>%
  filter(as.numeric(comp1_padj) < 0.05, as.numeric(comp2_padj) < 0.05, as.numeric(comp3_padj) < 0.05 )
conditions<- c("comp1_log2","comp2_log2","comp3_log2")
matrixFull <- as.matrix(matrixFull[conditions])

#OR we can generate it to be all the significant genes from comp1, or all the genes from comp2. ie, they will get printed if only one data point of the two comps crossed the sig. threshold.
matrixTest  <- matrix %>%
  filter( comp2_padj != "#N/A", comp1_padj != "#N/A" ) %>%
  filter(as.numeric(comp1_padj) < 0.05| as.numeric(comp2_padj) < 0.05)
conditions<- c("comp1_log2","comp2_log2")
matrixTest <- as.matrix(matrixTest[conditions])

#OR we can generate it to be all the fc genes from comp1, or all the fc genes from comp2. ie, they will get printed if only one data point of the two comps crossed the fc threshold.
matrixTest2  <- matrix %>%
  filter(comp2_padj != "#N/A", comp1_padj != "#N/A") %>%
  filter(as.numeric(comp1_padj) < 0.05, as.numeric(comp2_padj) < 0.05 ) %>%
  filter((comp1_log2 > f |  comp1_log2 < -f) | (comp2_log2 > f |  comp2_log2 < -f))
conditions<- c("comp1_log2","comp2_log2")
matrixTest2 <- as.matrix(matrixTest2[conditions])

#or we can just do absolutely everything - this is meaningless, it's 18000 genes.
conditions<- c("comp1_log2","comp2_log2")
matrix1 <- as.matrix(matrix[conditions])

#here we'll just take those genes for the unique sides of the venn diagram. The unique lists have already been filtered for
#significance and fold change when making the venn diagram.
matrix_comp1_uniq <- matrix[match(comp1unique,row.names(matrix)), ]
conditions<- c("comp1_log2","comp2_log2")
matrix_comp1_unique <- as.matrix(matrix_comp1_uniq[conditions])

matrix_comp2_uniq <- matrix[match(comp2unique,row.names(matrix)), ]
conditions<- c("comp1_log2","comp2_log2")
matrix_comp2_unique <- as.matrix(matrix_comp2_uniq[conditions])

### To generate Heat map
generateHM <- function(matrix, title, vennOption,outFileName){
  
  pHM<-pheatmap(matrix, show_rownames = T, 
                #color = cols,
                fontsize_row = 0.5,
                scale = "none", border_color = NA,
                fontsize_col = 10,treeheight_col = 0,
                labels_col= c(comparison1,comparison2,comparison3),height=30,width=10, main = title)
  
  if( vennOption == "YesVenn"){
    vennpHM <- ggarrange(plotlist= list(venn,pHM$gtable), ncol= 1, nrow=2, heights = c(1, 2))
    vennpHM <- annotate_figure(vennpHM, top = text_grob(paste("log2 fold change = ",f,sep=""), 
                                                        color = "black", face = "bold", size = 16) )}
  else {
    vennpHM <- pHM
  }
  print(vennpHM)
  ggsave(paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/fc",f,outFileName,".pdf",sep=""),vennpHM, height=12,width=6)
  write.csv(rownames(matrix),paste(currentWD,"/GOTermAnalysis/",comparison1,"_",comparison2,"_",comparison3,"/",f,"FC",outFileName,".csv",sep=""))
  
}
#generateHM(matrixTest2, paste("Only one side needs to cross fc",f,"and all significant"), "NoVenn")
#generateHM(matrixTest, "Only one side needs to be significant", "NoVenn")
generateHM(matrixFull, paste("All genes of all comparisons pass significance threshold,","\n","no fold change filter"), "NoVenn",paste(comparison1,comparison2,comparison3,"FullMatrix",sep=""))
generateHM(matrixIntersect, "All genes filtered for significance and fold change", "YesVenn",paste(comparison1,comparison2,comparison3,"intersect",sep=""))
#generateHM(matrix_comp1_unique,paste("The",nrow(matrix_comp1_unique),"genes unique to",comparison1,"comparision",sep=" "), "NoVenn", paste(comparison1,"unique",sep=""))
#generateHM(matrix_comp2_unique,paste("The",nrow(matrix_comp2_unique),"genes unique to",comparison2,"comparision",sep=" "), "NoVenn",paste(comparison2,"unique",sep=""))




#jo and derrys test for specifying heatmap colours. 
#needs  color = cols, putting back into the heatmap function above.

# m <- matrix(c(rnorm(1000)), ncol=100)
# distmat <- dist(t(m))
# 
# 
# makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
# {
#   stopifnot(length(colors) == 4)
#   ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
#   ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
#   return(c(ramp1, ramp2))
# }
# 
# cutoff.distance <- 3  
# cols <- makeColorRampPalette(c("white", "red",    # distances 0 to 3 colored from white to red
#                                "green", "black"), # distances 3 to max(distmat) colored from green to black
#                              cutoff.distance / max(distmat),
#                              100)

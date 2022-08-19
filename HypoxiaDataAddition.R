#Figuring out how to make three way venn diagrams

#libraries
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggvenn)

### data input

#What fold change do you want? ie log2(4) = 2
#f <- as.numeric(args[1])
f <- 1.5

#comparison1<-args[2]
#comparison2<-args[3]
comparison1<-"cWT_vs_cPRT6"
comparison2<-"cWT_vs_vWT"
comparison3<-"WT-HS_vs_WT-NS"
comparison4<-"prt6HS_vs_prt6NS"


#comp1 <- read.csv(args[4])
#comp2 <- read.csv(args[5])
comp1 <- read.csv("/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_WT_mutants_control/CWTvsCprt6_0FC.csv")
comp2 <- read.csv('/Users/joannachustecki/Documents/RNAseqAnalysis/VRN2_RNA_seq_data_Untreated_vs_Treated/CWTvsVWT_0FC.csv')
comp3 <- read.csv("/Users/joannachustecki/Documents/RNAseqAnalysis/withHypoxiaData/Homeostatic_response_paper_ESM/TranscriptomeData-Table1.csv",header=TRUE)

currentWD <- getwd()

#we'll split the hypoxia data frame into the comparisons we want, that is WTHS_vs_WTNS and prt6HS_vs_prt6NS
comp3.a.WTHS_vs_WTNS <- data.frame(comp3$locus,comp3$WTHSvsWTNS_Log2signal)
comp3.b.prt6HS_vs_prt6NS <- data.frame(comp3$locus,comp3$prt6HSvsprt6NS_Log2signal)

colnames(comp3.a.WTHS_vs_WTNS) <- c("geneName","log2Signal")
colnames(comp3.b.prt6HS_vs_prt6NS) <- c("geneName","log2Signal")

#How to better label these for the future? Might just need this to be a specific script for the hypoxia data.

#####venn diagram

y <- list(comp1x = comp1 %>%
            filter(padj != "#N/A" ) %>% 
            filter(as.numeric(padj) < 0.05) %>% 
            filter(log2FoldChange > f |  log2FoldChange < -f),
          comp2x = comp2 %>%
            filter( padj != "#N/A" ) %>% 
            filter(as.numeric(padj) < 0.05) %>% 
            filter(log2FoldChange > f |  log2FoldChange < -f),
          comp3.a = comp3.a.WTHS_vs_WTNS %>%
            filter(log2Signal > f |  log2Signal < -f),
          comp3.b = comp3.b.prt6HS_vs_prt6NS %>%
            filter(log2Signal > f |  log2Signal < -f))
# comp4x = comp4x %>%
#   filter(padj < 0.05, padj != "#N/A") %>% 
#   filter(log2FoldChange > f |  log2FoldChange < -f)


#Just want the cWT and hWT for now, not hPRT6.
y2 <- list(unlist(y$comp1x[,1]), unlist(y$comp2x[,1]),unlist(y$comp3.a[,1]))
#,
#cPRT6_vs_vPRT6 = unlist(y$comp3x[,1]),cDM_vs_vDM = unlist(y$comp4x[,1]))
names(y2) <- c(comparison1, comparison2,comparison3)

venn <- ggvenn(
  y2, 
  fill_color = c("plum1", "#EFC000FF","Turquoise"),
  stroke_size = 0.5, set_name_size = 4, text_size = 4) +
  ggtitle(paste("Microarray (hypoxia) and RNAseq (vernalisation) data,","\n", "filtered for log2 foldchange +/- ",f,"",sep=""))
venn
ggsave(paste(currentWD,"/withHypoxiaData/",comparison1,"_",comparison2,"_",comparison3,"_FC",f,".pdf",sep=""),venn)


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
c1c2_intersect<-intersect(y2[[1]],y2[[2]])
c2c3_intersect<-intersect(y2[[2]],y2[[3]])
c3c1_intersect<-intersect(y2[[3]],y2[[1]])

#need to be careful if i'm ever doing heatmaps with this data, that the fold change 
#is currently hypoxia vs normal. So if wanted normal vs hypoxia, would just need to switch all
#the fold change values over to -ve or +ve


#Take these forward for GO term analysis
write.csv(comp1only,paste(currentWD,"/withHypoxiaData/geneList-",comparison1,"Only_",f,"FC.csv",sep=""))
write.csv(comp2only,paste(currentWD,"/withHypoxiaData/geneList-",comparison2,"Only_",f,"FC.csv",sep=""))
write.csv(comp3only,paste(currentWD,"/withHypoxiaData/geneList-",comparison3,"Only_",f,"FC.csv",sep=""))

write.csv(comp1unique,paste(currentWD,"/withHypoxiaData/geneList-",comparison1,"unique_vs_",comparison2,"and",comparison3,"_",f,"FC.csv",sep=""))
write.csv(comp2unique,paste(currentWD,"/withHypoxiaData/geneList-",comparison2,"unique_vs_",comparison1,"and",comparison3,"_",f,"FC.csv",sep=""))
write.csv(comp3unique,paste(currentWD,"/withHypoxiaData/geneList-",comparison3,"unique_vs_",comparison1,"and",comparison2,"_",f,"FC.csv",sep=""))

write.csv(threeways_intersect,paste(currentWD,"/withHypoxiaData/geneList-ALLintersect",comparison1,"_",comparison2,"_",comparison3,"_",f,"FC.csv",sep=""))

write.csv(c1c2_intersect,paste(currentWD,"/withHypoxiaData/geneList-intersectonly",comparison1,"_",comparison2,"_",f,"FC.csv",sep=""))
write.csv(c2c3_intersect,paste(currentWD,"/withHypoxiaData/geneList-intersectonly",comparison2,"_",comparison3,"_",f,"FC.csv",sep=""))
write.csv(c3c1_intersect,paste(currentWD,"/withHypoxiaData/geneList-intersectonly",comparison3,"_",comparison1,"_",f,"FC.csv",sep=""))


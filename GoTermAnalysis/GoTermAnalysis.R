#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


library(stringr)
library(reshape2)
library(ggplot2)

l <- unlist(str_split(args[1],"/"))
#Name  <- "geneList-cPRT6_vs_vPRT6unique_vs_cWT_vs_vWT_2FC.csv"
Name <- l[4]
#Name <-"cWTvWTcustom_4"

#dirName<-"cPRT6_vs_vPRT6unique_vs_cWT_vs_vWT"
dirName <- l[3]
#dirName<-"heatMapRegionalSplits"

baseWD<- getwd()

setwd(paste(l[1],l[2],sep="/"))
#setwd("GOTermAnalysis")


#setwd(paste('/Users/joannachustecki/Documents/RNAseqAnalysis/GOTermAnalysis',dirName,sep=""))
#DEGs <- read.csv('geneList-cPRT6_vs_vPRT6intersectcWT_vs_vWT_1.5FC.csv', header = TRUE) ##### This is your list/table of DEGs
DEGs <- read.csv(paste(dirName,"/",Name,sep=""), header = FALSE) ##### This is your list/table of DEGs


DEGenes <- unlist(list(DEGs))
#DEGenes <- DEGs[,2]
DEGenes <- unique(DEGenes)
#directory<-paste("/Users/lfulfy/Documents/MIBTP/PhD Year 1/Rory/GO_Term_Analysis/",Name,sep="")
#dir.create(directory); savedir<-(directory); setwd(savedir)
GenID_coef <- read.csv(paste(baseWD,l[1],"All_Arabidopsis_Genes.csv",sep="/"), header = FALSE) #universe of genes for enrichment analysis

GenID_coef <- unique(GenID_coef)

BiocManager::install('topGO')
BiocManager::install('org.At.tair.db')

### End Rory input data

# Original algorithm for Anna, with directory paths commented out to use for Rory
AnalyseGOTerms<-function(DEGenes,Name){
  algorithm<-"classic"
  #directory<-paste("/Users/lfulfy/Documents/MIBTP/PhD Year 1/Anna/Microarray_Analysis/Initial_Analysis/","CrazyGOTermTesting","_",Sys.Date(),"/GOTermResults","_",Name,sep="")
  #dir.create(directory); savedir<-(directory); setwd(savedir)
  for(ont in c("BP","MF","CC")){
    print (paste("testing",ont))
    ## get average expressions; source("https://bioconductor.org/biocLite.R")
    #biocLite('topGO')
    library(topGO)
    #Pre-load the Arabidopsis GO annotation
    #biocLite('org.At.tair.db')
    library(org.At.tair.db)
    #allGO2genes contain a list of GO terms Biological processes (BP). Can remake using (MF) for metabolic function etc...
    allGO2genes = annFUN.org(whichOnto=ont, feasibleGenes = NULL,
                             mapping="org.At.tair.db") #ID = SYMBOL for Anna analysis


    #Generate a "universe" of genes to use as a background for the re-sampling enrichment analysis. Here
    # the exact same gene list used for DE has been used, which is the most appropriate.
    all.genes<- as.vector(GenID_coef$V1) # WAS as.vector(droplevels(GenID_coef$GENE_SYMBOL)) for Anna analysis

    #geneList variable contains 0's and 1's against each gene where a 1 means a gene has been DE.
    geneList <- factor(as.integer(all.genes %in% DEGenes))
    #If statement handles the special case of having NO DE genes - levels is the amount of info in the returned geneList vector,
    #where only one level is present in the case of gene names with only 0 values attached (0=noDE, 1=DE)
    if(length(levels(geneList))>1){
      #Reattaches gene names to the 0/1 label of DE
      names(geneList) <- all.genes

      # first run the ‘new’ function which will create a topGO object.
      GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                    annot = annFUN.GO2genes, GO2genes = allGO2genes)
      ### Enrichment Test ###
      # Then run the  'runTest' function. Can choose many different statistical tests.  Here, algorithm = "elim", statistic = "ks".
      results.ks <- runTest(GOdata, algorithm=algorithm, statistic = "ks")
      results.fish <- runTest(GOdata, algorithm=algorithm, statistic = "fisher")
      # The GenTable function returns a data frame containing the top topNodes GO terms identified by the elim
      # algorithm, s
      ############################################################
      #ResultsTableGO <- GenTable( GOdata, ks.weight = results.ks,
      #fish.weight = results.fish,
      #orderBy = "fish.weight",topNodes = 40)
      ResultsTableGO <- GenTable( GOdata, ks.weight = results.ks,
                                  fish.weight = results.fish,
                                  orderBy = "fish.weight",topNodes = 500, numChar=50)
      ResultsTableGO<-ResultsTableGO[which(ResultsTableGO$fish.weight<0.05),]
      ###################################################################
      write.table(ResultsTableGO,file=paste(dirName,"/result_GO-",Name,"_",ont,"_","Fisher+KS",algorithm,"_test.txt",sep=""),sep="\t")

      #####################EXTRACT GENES FROM 5 MOST SIG ENRICHED LISTS########
      # retrieve genes2GO list from the "expanded" annotation in GOdata
      allGO = genesInTerm(GOdata)
      noGOTerms<-nrow(ResultsTableGO)
      SAM_ANOTATION = lapply(allGO,function(x) x[x %in% DEGenes] )
      GOCATS<-c(ResultsTableGO[1:noGOTerms,1])
      GOTERM<-c(ResultsTableGO[1:noGOTerms,2])
      GO_Gene_results<-matrix(nrow = noGOTerms,ncol = 3,data = 0)
      for (i in 1:noGOTerms){
        GO_Gene_results[i,1]<-GOCATS[i]
        GO_Gene_results[i,2]<-GOTERM[i]
        numberofgenes<-paste(SAM_ANOTATION[[GOCATS[i]]],collapse=", ")
        GO_Gene_results[i,3]<-numberofgenes
      }
      TableName<-paste(dirName,"/result_GO-Individual_Genes-",Name,"_",ont,"_","Fisher+KS",algorithm,"_test.txt",sep="")
      write.table(GO_Gene_results,file=TableName,sep="\t")
      # AND WRITE TO FOLDER
    }
  }
}



AnalyseGOTerms(DEGenes,Name) #run function below to generate AnalyseGOTerms function used here

result_GO <- read.delim(paste(dirName,"/result_GO-",Name,"_BP_","Fisher+KSclassic_test.txt",sep=""), sep="\t")


df1 <- result_GO[,c("Term","Significant","Expected")]
#MELT! Such a useful tool!!!!
df2 <- melt(df1, id.vars='Term')

#notice that the fishers test weighting is transformed by -log10- so the higher the number, the lower the pvalue = more significant.
fishweightBarplot<- ggplot(data = df2, mapping = aes(x = value, y = value)) + 
                    geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
                    coord_flip() +
                    ggtitle(Name) +
                    geom_text(aes( label=c(round(-log10(as.numeric(result_GO$fish.weight)), 3),
                                           rep(NA, length(result_GO$fish.weight))) ),
                              nudge_y = 4, size = 5 ) 
fishweightBarplot

#ggsave(paste(dirName,"/",Name,"_sigFishweightBarPlot.pdf",sep=""),fishweightBarplot,
#       height= 20, width = 10)

fishweightBarplotOrdered<- ggplot(data = df2, mapping = aes(x = reorder(Term, value), y = value)) + 
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
  coord_flip() +
  ggtitle(Name) +
  labs(x = "Term") +
  geom_text(aes( label=c(round(-log10(as.numeric(result_GO$fish.weight)), 3),
                         rep(NA, length(result_GO$fish.weight))) ),
            nudge_y = 4, size = 5 ) +
  theme(axis.text=element_text(size=10))

ggsave(paste(dirName,"/",Name,"_sigFishweightBarPlotOrdered.pdf",sep=""),fishweightBarplotOrdered,
       height= 25, width = 10)



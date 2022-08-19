library(DESeq2)

normReadCounts<-read.csv("VRN2_RNA_seq_data_AllTogether/vWT_vs_vVRN2_0FC.csv")


for(i in 16:ncol(normReadCounts)){
  normReadCounts[,i]<-as.integer(normReadCounts[,i])
}

trim<-normReadCounts[,16:ncol(normReadCounts)]

row.names(trim) <- normReadCounts[,1]

trim<- as.matrix(trim)

#Plot PCA using  DEGs
#this function applies log2 transform across the count data, to normalize it for clustering
rld <- rlog(trim)

pca<-prcomp(t(rld))

#DESeqDataSetFromMatrix(countData = trim,
                  #     colData = coldata,
                   #    design = ~ condition)

nonNumericConditons<- substr(rownames(pca$x),1,nchar(rownames(pca$x))-1)

#if you want to plot any other principle components apart from 1 and 2, just 
#specify which in both the e.g. plot(pca$x[,2],pca$x[,3],) and text() command
#default plot(pca$x will just plot PC1 and 2)
pdf("PCAplotForAllSamples.pdf", width=7.5, height=6)
plot(
  pca$x,
  col = as.factor(nonNumericConditons ),
  cex = 3.0)
text(pca$x[,"PC1"], pca$x[,"PC2"], labels=nonNumericConditons, cex = 0.4)
title(paste("PCA showing different expression profiles for untreated/treated samples","\n","(normalised read counts are log2 tranformed prior to PCA)", sep=""))

par(xpd=TRUE)
legend( 31.5,29,
  bty = "n",
  c(nonNumericConditons),
  fill = as.factor(nonNumericConditons),
  cex = 0.45)
dev.off()

#percentage of the total variation in the dataset can be explained by each principal component?
#we want the first 2 or 3 PCs in the scree plot to cover most of the variation in the data 

#We can compute the proportion of variation for each PC using std, i.e standard deviation values from the PCA result
var_explained_df <- data.frame(PC=colnames(pca$x)[1:9],
                               var_explained=(pca$sdev[1:9])^2/sum((pca$sdev[1:9])^2))
library(dplyr)
library(ggplot2)
ve<-var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained))+
  geom_col()+
  labs(title="Scree plot: PCA on scaled data (only PC1:9 plotted)")
ggsave( "here.pdf",ve, width=6,height=4)


  
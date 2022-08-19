#!/bin/bash

#This script takes two of the control vs treatment gives you lists of significant DEGs at your specified fold changes 
#that unique or shared between them and builds a venn digram and a heatmap. You also get an exported list of DE Gene names.

#Let's start by nominating our comparisons we'd like to run, we'll use these throughout the rest of the script
#To make directories, and communicate to R which files to use. 
#if you want to specify these here, use:

## comparison1="cWT_vs_cVRN2"
## comparison2="cDM_vs_cPRT6"

#else, use the wrapper format:
comparison1=${1}
comparison2=${2}


#Then later we need some column names of our dataframes, these need to be the same as in your original dataframes. HC = header comparison
#if you want to specify these here, use:
## HC1="WT.C"
## HC2="vrn2.C"
## HC3="DM.C"
## HC4="prt6.C"
#else, use the wrapper format:
HC1=${3}
HC2=${4}
HC3=${5}
HC4=${6}


#Then let's use these to make some directories for storing all the sec tions of our comaprion venn diagrams. 
#Notation, 
    #Intersect is all those genes in the middle of the venn diagram
    #Only are all those genes that are within a full circle, so all the genes associated with a comparison whether they overlap with the other comparison not.
    #Unique are the genes either side of the venn diagram, those that are for one comparison but are not shared at all with the other. 


mkdir GOTermAnalysis/${comparison1}_${comparison2}

mkdir GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}intersect${comparison1}
mkdir GOTermAnalysis/${comparison1}_${comparison2}/${comparison1}Only
mkdir GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}unique_vs_${comparison1}
mkdir GOTermAnalysis/${comparison1}_${comparison2}/${comparison1}unique_vs_${comparison2}
mkdir GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}Only


#Here we loop over 'f', our chosen fold change values we want to analyse at. I chose 1, 1.5 and 2 fold change here.
#It is within this script that we filter genes by adjusted p-values, and write the DEG lists out into .csv format. 
for f in 0.25 0.58 1 1.5 2

do
    #variables are: fold-change comparison1 comparison2 file-comparison1 file-comparison2
    #Rscript jo_heatmap_v2.r $f cPRT6_vs_vPRT6 cDM_vs_vDM VRN2_RNA_seq_data_Untreated_vs_Treated/cPRT6vsvPRT6_0FC.csv VRN2_RNA_seq_data_Untreated_vs_Treated/CDMvsVDM_0FC.csv
    Rscript jo_heatmap_v2.r $f $comparison1 $comparison2 VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv
done 

for f in 0.25 0.58 1 1.5 2

do
    Rscript GOTermAnalysis/GoTermAnalysis.R GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}intersect${comparison1}/geneList-${comparison2}intersect${comparison1}_${f}FC.csv
    Rscript GOTermAnalysis/GoTermAnalysis.R GOTermAnalysis/${comparison1}_${comparison2}/${comparison1}Only/geneList-${comparison1}Only_${f}FC.csv
    Rscript GOTermAnalysis/GoTermAnalysis.R GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}Only/geneList-${comparison2}Only_${f}FC.csv
    Rscript GOTermAnalysis/GoTermAnalysis.R GOTermAnalysis/${comparison1}_${comparison2}/${comparison2}unique_vs_${comparison1}/geneList-${comparison2}unique_vs_${comparison1}_${f}FC.csv
    Rscript GOTermAnalysis/GoTermAnalysis.R GOTermAnalysis/${comparison1}_${comparison2}/${comparison1}unique_vs_${comparison2}/geneList-${comparison1}unique_vs_${comparison2}_${f}FC.csv
done 


#done

fc="0.58" 

Rscript geneExpGraphsFromGOTerms.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv ${comparison1}_${comparison2} $HC1 $HC2 $HC3 $HC4 ${comparison2}Only ${fc}
Rscript geneExpGraphsFromGOTerms.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv ${comparison1}_${comparison2} $HC1 $HC2 $HC3 $HC4 ${comparison1}Only ${fc}
Rscript geneExpGraphsFromGOTerms.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv ${comparison1}_${comparison2} $HC1 $HC2 $HC3 $HC4 ${comparison1}unique_vs_${comparison2} ${fc}
Rscript geneExpGraphsFromGOTerms.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv ${comparison1}_${comparison2} $HC1 $HC2 $HC3 $HC4 ${comparison2}unique_vs_${comparison1} ${fc}
Rscript geneExpGraphsFromGOTerms.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv ${comparison1}_${comparison2} $HC1 $HC2 $HC3 $HC4 ${comparison2}intersect${comparison1} ${fc}





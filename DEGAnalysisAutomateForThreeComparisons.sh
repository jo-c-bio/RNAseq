#!/bin/bash

#This script takes two of the control vs treatment gives you lists of significant DEGs at your specified fold changes 
#that unique or shared between them and builds a venn digram and a heatmap. You also get an exported list of DE Gene names.

#Let's start by nominating our comparisons we'd like to run, we'll use these throughout the rest of the script
#To make directories, and communicate to R which files to use. 

#if you want to specify these here, use:
## comparison1="cWT_vs_vWT"
## comparison2="cPRT6_vs_vPRT6"
## comparison3="cDM_vs_vDM"
#else, use these for the wrapper:
comparison1=${1}
comparison2=${2}
comparison3=${3}

#Then later we need some column names of our dataframes, these need to be the same as in your original dataframes. CH = comparison header 
#if you want to specify these here, use:
## CH1="WT.C"
## CH2="WT.V"
## CH3="prt6.C"
## CH4="prt6.V"
## CH5="DM.C"
## CH6="DM.V"
#else, use these for the wrapper:
CH1=${1}
CH2=${2}
CH3=${3}
CH4=${4}
CH5=${5}
CH6=${6}


#Then let's use these to make some directories for storing all the sections of our comaprion venn diagrams. 
#Notation, 
    #Intersect is all those genes in the middle of the venn diagram
    #Only are all those genes that are within a full circle, so all the genes associated with a comparison whether they overlap with the other comparison not.
    #Unique are the genes either side of the venn diagram, those that are for one comparison but are not shared at all with the other. 


mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}

mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/ALLintersect${comparison1}_${comparison2}_${comparison3}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison1}_${comparison2}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison2}_${comparison3}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison3}_${comparison1}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison1}Only
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison2}Only
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison3}Only
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison1}unique_vs_${comparison2}and${comparison3}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison2}unique_vs_${comparison1}and${comparison3}
mkdir GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison3}unique_vs_${comparison1}and${comparison2}


#Here we loop over 'f', our chosen fold change values we want to analyse at. I chose 1, 1.5 and 2 fold change here.
for f in 1 1.5 2

do
    #variables are: fold-change comparison1 comparison2 file-comparison1 file-comparison2
    #Rscript jo_heatmap_v2.r $f cPRT6_vs_vPRT6 cDM_vs_vDM VRN2_RNA_seq_data_AllTogether/cPRT6vsvPRT6_0FC.csv VRN2_RNA_seq_data_AllTogether/CDMvsVDM_0FC.csv
    Rscript jo_heatmap_for3waycomparison.R $f $comparison1 $comparison2 $comparison3 VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv
done 

for f in 1 1.5 2

do
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/ALLintersect${comparison1}_${comparison2}_${comparison3}/geneList-ALLintersect${comparison1}_${comparison2}_${comparison3}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison1}_${comparison2}/geneList-intersectonly${comparison1}_${comparison2}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison2}_${comparison3}/geneList-intersectonly${comparison2}_${comparison3}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/intersectonly${comparison3}_${comparison1}/geneList-intersectonly${comparison3}_${comparison1}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison1}Only/geneList-${comparison1}Only_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison2}Only/geneList-${comparison2}Only_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison3}Only/geneList-${comparison3}Only_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison1}unique_vs_${comparison2}and${comparison3}/geneList-${comparison1}unique_vs_${comparison2}and${comparison3}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison2}unique_vs_${comparison1}and${comparison3}/geneList-${comparison2}unique_vs_${comparison1}and${comparison3}_${f}FC.csv
         Rscript GOTermAnalysis/GoTermAnalysis.R  GOTermAnalysis/${comparison1}_${comparison2}_${comparison3}/${comparison3}unique_vs_${comparison1}and${comparison2}/geneList-${comparison3}unique_vs_${comparison1}and${comparison2}_${f}FC.csv
done 

Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ALLintersect${comparison1}_${comparison2}_${comparison3}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison1}_${comparison2}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison2}_${comparison3}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison3}_${comparison1}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison1}Only
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison2}Only
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison3}Only
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison1}unique_vs_${comparison2}and${comparison3}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison2}unique_vs_${comparison1}and${comparison3}
Rscript geneExpGraphsFromGoTerms3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison3}unique_vs_${comparison1}and${comparison2}

Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ALLintersect${comparison1}_${comparison2}_${comparison3}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison1}_${comparison2}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison2}_${comparison3}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 intersectonly${comparison3}_${comparison1}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison1}Only
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison2}Only
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison3}Only
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison1}unique_vs_${comparison2}and${comparison3}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison2}unique_vs_${comparison1}and${comparison3}
Rscript geneExpGraphsFromCustomList3way.r VRN2_RNA_seq_data_AllTogether/${comparison1}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison2}_0FC.csv VRN2_RNA_seq_data_AllTogether/${comparison3}_0FC.csv ${comparison1}_${comparison2}_${comparison3} $CH1 $CH2 $CH3 $CH4 $CH5 $CH6 ${comparison3}unique_vs_${comparison1}and${comparison2}

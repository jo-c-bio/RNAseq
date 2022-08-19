# Broad scoping analysis of RNA seq data in WT, vrn2, prt6 and vrn2prt6 mutants of Arabidopsis thaliana, and vernalised (cold-treated) samples of these genotypes.

## Download

Retrieve a local copy with

```sh
git clone https://github.com/jo-c-bio/RNAseq.git
```

## Requirements

Before running, ensure you have access to:
- R (Version 4.0.1 or later)
  - Required packages:
    - ggplot2
    - DESeq2
    - dplyr
    - ggpubr
    - ggvenn
    - pheatmap
    - reshape2
    - stringr

## Directory set-up

All main scripts live within the top directory, except GoTermAnalysis.R, which is in the GoTermAnalysis file, where most of the output will go. 

## General analysis

The main script is `DEGAnalysisAutomate.sh` which takes as input two comparisons between genotypes. These are also the names of two DEseq2 files that are the raw data for this analysis. For example, 
`comparison1="cWT_vs_cVRN2" comparison2="cDM_vs_cPRT6"`
You can modify these variables in the script.

Also, we'll change the column heading variable for these comparisons, for example,
`HC1="WT.C"`
`HC2="vrn2.C"`
`HC3="DM.C"`
`HC4="prt6.C"`
mutants are lowercase, DM and WT upper. .C or .V is control or vernalised.

For the upcoming scripts, I've implemented a naming system for the gene list outputs. These are:

- `${comparison2}intersect${comparison1}` The genes differentially expressed in both comparisons, either up or down.
- `${comparison1}Only` All genes differentially expressed in comparison 1, whether they are also up or down in the other comparison or not.
- `${comparison2}Only` All genes differentially expressed in comparison 2, whether they are also up or down in the other comparison or not.
- `${comparison2}unique_vs_${comparison1}` The differentially expressed genes that are unique to comparison 2, not shared or present in/with comparison 1
- `${comparison1}unique_vs_${comparison2}` The differentially expressed genes that are unique to comparison 1, not shared or present in/with comparison 2
(note: all genes should be significant, and pass either above or below the log2 fold change threshold (the gene lists are not split into up or down, but the `jo_heatmap_v2.r` code is editable to do so))

The script runs three .R scripts, the first being `jo_heatmap_v2.r`. It takes the raw comparison files, comparison names and a fold change value. The script filters through only genes from the data files that pass a threshold of padj < 0.05, and a Log2 fold change + or - the fold change (that has been input as a variable). It then creates gene lists of those that are either overlapping or specific to the comparisons. These gene lists can be altered, but the output is used in the two subsequent scripts. It then makes heatmaps by plotting the fold change values of the intersecting gene list of the two comparions, or those gene lists that are unique to either comparison. 

The next script is `GoTermAnalysis/GoTermAnalysis.R` which takes any gene list, and runs GO Term analysis on it. The `DEGAnalysisAutomate.sh` script tells it to run the analysis on all 5 types of gene list we have generates. The script outputs three file types. These are GO terms for each catergory of Biological Process, Cellular Compartment and Molecular Function GO Term types. For each type, we have a file with the GO Term and any genes associated with it from our input file (`result_GO-Individual_Genes`), and the other is the statistics and number of significant genes for each GO Term (`result_GO-geneList`). 

The last script is `geneExpGraphsFromGOTerms.r`. This takes the `result_GO-Individual_Genes` lists for any terms that have over a specified number of significant genes (this is currently set to 10, but editable in the script). It is also possible to specify which FC filetred gene list you use here, change the fc variable in `DEGAnalysisAutomate.sh`. It then outputs one graph per GO term, with the genes associated with that GO Term on the x-axis, and the mean +-SD (across biological replicates) of the normalised read count data from the raw comparison files on the y-axis. 




To run the wrapper script, use:

Rscript `DEGAnalysisAutomate.sh` [comparison 1 name] [comparison 2 name] [column head comp1part1] [column head comp1part2] [column head comp2part1] [column head comp1part2]

For example 

```sh
./DEGAnalysisAutomate.sh cWT_vs_cVRN2 cDM_vs_cPRT6 WT.C vrn2.C DM.C prt6.C
```

would run the three scripts to show DEGs shared and unique to vrn2 compared to col-0, and the prt6 mutant compared to the double mutant

## Other Scripts

There is a series of other useful scripts available.

### `DEGAnalysisAutomateForThreeComparisons.sh`
If you'd like to make a 3-way comparison instead, you can use the scripts:
    - `DEGAnalysisAutomateForThreeComparisons.sh`
    which runs:
        -`jo_heatmap_for3waycomparison.R`
        -the original `GoTermAnalysis.R`
        -`geneExpGraphsFromGoTerms3way.r`

it runs in a very similar way to `DEGAnalysisAutomate.R`, for example:
```sh
    ./DEGAnalysisAutomateForThreeComparisons.sh cWT_vs_vWT cPRT6_vs_vPRT6 cDM_vs_vDM WT.C WT.V prt6.C prt6.V DM.C DM.V
```

### `PCA_intitial.R`

Generates PCA and scree plots for looking at all biological reps and all treatments to judge the similarity in variance across them. Run manually. 

### `HypoxiaDataAddition.R`
We incorporated DEGs from a microarray experiment [1] into our comparison, and the manual code is available in `HypoxiaDataAddition.R`. We didn;t take it further than gene lists as Mircoarray fold changes are not directly comparable to RNAseq fold changes because of the methodology. 

### `DE_genes_automate.R`

A simpler script that will give you gene lists of UP and DOWN DEGs in terms of each simple comparison at a time. 
Note: it filters by log2 fold change but the FILE NAME reports the actual fold change. Also, it filters by 0.05 padj values and takes out and padj values with #N/A.
run using:

e.g log2 fold change filter on genes up in vrn2 mutant  -> "cWT_vs_cVRN2_UPby4.csv"

```R
Rscript DE_genes_automate.R cWT_vs_vWT VRN2_RNA_seq_data_AllTogether/cWT_vs_vWT_0FC.csv 
``` 
## References 

[1] https://www.nature.com/articles/nature10534#Sec14


CWTvsVWT
Cvrn2vsVvrn2
Cprt6vsVprt6
CDMvsVDM
VWTvsVDM
VWTvsVprt6
VWTvsVvrn2
Rscript DE_genes_automate.R cWT_vs_vWT VRN2_RNA_seq_data_AllTogether/cWT_vs_vWT_0FC.csv 
Rscript DE_genes_automate.R cVRN2_vs_vVRN2 VRN2_RNA_seq_data_AllTogether/cVRN2_vs_vVRN2_0FC.csv 
Rscript DE_genes_automate.R cPRT6_vs_vPRT6 VRN2_RNA_seq_data_AllTogether/cPRT6_vs_vPRT6_0FC.csv 
Rscript DE_genes_automate.R cDM_vs_vDM VRN2_RNA_seq_data_AllTogether/cDM_vs_vDM_0FC.csv
Rscript DE_genes_automate.R vWT_vs_vDM VRN2_RNA_seq_data_AllTogether/vWT_vs_vDM_0FC.csv 
Rscript DE_genes_automate.R vWT_vs_vPRT6 VRN2_RNA_seq_data_AllTogether/vWT_vs_vPRT6_0FC.csv 
Rscript DE_genes_automate.R vWT_vs_vVRN2 VRN2_RNA_seq_data_AllTogether/vWT_vs_vVRN2_0FC.csv 


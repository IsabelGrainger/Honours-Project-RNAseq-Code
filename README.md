# Honours-Project-RNAseq-Code
This code was used to perform differential gene expression analysis of three publicly available datasets comparing _Arabidopsis_ gene expression between lines active and defective in a component of the m6A writer complex.

The analysis of each dataset was performed separately, and the resulting logFC values and volcano plots were interpreted in my project "Why is the _Arabidopsis thaliana FLC_ locus so intensely studied?"

The datasets were converted into salmon .quant.sf files by Dr Matthew Parker, and then I imported them into RStudio for analysis using tximport.

For this code to work, I first set up the correct R environment using the following commands:

setwd("C:/Users/Isabel/OneDrive - University of Dundee/Documents/Honours Project/Hon Analysis/Bioinformatics/data")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

require(edgeR)

require(EDASeq)

require(ggplot2)

Then the rest of the code in the text file "m6A_DGE_script" can be run sequentially.

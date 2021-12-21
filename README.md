# Honours-Project-RNAseq-Code
Code used to perform differential gene expression analysis of three publicly available datasets comparing Arabidopsis gene expression between lines active and defective in a component of the m6A writer complex.

The analysis of each dataset was performed separately, and the resulting logFC values and volcano plots were interpreted in my project "Why is the _Arabidopsis thaliana FLC_ locus so intensely studied?"

### Setting up the environment  
setwd("C:/Users/Isabel/OneDrive - University of Dundee/Documents/Honours Project/Hon Analysis/Bioinformatics/data")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

require(edgeR)
require(EDASeq)
require(ggplot2)

### import mta
library(tximport)
dir <- "C:/Users/Isabel/OneDrive - University of Dundee/Documents/Honours Project/Hon Analysis/Bioinformatics/data"
list.files(dir)

files <- list.files(file.path(dir, "mta_col_quants"))
mtafiles <- file.path(dir, "mta_col_quants", files)

### define gene / transcript ID column
library("tidyverse")

col0_1_quant <- read_delim("mta_col_quants/col0_1.quant.sf", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

col0_1_quant$geneid <- substring(col0_1_quant$Name, 1,9)

### create mta data files
TLmtadata <- tximport(mtafiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")],txOut = TRUE)
GLmtadata <- tximport(mtafiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")])

### TL mta matrix
TLmtaMatrix<- as.matrix(TLmtadata$counts)
dimnames(TLmtaMatrix)[[2]]<-gsub(".quant.sf","",files)
dimnames(TLmtaMatrix)
TLmtaMatrix[rownames(TLmtaMatrix)=="AT5G10140.1",]

### mta FLC transcripts
TLmtaMatrix[substring(rownames(TLmtaMatrix),1,9)=="AT5G10140",]

#### mta COOLAIR transcripts
TLmtaMatrix[substring(rownames(TLmtaMatrix),1,9)=="AT5G01675",]

### GL mta matrix
GLmtaMatrix<- as.matrix(GLmtadata$counts)
dimnames(GLmtaMatrix)[[2]]<-gsub(".quant.sf","",files)
dimnames(GLmtaMatrix)
GLmtaMatrix[rownames(GLmtaMatrix)=="AT5G10140",]

### sanity check
plotPCA(TLmtaMatrix)
plotPCA(GLmtaMatrix)

### GL mta dge###
cond <- as.factor(substr(colnames(GLmtaMatrix),1,1))
dge <- DGEList(GLmtaMatrix,group= cond)
dge <- calcNormFactors(dge, method="TMM")
design <- model.matrix( ~dge$sample$group )
dge <- estimateDisp(dge, design)

fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=2)

### extracting logFC values
res_LRT_tb <- lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_LRT_tb
res_LRT_tb[res_LRT_tb$gene == "AT5G10140","logFC"]

res_LRT_tb
res_LRT_tb[res_LRT_tb$gene == "AT5G01675","logFC"]

### volcano plot
tt <- topTags(lrt,n=32310)$table
names <- rownames(tt)
rownames(tt) <- NULL
tt2<- cbind(names,tt)

FLCtt2 <- tt2 %>% filter(names == "AT5G10140")
COOLAIRtt2 <- tt2 %>% filter(names == "AT5G01675") 
ggplot(tt2, aes(x=logFC,y=-log(FDR),color=logCPM)) + 
  scale_colour_gradientn(colours=c("#0033CC" ,"#66CCFF", "#FFFFFF" )) +
  geom_point() + 
  geom_point(data = FLCtt2, col = "red", size = 3) +
  geom_point(data = COOLAIRtt2, col = "orange", size = 3)


### importing vir
vfiles <- list.files(file.path(dir, "vir_VIRc_quants"))
virfiles <- file.path(dir, "vir_VIRc_quants", vfiles)
GLvirdata <- tximport(virfiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")])
TLvirdata <- tximport(virfiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")],txOut = TRUE)

### TL vir matrix
TLvirMatrix<- as.matrix(TLvirdata$counts)
dimnames(TLvirMatrix)[[2]]<-gsub(".quant.sf","",vfiles)
dimnames(TLvirMatrix)
TLvirMatrix[rownames(TLvirMatrix)=="AT5G10140.1",]

### vir FLC transcripts
TLvirMatrix[substring(rownames(TLvirMatrix),1,9)=="AT5G10140",]

### vir COOLAIR transcripts
TLvirMatrix[substring(rownames(TLvirMatrix),1,9)=="AT5G01675",]

### GL vir matrix
GLvirMatrix<- as.matrix(GLvirdata$counts)
dimnames(GLvirMatrix)[[2]]<-gsub(".quant.sf","",vfiles)
dimnames(GLvirMatrix)
GLvirMatrix[rownames(GLvirMatrix)=="AT5G10140",]

###sanity check
plotPCA(TLvirMatrix)
plotPCA(GLvirMatrix)

### GL vir dge
cond <- as.factor(substr(colnames(GLvirMatrix),1,1))
dge <- DGEList(GLvirMatrix,group= cond)
dge <- calcNormFactors(dge, method="TMM")
design <- model.matrix( ~dge$sample$group )
dge <- estimateDisp(dge, design)

fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=2)

### extracting logFC values
res_LRT_tb <- lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_LRT_tb
res_LRT_tb[res_LRT_tb$gene == "AT5G10140","logFC"]

res_LRT_tb[res_LRT_tb$gene == "AT5G01675","logFC"]

### volcano plot
topTags(lrt)
tt <- topTags(lrt,n=32310)$table
names <- rownames(tt)
rownames(tt) <- NULL
tt2<- cbind(names,tt)

FLCtt2 <- tt2 %>% filter(names == "AT5G10140")
COOLAIRtt2 <- tt2 %>% filter(names == "AT5G01675")
ggplot(tt2, aes(x=logFC,y=-log(FDR),color=logCPM)) + 
  scale_colour_gradientn(colours=c("#0033CC" ,"#66CCFF", "#FFFFFF" )) +
  geom_point() + 
  geom_point(data = FLCtt2, col = "red", size = 3) +
  geom_point(data = COOLAIRtt2, col = "orange", size = 3)


### importing fip37
ffiles <- list.files(file.path(dir, "fip37_col_quants"))
fipfiles <- file.path(dir, "fip37_col_quants", ffiles)
TLfipdata <- tximport(fipfiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")],txOut = TRUE)
GLfipdata <- tximport(fipfiles, type="salmon", countsCol = "NumReads", tx2gene = col0_1_quant[,c("Name","geneid")])

### TL fip37 matrix
TLfipMatrix<- as.matrix(TLfipdata$counts)
dimnames(TLfipMatrix)[[2]]<-gsub(".quant.sf","",ffiles)
dimnames(TLfipMatrix)
TLfipMatrix[rownames(TLfipMatrix)=="AT5G10140",]

### FIP37 / fip37 FLC transcripts
TLfipMatrix[substring(rownames(TLfipMatrix),1,9)=="AT5G10140",]

#### FIP37 / fip37 COOLAIR transcripts
TLfipMatrix[substring(rownames(TLfipMatrix),1,9)=="AT5G01675",]

### GL fip37 matrix
GLfipMatrix<- as.matrix(GLfipdata$counts)
dimnames(GLfipMatrix)[[2]]<-gsub(".quant.sf","",ffiles)
dimnames(GLfipMatrix)
GLfipMatrix[rownames(GLfipMatrix)=="AT5G10140",]

### sanity check
plotPCA(TLfipMatrix)
plotPCA(GLfipMatrix)

### fip37 dge###
cond <- as.factor(substr(colnames(GLfipMatrix),1,1))
dge <- DGEList(GLfipMatrix,group= cond)
dge <- calcNormFactors(dge, method="TMM")
design <- model.matrix( ~dge$sample$group )
dge <- estimateDisp(dge, design)

fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=2)

### extracting logFC values
res_LRT_tb <- lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_LRT_tb
res_LRT_tb[res_LRT_tb$gene == "AT5G10140","logFC"] 

res_LRT_tb[res_LRT_tb$gene == "AT5G01675","logFC"]

### volcano plot
tt <- topTags(lrt,n=10000)$table
names <- rownames(tt)
rownames(tt) <- NULL
tt2<- cbind(names,tt)

FLCtt2 <- tt2 %>% filter(names == "AT5G10140")
COOLAIRtt2 <- tt2 %>% filter(names == "AT5G01675")
ggplot(tt2, aes(x=logFC,y=-log(FDR),color=logCPM)) + 
  scale_colour_gradientn(colours=c("#0033CC" ,"#66CCFF", "#FFFFFF" )) +
  geom_point() + 
  geom_point(data = FLCtt2, col = "red", size = 3) +
  geom_point(data = COOLAIRtt2, col = "orange", size = 3)

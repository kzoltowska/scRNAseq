#Loading packages needed for the analysis ----
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)

#Loading data & building Seurat objects ----

mutant<-Read10X("~/Single_cell_RNAseq_scripts/data/MUT/")
wt<-Read10X("~/Single_cell_RNAseq_scripts/data/WT/")
se_mutant<-CreateSeuratObject(counts=mutant,min.cells = 3, min.features = 200)
se_wt<-CreateSeuratObject(counts=wt, min.cells = 3, min.features = 200)
se_mutant$Condition<-"mutant"
se_wt$Condition<-"wt"

#Assembling Seurat object list for easy processing ----
se_list<-list("WT"=se_wt, "Mutant"=se_mutant)

#Looking at mitochondria transcripts that could indicate cell death ----
se_list<-lapply(se_list, function(x){
  x$mitochondrial<-PercentageFeatureSet(x, pattern="^mt")
  return(x)
})

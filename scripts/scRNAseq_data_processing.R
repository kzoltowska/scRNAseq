#Loading packages needed for the analysis ----
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)

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

df_mitochondrial<-rbind(se_list$WT@meta.data %>% dplyr::select(Condition, mitochondrial),
                        se_list$Mutant@meta.data %>% dplyr::select(Condition, mitochondrial))
df_mitochondrial$Condition<-as.factor(df_mitochondrial$Condition)

ggplot(df_mitochondrial, aes(x=Condition, y=mitochondrial)) +geom_violin(fill="cornflowerblue") + theme_pubr()+
  ylab("% mitochondrial transcript")

#Filtering cells with high percentage (>25%) of mitochondrial transcripts----

se_list_filtered<-lapply(se_list, function(x){
  x<-subset(x,mitochondrial < 25)
  return(x)
})

#Cells with low feature count can cluster differently thus it is worth assigning them to a seperate category for visualisation ----
#Here low feature is considered as 300 - note that while making Seurat object 200 was considered as min.features
se_list_filtered<-lapply(se_list_filtered, function(x){
  x$low_feature<-x$nFeature_RNA<300
  return(x)
})

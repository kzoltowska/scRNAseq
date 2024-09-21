#Loading packages needed for the analysis ----
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(clustree)

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

se_list<-lapply(se_list, function(x){
  x<-subset(x,mitochondrial < 25)
  return(x)
})

#Cells with low feature count can cluster differently thus it is worth assigning them to a seperate category for visualisation ----
#Here low feature is considered as 300 - note that while making Seurat object 200 was considered as min.features
se_list<-lapply(se_list, function(x){
  x$low_feature<-x$nFeature_RNA<300
  return(x)
})

#Data preprocessing - normalisation, finding variable features and scaling----
#Normalisation:“LogNormalize”: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
#This is then natural-log transformed using log1p
se_list_filtered<-se_list
se_list_filtered<-lapply(se_list_filtered, NormalizeData)

#Find variable features----
#by default 2000 features are found
se_list_filtered<-lapply(se_list_filtered, FindVariableFeatures)

#Scaling the data ----
se_list_filtered<-lapply(se_list_filtered, ScaleData)

#Run PCA----
#By default 50 PCs are calculated
#The data will be stored in @reductions under a reduction name pca
se_list_filtered <- lapply(se_list_filtered, RunPCA)

#There are 2 strategies to know how many PCs to include in the further analysis:
#heatmaps - seem to be more intuitive
#elbow plots
#DimHeatmap function is used to draw the heatmaps
#Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. 
#Allows for nice visualization of sources of heterogeneity in the dataset.

lapply(se_list_filtered, function(x){
  DimHeatmap(x, dims=1:15, cells=500)
  })

lapply(se_list_filtered, ElbowPlot)

#Considering the output of the heatmap and elbow plot 14 PCs will be used further
#Next step aims at finding nearest neighbours

se_list_filtered<-lapply(se_list_filtered, function(x){
  x<-FindNeighbors(x, dims=1:14, graph.name=c("RNA_nn", "RNA_snn"))
})

#Finding clusters
#Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
se_list_filtered<-lapply(se_list_filtered, function(x){
  x<-FindClusters(x, resolution=0.8, graph.name="RNA_snn")
  x<-FindClusters(x, resolution=0.4, graph.name="RNA_snn")
  x<-FindClusters(x, resolution=0.2, graph.name="RNA_snn")
  x<-FindClusters(x, resolution=0.1, graph.name="RNA_snn")
})

#looking at clustering based on resolution
#clustree gives nice visual insight how the assignment to clusters change based on resolution
lapply(se_list_filtered, function(x){
  clustree(x@meta.data, prefix="RNA_snn_res.")
}
)

#calculate UMAP and plot clusters ----
#it is worth taking a look at different resolutions
se_list_filtered<-lapply(se_list_filtered, function(x){
  x<-RunUMAP(x, dims=1:14, n.neighbors=20)
})

lapply(se_list_filtered,function(x){
DimPlot(x, reduction="umap", group.by = "RNA_snn_res.0.4")
}
)

#merge the data into single merged suerat object and re-perform the analysis
#perform all the steps as before but now on merge
SOM<-merge(se_list$WT, se_list$Mutant,
           add.cell.ids=c(1,2))

SOM<-NormalizeData(SOM)

SOM <- FindVariableFeatures(SOM)

SOM <- ScaleData(SOM)

SOM<-RunPCA(SOM)

DimHeatmap(SOM, dims=1:15, cells=500)
 
SOM<-FindNeighbors(SOM, dims=1:14, graph.name = c("nn_som", "snn_som"))

SOM<-FindClusters(SOM, graph.name = "snn_som", resolution = 0.1)
SOM<-FindClusters(SOM, graph.name = "snn_som", resolution = 0.2)
SOM<-FindClusters(SOM, graph.name = "snn_som", resolution = 0.4)
SOM<-FindClusters(SOM, graph.name = "snn_som", resolution = 0.8)

clustree(SOM, prefix="snn_som_res.")

SOM<-RunUMAP(SOM, dims=1:14)

DimPlot(SOM, reduction = "umap", group.by = "Condition")
DimPlot(SOM, reduction = "umap", group.by = "snn_som_res.0.8")
DimPlot(SOM, reduction = "umap", group.by = "snn_som_res.0.4")
DimPlot(SOM, reduction = "umap", group.by = "snn_som_res.0.2")
DimPlot(SOM, reduction = "umap", group.by = "snn_som_res.0.1")

#Clusters need to be assigned before proceeding----
Idents(SOM)<-SOM$snn_som_res.0.1

#find markers ----
SOM_joined<-JoinLayers(SOM)
markers<-FindAllMarkers(SOM_joined, min.diff.pct = 0.5, only.pos = TRUE,
                 max.cells.per.ident = 100)

#Look at the markers using dotplot ----
genes <- slice_head(group_by(markers,cluster),n=3)$gene
DotPlot(SOM_joined,features=genes,group.by="snn_som_res.0.1",
        dot.scale=4,cols="RdYlBu") +
  theme(axis.text.x=element_text(angle=45,hjust=1))



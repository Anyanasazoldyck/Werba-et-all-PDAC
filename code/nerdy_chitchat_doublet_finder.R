
library(Seurat)
#install doublet finder 
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)
library(tidyr)
seu=objectsseurat_list_no_liver$P03
#filter low quality cells 
seu[["percent.mt"]] = PercentageFeatureSet(seu, pattern = "^MT")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu= subset(seu, nFeature_RNA >500 & nCount_RNA >1500 & percent.mt <15)


#normalize 
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)
ElbowPlot(seu)

pcs_use=1:10

# ---- chose PC "Squid Stat" ----
ElbowPlot(seu)

#find seurat clusters  
seu <- FindNeighbors(seu, dims = pcs_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu <- RunUMAP(seu, dims =pcs_use, spread = 8, min.dist = 0.2)
DimPlot(seu, reduction = "umap", pt.size = 0.5)

# PK
sweep.res <- paramSweep(seu, PCs = pcs_use, sct = F) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
best.pk <- find.pK(sweep.stats )
ggplot(best.pk, aes(x=pK, y=BCmetric))+
  geom_point()+
  geom_line()
my_best_pk= 0.15
#estimate the multiplate rate 
length(colnames(seu))
multiplet_rate= 0.08
 

# ---- Homotypic adjustment ----
annotations <- seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(multiplet_rate * nrow(seu@meta.data)) # multiply by number of cells to get the number of expected multiplets 
nExp.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets


# ---- Run DoubletFinder ----
print(">> Running DoubletFinder/mnt/d/Command Line Tutorials/raw_data/data_d4.")
seu <- DoubletFinder::doubletFinder(
  seu = seu,
  PCs = pcs_use, 
  pN = 0.25, 
  pK = my_best_pk,
  nExp = nExp.adj, 
  sct = FALSE
)



meta_unfiltered = seu@meta.data
seu <- subset(seu, cells = colnames(seu)[seu@meta.data[["DF.classifications_0.25_0.15_1145"]] == "Singlet"])



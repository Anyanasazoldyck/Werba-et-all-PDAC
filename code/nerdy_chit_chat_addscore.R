setwd("D:/PDAC_SCRNA")


#libraries##############
library(Seurat)
library(readxl)



###############################################################################################
#----------------moffitt----------------------
###############################################################################################
#subset the epithelial cells-----

sc.subset <- subset(all_samples, subset =cell.types %in% c("Epithelial", "Prolif_Epithelial"))
#recluster


sc.subset <- NormalizeData(sc.subset, normalization.method = "LogNormalize")
sc.subset <- FindVariableFeatures(sc.subset, selection.method = "vst", nfeatures = 2000)
sc.subset <- ScaleData(sc.subset)
sc.subset <- RunPCA(sc.subset, features = VariableFeatures(object = sc.subset))
ElbowPlot(sc.subset)
dims_to_use <- 1:20
sc.subset <- FindNeighbors(sc.subset, dims = dims_to_use)
sc.subset <- FindClusters(sc.subset, resolution = 0.5)
sc.subset <- RunUMAP(sc.subset, dims = dims_to_use, spread = 3, min.dist = 0.05)

DimPlot(sc.subset, group.by = "orig.ident")



#moffitt signiture--------
moffitt <- as.data.frame(read_excel("D:/PDAC_SCRNA/csv/moffitt.xlsx"))
colnames(moffitt)<- c("Classical","Basal")
moffitt <- moffitt[-c(1),]


gs <- list(Classical = moffitt$Classical,
           Basal     = moffitt$Basal)

sc.subset <- AddModuleScore(sc.subset, features = gs, ctrl = 100,
                            name = c("Classical","Basal"))

sc.subset$Classical_minus_Basal <- sc.subset$Classical1 - sc.subset$Basal2
sc.subset$Subtype <- ifelse(sc.subset$Classical_minus_Basal > 0, "Classical","Basal")


DimPlot(sc.subset, group.by = "Subtype")


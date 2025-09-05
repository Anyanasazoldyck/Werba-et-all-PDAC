setwd("D:\\PDAC_SCRNA")



library(Seurat)
library(readxl)
library(stringr)
library(ggplot2)
library(Azimuth)
#set files ####
plots="plots"
dir.create(plots, recursive = T)
#-------------------------------
#-------Load data --------

#load the data ####



ss <- read_xlsx("ss.xlsx")
ss <- ss[!ss$Procedure == "liver", ]  

# extract the file names
all_files <- list.files("GSE205013_RAW", full.names = TRUE)
sample_ids <- str_extract(all_files, "P[0-9]+")
unique_samples <- unique(sample_ids)


seurat_list <- list()

for (p in unique_samples) {
  
  feat_file = all_files[grepl(paste0(p,".*features.tsv.gz"),all_files)]
  bc_file = all_files[grepl(paste0(p,".*barcodes.tsv.gz"),all_files)]
  mtx_file = all_files[grepl(paste0(p,".*matrix.mtx.gz"),all_files)]

  data <- ReadMtx  (
    mtx  = mtx_file,
    features = feat_file,
    cells = bc_file
  )
  sc_obj= CreateSeuratObject(counts = data,project = p,assay = "RNA")
  
  #add metadata 
  meta_row = ss[ss$Sample_ID==p,] # the reason is to match meta to p
  
  sc_obj$Treatment = meta_row$Treatment
  sc_obj$Moffitt = meta_row$Moffitt
  sc_obj$Stage = meta_row$Stage
  
  
  seurat_list [[p]] <- sc_obj

}


seurat_list_no_liver = seurat_list[ss$Sample_ID] #I wanted to filter out liver 


#merge 
all_samples = merge(x=seurat_list_no_liver[[1]],y=seurat_list_no_liver[-1])


#---------------quality control ---------------------------

all_samples[["percent.mt"]] <- PercentageFeatureSet(all_samples, pattern = "^MT-")
p= VlnPlot(all_samples, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3, pt.size = 0)
png(paste0(plots,"/PreQC.png"), res=300, height = 6*300, width = 15*300)
p
dev.off()

plot1 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
png(paste0(plots,"/PreQC_scatter.png"), res=300, height = 6*300, width = 15*300)
plot1 + plot2
dev.off()



all_samples =subset(all_samples, subset = nFeature_RNA > 500  & nCount_RNA >1500 & percent.mt <15)
p= VlnPlot(all_samples, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3, pt.size = 0)
png(paste0(plots,"/PostQC.png"), res=300, height = 6*300, width = 15*300)
p
dev.off()

plot1 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
png(paste0(plots,"/PostQC_scatter.png"), res=300, height = 6*300, width = 15*300)
plot1 + plot2
dev.off()



#----------doublet_finder---------------------------------
res_list <- lapply(names(samp_split), function(nm) {
  out <- run_doubletfinder_lognorm(samp_split[[nm]])
  out
})
names(res_list) <- names(samp_split)
saveRDS(res_list,"res_doublet_filtered.rds")

library(Seurat)



all_samples = merge(x=seu_list[[1]], y=seu_list[-1])

meta=all_samples@meta.data 
write.csv(meta,"metadata_singletpredics.csv")
col=colnames(meta)
unwantedcols1=grepl("DF.classifications|pANN", col)
meta <- meta[, !unwantedcols1]
all_samples@meta.data =meta

saveRDS(all_samples,"all_samples_no_doublets.rds")
#------------normalize data ------------------------------

all_samples <- NormalizeData(all_samples, normalization.method = "LogNormalize", scale.factor = 10000)
all_samples <- FindVariableFeatures(all_samples, selection.method = "vst", nfeatures = 2000)
all_samples <- ScaleData(all_samples)



#-----------------PCA---------------------------------------
all_samples <- RunPCA(all_samples, reduction.name = "unintegrated_PCA")


png("plots/PCA_elbow.png")
ElbowPlot(all_samples, reduction = "unintegrated_PCA")
dev.off()

stdev= all_samples@reductions$unintegrated_PCA@stdev
percent_var= stdev^2 / sum(stdev^2)*100

png("plots/PCA_Heatmaps.png", res = 300, width = 300*8, height = 300*8)
DimHeatmap(all_samples, reduction = "unintegrated_PCA", dims = 1:10)
dev.off()


dims_to_use= 1:5

#pre integration umap -------------------------------------
all_samples= FindNeighbors(all_samples, dims_to_use, reduction = "unintegrated_PCA")
all_samples <- FindClusters(all_samples, resolution = 0.2)
all_samples <- RunUMAP(
  all_samples,
  reduction = "unintegrated_PCA",
  dims = dims_to_use,                
  reduction.name = "unintegrated_umap",
  min.dist = 0.2,
  spread = 8,
  verbose = TRUE)
  
umap_theme = theme(plot.title = element_text(hjust = 0.5), 
                   legend.position = "none",  
                   axis.text = element_blank(), 
                   axis.title = element_text(size = 12), 
                   text = element_text(size = 12, family = "ArialMT"),
                   axis.ticks = element_blank())
png("plots/unintegrated_umap.png", res = 300, width = 300*6, height = 300*6)
DimPlot(all_samples, reduction = "unintegrated_umap")
dev.off()
  
  
#--------------------integration------------------------------


all_samples_integrated <- IntegrateLayers(
  object = all_samples, method = CCAIntegration,
  orig.reduction = "unintegrated_PCA", new.reduction = "integrated_PCA",
  verbose = T
)

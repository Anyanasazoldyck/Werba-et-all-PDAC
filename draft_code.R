setwd("D:\\PDAC_SCRNA")



library(Seurat)
library(readxl)
library(stringr)
#set files ####
plots="plots"
dir.create(plots, recursive = T)

objects_file= "objects"
dir.create(objects_file, recursive = T)

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


#---------------merge----------------------------------- 
all_samples = merge(x=seurat_list_no_liver[[1]],y=seurat_list_no_liver[-1])
saveRDS(all_samples, paste0(objects_file,"/all_samples.rds"))

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

#------------filtering blood cells ---------------------
all_samples[["percent.eryth"]] <- PercentageFeatureSet(all_samples, features  = c("HBA1", "HBA2", "HBB", "HBM",  "ALAS2"))
p= VlnPlot(all_samples, features = "percent.eryth", pt.size = 1)
png(paste0(plots,"/PretQC_bloodcellexp.png"), res=300, height = 6*300, width = 15*300)
p
dev.off()
all_samples=subset(all_samples,subset = percent.eryth <1)

saveRDS(all_samples, paste0(objects_file,"/all_samples_filtered.rds"))


#------------normalize data ------------------------------

all_samples <- NormalizeData(all_samples, normalization.method = "LogNormalize", scale.factor = 10000)
all_samples <- FindVariableFeatures(all_samples, selection.method = "vst", nfeatures = 2000)
all_samples <- ScaleData(all_samples, vars.to.regress = c("percent.mt"), features = VariableFeatures(all_samples) )
saveRDS(all_samples, paste0(objects_file,"/all_samples_scaled.rds"))





setwd("D:\\PDAC_SCRNA")


library(ggplot2)
library(Seurat)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library (ggplot2)
library(patchwork)
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

saveRDS(seurat_list_no_liver,paste0(objects_file,"seurat_list_no_liver.rds"))
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

#-----------------doublet finder for each sample individually-----------------

samp_split <- SplitObject(all_samples)

res_list <- lapply(names(samp_split), function(nm) {
  out <- run_doubletfinder_lognorm(samp_split[[nm]])
  out
})
names(res_list) <- names(samp_split)
saveRDS(res_list,"res_doublet_filtered.rds")




#------------normalize data ------------------------------

all_samples <- NormalizeData(all_samples, normalization.method = "LogNormalize", scale.factor = 10000)
all_samples <- FindVariableFeatures(all_samples, selection.method = "vst", nfeatures = 2000)
all_samples <- ScaleData(all_samples, vars.to.regress = c("percent.mt"), features = VariableFeatures(all_samples) )
saveRDS(all_samples, paste0(objects_file,"/all_samples_scaled.rds"))



#pre integration vizulization (This part was not used in the final run )-------

#-------------pre integration PCA---------------------------------------
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
                   
                   axis.text = element_blank(), 
                   axis.title = element_text(size = 12), 
                   text = element_text(size = 12, family = "ArialMT"),
                   axis.ticks = element_blank())
png("plots/unintegrated_umap.png", res = 300, width = 300*6, height = 300*6)
DimPlot(all_samples, reduction = "unintegrated_umap")
dev.off()





#integration --------------------------------
samples_list <- SplitObject(all_samples, split.by = "orig.ident")
#rename projects 
# Give each Seurat object a unique project name
for (i in seq_along(samples_list)) {
  samples_list[[i]]@project.name <- names(samples_list)[i]
}
#normalize 
# normalize and identify variable features for each dataset independently
samples_list <- lapply(samples_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#features that are most variable across all samples 
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 1000)


#use these features to scale and dimentionally reduce the data 
samples_list <- lapply(samples_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#check if there are doublicated cells name 
# for a list of objects
sapply(samples_list, function(x) anyDuplicated(Cells(x)))
#no doublicated cells name 

#integrate using anchores
anchores <- FindIntegrationAnchors(
  object.list = samples_list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:20
)
all_samples_integrated <- IntegrateData(anchorset = anchores, dims = 1:20)
saveRDS(anchores,"my_anchors.rds")
saveRDS(all_samples_integrated,"all_samples_integrated.rds")




#------------post integration processing -------------
all_samples_integrated <- readRDS("D:/PDAC_SCRNA/objects/all_samples_integrated.rds")


DefaultAssay(all_samples_integrated) <- "integrated"
all_samples_integrated <- ScaleData(all_samples_integrated)
all_samples_integrated <- RunPCA(all_samples_integrated)
all_samples_integrated <- FindNeighbors(all_samples_integrated, dims = 1:5)
all_samples_integrated <- FindClusters(all_samples_integrated, resolution = 0.2)
all_samples_integrated <- RunUMAP(all_samples_integrated, dims = 1:5)

p1= DimPlot(all_samples_integrated, reduction = "umap",group.by = "integrated_snn_res.0.2")
p2= DimPlot(all_samples_integrated, reduction = "umap", group.by = "orig.ident")

p1
png("plots/integrated_umap_res_02.png", res = 300, width = 300*10, height = 300*5)
p1+p2
dev.off()

png("plots/integrated_umap_Treatment.png", res = 300, width = 300*5, height = 300*5)
DimPlot(all_samples_integrated, reduction = "umap",group.by = "Treatment")
dev.off()


#find markers for each cluster###################
DefaultAssay(all_samples_integrated) <- "RNA"
all_samples_integrated <- JoinLayers(all_samples_integrated, assay = "RNA")
?JoinLayers
markers <- FindAllMarkers(all_samples_integrated,only.pos = T)
saveRDS(all_samples_integrated,"all_samples_integrated_dim_reduction.rds")

write.csv(markers,"csv/cluster_markers.csv")


#filter markers ------------------------
markers.sig <- markers[,"p_val_adj"<0.05 & "avg_log2FC">1]
topn=markers %>% dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, .by_group = T) %>% dplyr::slice_head(n=5)
topn$gene
p=DoHeatmap(all_samples_integrated, features = topn$gene)

png("plots/heatmap.png", res = 300, width = 300*8, height = 300*8)
p
dev.off()



#make a complex hm
DefaultAssay(all_samples_integrated)<-"RNA"
genes <- unique(topn$gene)
df <- FetchData(all_samples_integrated, vars = genes)
df$cluster <- Idents(all_samples_integrated)[rownames(df)]
colnames(df)
head(df)

df_longer<- pivot_longer(df,cols = !cluster,names_to = "genes", values_to = "expression")


df_longer<-df_longer %>% group_by(cluster,genes) %>% summarise(mean_expression =mean(expression),.groups = "drop")

hm_mtx <- as.data.frame(pivot_wider(df_longer,names_from = cluster, values_from = mean_expression)) 
rownames(hm_mtx)<-hm_mtx$genes
hm_mtx$genes<-NULL
mat <- t(scale(t(as.matrix(hm_mtx))))




library(ComplexHeatmap)
#define color pal 
cols = colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
?colorRampPalette
p=Heatmap(mat, name = "Z-score",
        cluster_rows = T, cluster_columns = F,
        show_row_names = TRUE, show_column_names = TRUE, color_space = cols)
png("plots/complexheatmap.png", res = 300, width = 300*8, height = 300*10)
p
dev.off()

write.csv(topn,"csv/top_cluster_markers.csv")


#Identifying the prolifrating cells 
ki_67 =FeaturePlot(all_samples_integrated, features = "MKI67", pt.size = 1, order = T)&umap_theme


png("plots/ki-67.png", res = 300, width = 300*5, height = 300*5)
ki_67
dev.off()


#identifying each of the cell type 
#Identifying immune related calls on the umap 
T_cell_markers = c("TRBC2","IL7R","CCL5")

t= FeaturePlot(all_samples_integrated, features = T_cell_markers, pt.size = 1, order = F)&umap_theme

png("plots/T_Cell.png", res = 300, width = 300*5, height = 300*5)
t
dev.off()

NK= c("NKG7", "GZMB")
nk= FeaturePlot(all_samples_integrated, features = NK, pt.size = 1, order = F)&umap_theme
png("plots/NK.png", res = 300, width = 300*5, height = 300*5)
nk
dev.off()

Myeloid = c("CD68","CCL12")
m= FeaturePlot(all_samples_integrated, features = Myeloid, pt.size = 1, order = F)&umap_theme

png("plots/Meyloid.png", res = 300, width = 300*5, height = 300*5)
m
dev.off()

#the rest of the process was done on pangloDP


#set the new ident
ss <- read_xlsx("csv/ss.xlsx",sheet = "Sheet2")



new.cluster.ids = ss$Type

names(new.cluster.ids) = levels(all_samples_integrated)
all_samples_integrated = RenameIdents(all_samples_integrated, new.cluster.ids)
all_samples_integrated = AddMetaData(all_samples_integrated, all_samples_integrated@active.ident, col.name = "cell.types")




p=DimPlot(all_samples_integrated, label = T, label.size = 6, repel = T,pt.size = 1)+umap_theme

png("plots/labeled_umap.png", res = 300, width = 300*10, height = 300*10)
p
dev.off()
#cell proportions-----------------
df <- as.data.frame(table(all_samples_integrated$orig.ident, all_samples_integrated$cell.types))
colnames(df)<-c("Sample","Cell_Type","Frequency")
df=df %>% group_by(Sample,Cell_Type) %>% mutate(total_cell_type_per_sample=sum(Frequency))
df<- df %>% group_by(Sample) %>% mutate(total_cells=sum(total_cell_type_per_sample))
df <- df %>% group_by(Sample,Cell_Type) %>% mutate(prportion=total_cell_type_per_sample/total_cells*100)
write.csv(df,"csv/proporion_of_cells_per_sample.csv")

p=ggplot(df,aes(x=Sample,y=prportion,fill = Cell_Type)) +geom_col() +
  xlab("Sample_ID")+
  ylab("Proportion")+my_theme+theme(axis.text.x = element_text(angle = 90))

png("plots/cell_type_proportion.png", res = 300, width = 300*5, height = 300*5)
p
dev.off()

#redraw the heatmap with labeled clusters

#make a complex hm
genes=top_cluster_markers$gene
#make a complex hm
DefaultAssay(all_samples_integrated)<-"RNA"
df <- FetchData(all_samples_integrated, vars = genes)
df$cluster <- Idents(all_samples_integrated)[rownames(df)]
colnames(df)
head(df)

df_longer<- pivot_longer(df,cols = !cluster,names_to = "genes", values_to = "expression")


df_longer<-df_longer %>% group_by(cluster,genes) %>% summarise(mean_expression =mean(expression),.groups = "drop")

hm_mtx <- as.data.frame(pivot_wider(df_longer,names_from = cluster, values_from = mean_expression)) 
rownames(hm_mtx)<-hm_mtx$genes
hm_mtx$genes<-NULL
mat <- t(scale(t(as.matrix(hm_mtx))))




library(ComplexHeatmap)
#define color pal 
cols = colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
?colorRampPalette
p=Heatmap(mat, name = "Z-score",
          cluster_rows = T, cluster_columns = T,
          show_row_names = TRUE, show_column_names = TRUE, color_space = cols)
png("plots/complexheatmap_annotated.png", res = 300, width = 300*8, height = 300*10)
p
dev.off()


#save the object without the integrated layer to reduce size #####

all_samples_no_int_asssay <- all_samples_integrated_labeled
all_samples_no_int_asssay@assays$integrated <- NULL  
saveRDS(all_samples_no_int_asssay,"objects/all_samples_no_integration_layer.rds")
###############################################################################################
#----------------moffitt----------------------
###############################################################################################
#subset the epithelial cells-----
gc()
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


p=DimPlot(sc.subset, group.by = "Subtype", cols = c("purple","pink")) +umap_theme
p1=DimPlot(sc.subset, group.by = "orig.ident",pt.size = 0.3) +umap_theme

png("plots/moffitt.png",width = 300*10, height = 300*5, res = 300)
p+p1
dev.off()


#---scaterplot-----
df <- sc.subset@meta.data %>% select (c("Classical1","Basal2", "Treatment", "orig.ident"))
df_treated <- df [df$Treatment=="treated",]
p1=ggplot(df_treated, aes(x = Basal2, y = Classical1, fill =Treatment )) +
  geom_point(size =0.8, color = "maroon") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +theme_classic()

p1_sample=ggplot(df_treated, aes(x = Basal2, y = Classical1)) +
  geom_point(size =0.8, aes(colour =orig.ident )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set2", name = "Sample") +
  theme_classic()
png("plots/treated_bysample.png",width = 300*5, height = 300*5, res = 300)
p1_sample
dev.off()

df_untreated <- df [df$Treatment=="naive",]
p2=ggplot(df_untreated, aes(x = Basal2, y = Classical1, fill =Treatment )) +
  geom_point(size =0.8, color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +theme_classic()


p2_sample=ggplot(df_untreated, aes(x = Basal2, y = Classical1 )) +
  geom_point(size =0.8, aes(colour =orig.ident)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Sample") +
  theme_classic()
png("plots/untreated_bysample.png",width = 300*5, height = 300*5, res = 300)
p2_sample
dev.off()

png("plots/classicalbasal.png",width = 300*10, height = 300*5, res = 300)
p2+p1
dev.off()

scale_color_manual(values = c(treated="red"))+
  
  ###############################################################################################
#-----------naive VS treated----------------------
###############################################################################################
library(clusterProfiler)
library(msigdbr)

sc.subset <- SetIdent(sc.subset,value = "Treatment" )
unique(sc.subset@active.ident)
treated_vs_naive <- FindMarkers(sc.subset, ident.1 ="treated", ident.2 = "naive"  )

#reactome <- msigdbr(species = "Homo sapiens", collection  = "C2", subcollection = "CP:REACTOME") |>
 # select(gs_name, gene_symbol) |> distinct()
hallmark <- msigdbr(species = "Homo sapiens", collection  = "H") %>%
   select(gs_name, gene_symbol) |> distinct()
treated_vs_naive$score <- sign(treated_vs_naive$avg_log2FC) * (-log10(treated_vs_naive$p_val_adj + 1e-10))
gene_list <- treated_vs_naive$score

names(gene_list) <- rownames(treated_vs_naive)

gene_list<-na.omit(gene_list)

BiocParallel::register(BiocParallel::SerialParam()) #fix serialize error
options(mc.cores = 1)
set.seed(1)

gene_list = sort(gene_list, decreasing = TRUE)
gsea_results=clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  verbose = FALSE
)
gsea_results <-gsea_results@result
gsea_results<- gsea_results[gsea_results$p.adjust<0.05,]
gsea_results$minus_logFDR= -log10(gsea_results$p.adjust)
p=ggplot(gsea_results, aes(
  x = reorder(ID, NES),
  y = NES,
  fill = NES,
  alpha=minus_logFDR,
  #size = `-logFDR`
)) +
  geom_point(shape = 21, colour = "black",size=10) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "red",
    midpoint = 0,
    name = "Enrichment"
  ) + scale_alpha(range = c(0.3, 1),
                  name = "LogFDR") +
  labs(
    x = "Pathway",
    y = "NES",
    title = paste0("Treated Vs. Naive")
  ) +
  coord_flip() +my_theme
png("plots/halmark.png", width = 400*6, height = 300*6 , res = 300)
write.csv(gsea_results,"csv/hallmark_ranked_by_score.csv")
p
dev.off()



#cell proportions by subtype------------
df <- as.data.frame(table(sc.subset$orig.ident, sc.subset$Subtype))
colnames(df)<-c("Sample","Cell_Type","Frequency")
df=df %>% group_by(Sample,Cell_Type) %>% mutate(total_cell_type_per_sample=sum(Frequency))
df<- df %>% group_by(Sample) %>% mutate(total_cells=sum(total_cell_type_per_sample))
df <- df %>% group_by(Sample,Cell_Type) %>% mutate(prportion=total_cell_type_per_sample/total_cells*100)
write.csv(df,"csv/proporion_of_classicalbasal_per_sample.csv")

p=ggplot(df,aes(x=Sample,y=prportion,fill = Cell_Type)) +geom_col() +
  scale_fill_manual(values = c(Classical = "purple", Basal = "pink")) +
  xlab("Sample_ID")+
  ylab("Proportion")+my_theme+theme(axis.text.x = element_text(angle = 90))

png("plots/tumor_type_proportion.png", res = 300, width = 300*5, height = 300*5)
p
dev.off()

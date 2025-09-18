setwd("D:/Anya/pdac")


library(Seurat)
library(rjags)
library(infercnv)
library(dplyr)
library(readr)
library(dplyr)

all_samples <- readRDS("D:/Anya/pdac/all_samples_no_integration_layer.rds")
sample_list <- SplitObject(all_samples, split.by = "orig.ident")

#set the one object p03 -----------------------
object <- sample_list$P03


annotation_df <- data.frame(
  cell_name  = colnames(object),
  cell_type  = as.character(object@active.ident),
  stringsAsFactors = FALSE
)

# define obs/ref
obs_types  <- c("Epithelial", "Prolif_Epithelial")          
ref_types  <- c("CAFs","TCell/NK","Prolif_Lymphoid","Myeloid")  
ref_mask <- annotation_df$group %in% ref_types
set.seed(1)
if (sum(ref_mask) > 1000) {
  keep_ref <- sample(which(ref_mask), 1000)
  keep_idx <- c(keep_ref, which(!ref_mask))
  annotation_df <- annotation_df[sort(keep_idx), , drop = FALSE]
}



write.table(
  annotation_df,
  file = "full_annotation.tsv",
  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

#I did not use endo as references 
keep_cells <- annotation_df$cell_name
object_count_sub <- GetAssayData(object, slot = "counts")[, keep_cells, drop = FALSE]
stopifnot(all(colnames(object_count_sub) == keep_cells))

# create inferCNV object
ref_group_names <- ref_types 
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = as.matrix(object_count_sub),
  annotations_file  = "full_annotation.tsv",
  delim             = "\t",
  gene_order_file   = "hg38_gencode_v27.tsv",
  ref_group_names   = ref_group_names,

)


infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,                 
  out_dir           = "P03",
  cluster_by_groups = TRUE,               
  denoise           = TRUE,
  HMM               = TRUE,
  analysis_mode     = "subclusters",
  leiden_resolution = 0.01,
  up_to_step        = 15                   
)


# which one had state more than 1 based on hmm obdervation
pred <- read_tsv("P03/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat",
                 show_col_types = FALSE)

grp_path <- file.path(out_dir, "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt")
grp_raw <- read.table(grp_path, sep=" ")
grp_raw <- grp_raw[,1:5]
grp_raw$barcode=rownames(grp_raw)
colnames(grp_raw)
grp <- grp_raw[,c("barcode","Dendrogram.Group")]

colnames(grp)<- c("cell_name","cluster")
sid <- function(x) sub(".*(s\\d+).*", "\\1", x)  

pred$cluster_key <- sid(pred$cell_group_name)
grp$cluster_key  <- sid(grp$cluster)
grp <- grp[,c("cell_name","cluster_key")]
mal_clusters <- pred %>%
  filter(state > 1) %>%             
  pull(cluster_key) %>%
  unique()


# 3) attach to Seurat object
object$infercnv_hmm_cluster <- grp$cluster[ match(colnames(object), grp$cell_name) ]
object$infercnv_hmm_call <- ifelse(object$infercnv_hmm_cluster %in% mal_clusters,
                                   "malignant","normal")

table(object$infercnv_hmm_call, useNA="ifany")


p03 <- DimPlot(object, group.by = "infercnv_hmm_call")
png ("hmm_pred.png", width = 300*6, height = 300*6, res = 300)
p03
dev.off()
#create a function that automates this for each sample ####
run_inferCNV_hmm <- function(object, sample_id, gene_loc_file,
                             out_root="infercnv_runs",
                             obs_types=c("Epithelial","Prolif_Epithelial"),
                             ref_types=c("CAFs","TCell/NK","Prolif_Lymphoid","Myeloid"),
                             max_refs=1000, cutoff=0.1) {
  
  dir_out <- file.path(out_root, sample_id)
  dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
  
  ann <- data.frame(
    cell_name = colnames(object),
    cell_type = as.character(object@active.ident),
    stringsAsFactors = FALSE
  )
  ref_counts <- table(ann$cell_type[ann$cell_type %in% ref_types])
  keep_ref_types <- names(ref_counts[ref_counts >= 2])
  
  ann <- ann[ann$cell_type %in% c(obs_types, keep_ref_types), , drop=FALSE]
  ref_group_names <- keep_ref_types
  if (length(ref_group_names) == 0) stop("No reference groups with ≥2 cells; supply pooled refs.")
  # keep all obs, cap refs
  is_obs <- ann$cell_type %in% obs_types
  is_ref <- ann$cell_type %in% ref_group_names
  set.seed(1)
  if (sum(is_ref) > max_refs) {
    keep_ref <- sample(which(is_ref), max_refs)
    keep_idx <- sort(c(which(is_obs), keep_ref))
    ann <- ann[keep_idx, , drop=FALSE]
  }
  
  # write inferCNV annotation (2 cols)
  ann_path <- file.path(dir_out, "annotation.tsv")
  write.table(ann[,c("cell_name","cell_type")], ann_path,
              sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # subset counts
  keep_cells <- ann$cell_name
  m <- GetAssayData(object, slot="counts")[, keep_cells, drop=FALSE]
  stopifnot(all(colnames(m) == keep_cells))
  
  # run
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = as.matrix(m),
    annotations_file  = ann_path,
    delim             = "\t",
    gene_order_file   = gene_loc_file,
    ref_group_names   = keep_ref_types
  )
  
  infercnv::run(
    obj,
    cutoff            = cutoff,
    out_dir           = dir_out,
    cluster_by_groups = TRUE,
    analysis_mode     = "subclusters",
    denoise           = TRUE,
    HMM               = TRUE,
    leiden_resolution = 0.01,
    up_to_step        = 17
  )
  
  return(dir_out)  # just return the output folder path
}




#debuggging--------------
object <- sample_list$P13
ref_types=c("CAFs","TCell/NK","Prolif_Lymphoid","Myeloid")

ann <- data.frame(
  cell_name = colnames(object),
  cell_type = as.character(object@active.ident),
  stringsAsFactors = FALSE
)
ref_counts <- table(ann$cell_type[ann$cell_type %in% ref_types])
keep_ref_types <- names(ref_counts[ref_counts >= 2])

ann <- ann[ann$cell_type %in% c(obs_types, keep_ref_types), , drop=FALSE]
ref_group_names <- keep_ref_types
if (length(ref_group_names) == 0) stop("No reference groups with ≥2 cells; supply pooled refs.")
# keep all obs, cap refs
is_obs <- ann$cell_type %in% obs_types
is_ref <- ann$cell_type %in% ref_group_names
set.seed(1)
if (sum(is_ref) > max_refs) {
  keep_ref <- sample(which(is_ref), max_refs)
  keep_idx <- sort(c(which(is_obs), keep_ref))
  ann <- ann[keep_idx, , drop=FALSE]
}

print(ref_types==unique(ann$cell_type))

print(unique(ann$cell_type))

#prof_lymphoid does not exisit in the anno of P13

print(keep_ref_types) #use keep refrence type




dir_out<- "infercnv_runs"
#attach mata data to each object--------------------
attach_infercnv_hmm_to_seurat <- function(object, dir_out) {
  library(readr); library(dplyr)
  
  # ---- find files inside this sample folder ----
  pred_path <- file.path(dir_out, pattern = "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat")
  grp_path  <- file.path(dir_out, pattern = "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt")
  
  if (length(pred_path) == 0 || length(grp_path) == 0) {
    stop(print(paste0("Missing inferCNV outputs",
                 dir_out, length(pred_path), length(grp_path))))
  }

  
  # ---- read predictions & define malignant clusters ----
  pred <- read_tsv(pred_path, show_col_types = FALSE)
  
  key <- function(x) sub(".*(s\\d+).*", "\\1", x)   # fix naming issue
  pred$cluster_key <- key(pred$cell_group_name)
  
  mal_keys <- pred %>% filter(state > 1) %>% pull(cluster_key) %>% unique()

  # ---- read groupings ----
  grp_raw <- read.table(grp_path, sep = "")
  grp_raw$barcode=rownames(grp_raw)
 
  grp_raw <- grp_raw[,c("barcode","Dendrogram.Group")]
  
  colnames(grp_raw)<- c("cell_name","cluster")
  grp_raw$cluster_key <- key(grp_raw$cluster)

  # ---- write into metadata ----
  object$infercnv_hmm_cluster <- grp_raw$cluster_key[ match(colnames(object), grp_raw$cell_name) ]
  object$infercnv_hmm_call <- ifelse(object$infercnv_hmm_cluster %in% mal_clusters,
                                     "malignant","normal")
  
  print(table(object$infercnv_hmm_call, useNA="ifany"))
  
  return(object)
}

for (sid in names(sample_list)) {
  sample_list[[sid]] <- attach_infercnv_hmm_to_seurat(sample_list[[sid]], dir_out = dirs[[sid]])
}

#function to rbind all meta datas
merged_df <- bind_rows(
  lapply(names(sample_list), function(sid) {
    md <- sample_list[[sid]]@meta.data
    md$.cell   <- rownames(md)      
    md$.sample <- sid                
    md
  })
)

rownames(merged_df) <- merged_df$.cell

main = all_samples
# ensure main cell names match merged_df rownames
stopifnot(colnames(all_samples)==rownames(merged_df))
print(colnames(all_samples)==rownames(merged_df))
# drop columns that already exist in main@meta.data to avoid clashes
new_cols <- setdiff(colnames(merged_df), colnames(main@meta.data))
to_add  <- merged_df[, new_cols, drop = FALSE]

main@meta.data <- merged_df
unique(main@active.ident)


DimPlot(main, group.by = "infercnv_hmm_call")


#subset the epithelial cells-----
sc.subset <- subset(main, subset =cell.types %in% c("Epithelial", "Prolif_Epithelial"))
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
p<- DimPlot(sc.subset, group.by = "infercnv_hmm_call", pt.size = 0.5)
p1<- DimPlot(sc.subset, group.by = "orig.ident", pt.size = 0.5)

png("malignant.png", res = 300, width = 300*10, height = 300*5)
p1+p
dev.off()
epi <-sc.subset
saveRDS(main, "all_seurat_complete_infercnv.rds")

#subset Tumor micro environment ----------------

main@active.ident
sc.subset <- subset(main, subset =cell.types %in% c("TCell/NK" , "Myeloid", "CAFs" , "Endothelial", "Prolif_Lymphoid", "Prolif_Myeloid"))
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
p<- DimPlot(sc.subset, group.by = "cell.types", pt.size = 0.5)
p1<- DimPlot(sc.subset, group.by = "orig.ident", pt.size = 0.5)

png("non_epi.png", res = 300, width = 300*10, height = 300*5)
p1+p
dev.off()

saveRDS(main, "all_seurat_complete_infercnv.rds")

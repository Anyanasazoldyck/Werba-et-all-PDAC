setwd("D:/Anya/pdac")


library(Seurat)
library(rjags)
library(infercnv)
library(dplyr)
library(readr)
library(dplyr)

all_samples <- readRDS("D:/Anya/pdac/all_samples_no_integration_layer.rds")
sample_list <- SplitObject(all_samples, split.by = "orig.ident")


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
  )
  
  return(dir_out)  # just return the output folder path
}








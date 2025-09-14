###############################################################################################
#------------------identifying malingant epithilal cells using inferCNV-------------------------
###############################################################################################
setwd("D:/PDAC_SCRNA")
#call INFERCNV and rJAGS####
library(infercnv)
library(rjags)
library(Seurat)
library(dplyr)

#split my seurat object into samples
samples_list <- SplitObject(all_samples_integrated_labeled, split.by = "orig.ident") #splitting by sample 

#download the genes/chromosome mapping file ####
download.file(
  url = "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt",
  destfile = "csv/hg38_gencode_v27.tsv"
)



# inferCNV sample per sample
P03 <- samples_list$P03
#count matrix
anno_df <- data.frame(cell_id = colnames(P03),
                      group   = as.character(Idents(P03)))



#refrence groups 
epi_patterns <- c("Epithelial", "Prolif_Epithelial") 
is_epithelial <- function(x) grepl(paste(epi_patterns, collapse="|"), x, ignore.case = TRUE)

epi_cells <- anno_df$cell_id[ is_epithelial(anno_df$group) ]
ref_cells_all <- anno_df$cell_id[ !is_epithelial(anno_df$group) ]

#get rid of cells above 1000
set.seed(1)
ref_cap <- 1000
ref_cells <- if (length(ref_cells_all) > ref_cap) sample(ref_cells_all, ref_cap) 
keep_cells <- union(epi_cells, ref_cells)
anno_df <- anno_df[match(keep_cells, anno_df$cell_id), , drop = FALSE]

#save annotation
write.table(anno_df, file = "csv/P03_annotations.tsv",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
ref_groups <- setdiff(unique(anno_df$group), c("NormalSpike", unique(anno_df$group[is_epithelial(anno_df$group)])))

#checking if there are duplicates                   
stopifnot(all(anno_df$cell_id %in% colnames(P03)))
stopifnot(!anyDuplicated(anno_df$cell_id))


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(P03,layer = "count"),
                                    annotations_file="csv/P03_annotations.tsv",
                                    delim="\t",
                                    gene_order_file="csv/hg38_gencode_v27.tsv",
                                    ref_group_names=ref_groups) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="P03_example", 
                             cluster_by_groups=TRUE, 
                             up_to_step = 10,
                             leiden_resolution = 0.01, 
                             HMM=TRUE,
                             plot_steps = T)


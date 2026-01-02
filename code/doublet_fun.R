

rm(list=ls())
gc()



#---------------doublet finder per sample------------------------

run_doubletfinder_lognorm <- function(
    seu,
    elbow_drop = 0.5,          
    target_cumvar = 90,        
    min_pc_floor = 5,         
    resolution = 0.1)
{
  # Basic preprocessing
  DefaultAssay(seu) <- "RNA"
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  
  # ---- chose PC "Squid Stat" ----
  stdev <- seu@reductions$pca@stdev
  percent_var <- (stdev^2 / sum(stdev^2)) * 100 #elbow plot values
  cumulative_var <- cumsum(percent_var)
  PC1 <- which(cumulative_var > target_cumvar)[1]
  dvar <- diff(percent_var)
  PC2 <- which(dvar < elbow_drop)[1] + 1
  min_pc <- max(min_pc_floor, min(PC1, PC2))  
  min_pc= as.numeric(min_pc)
  pcs_use <- 1:min_pc
  #find seurat clusters  
  seu <- FindNeighbors(seu, dims = pcs_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  
  # PK
  sweep.res <- paramSweep(seu, PCs = pcs_use, sct = F) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
  best.pK <- as.numeric(as.vector(find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]))
  
  
  
  # ---- Estimate multiplet rate from 10x lookup table ----
  print(">> Estimating multiplet rate from 10x table/mnt/d/Command Line Tutorials/raw_data/data_d4.")
  multiplet_rates_10x <- data.frame( 'Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                     'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000), 
                                     'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000) ) 
  multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu@meta.data)) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% 
    dplyr::select(Multiplet_rate) %>% 
    as.numeric(as.character())
  # ---- Homotypic adjustment ----
  annotations <- seu@meta.data$seurat_clusters
  # use the clusters as the user-defined cell types 
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(multiplet_rate * nrow(seu@meta.data)) # multiply by number of cells to get the number of expected multiplets 
  nExp.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  
  # ---- Run DoubletFinder ----
  print(">> Running DoubletFinder/mnt/d/Command Line Tutorials/raw_data/data_d4.")
  seu <- DoubletFinder::doubletFinder(
    seu = seu,
    PCs = pcs_use, 
    pN = 0.25, 
    pK = best.pK,
    nExp = nExp.adj, 
    sct = FALSE
  )
  
  
  
  df.col <- grep("^DF.classification", colnames(seu@meta.data), value = TRUE)[1]
  seu <- subset(seu, cells = colnames(seu)[seu@meta.data[[df.col]] == "Singlet"])
  
  
  # Return results + parameters for audit
  list(seu)
}



res_list <- lapply(names(samp_split), function(nm) {
  out <- run_doubletfinder_lognorm(samp_split[[nm]])
  out
})
names(res_list) <- names(samp_split)
saveRDS(res_list,"res_doublet_filtered.rds")

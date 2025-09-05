

rm(list=ls())
gc()



#---------------doublet finder per sample------------------------
#set a file dir
doublets= file.path(plots,"doublets")
dir.create(doublets)
library(DoubletFinder)
library(sctransform)
library(glmGamPoi)
library(Seurat)
process_sample <- function(sample) {
  # Basic preprocessing
  sample <- SCTransform(sample, verbose = FALSE)
  sample <- RunPCA(sample, verbose = FALSE)
  sample <- RunUMAP(sample, dims = 1:30)
  sample <- FindNeighbors(sample, dims = 1:30)
  sample <- FindClusters(sample, resolution = 0.2)
  
  
  #find the perfect p dim 
  percent_var <- (stdv^2/sum(stdv^2)) * 100
  cumulative_var <- cumsum(percent_var)
  PC1 <- which(cumulative_var > 90)[1]
  PC2 <- which(diff(percent_var) < 0.1)[1] + 1
  min_pc <- min(PC1, PC2)
  
  # Run parameter sweep to get best pK
  sweep.res <- paramSweep_v3(sample, PCs = 1:min_pc, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  best.pK <- as.numeric(as.vector(find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]))
  
  # Expected doublet rate (rule of thumb: ~1% per 1000 cells)
  multiplet_rates_10x <- data.frame(
    'Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    'Loaded_cells'  = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  
  multiplet_rate <- multiplet_rates_10x %>%
    dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>%
    dplyr::slice(which.max(Recovered_cells)) %>%
    dplyr::select(Multiplet_rate) %>%
    as.numeric(as.character())
  
  
  # Homotypic adjustment
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  # Run DoubletFinder
  sample <- doubletFinder_v3(sample,
                             PCs = 1:min_pc,
                             pN = 0.25,
                             pK = best.pK,
                             nExp = nExp.adj,
                             reuse.pANN = FALSE,
                             sct = TRUE)
  
  # Keep singlets only
  df.col <- grep("DF.classification", colnames(sample@meta.data), value = TRUE)[1]
  sample <- subset(sample, cells = rownames(sample@meta.data)[sample@meta.data[[df.col]] == "Singlet"])
  
  return(sample)
}





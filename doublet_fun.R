#---------------doublet finder per sample------------------------
#set a file dir
doublets= file.path(plots,"doublets")
dir.create(doublets)
library(DoubletFinder)


process_sample <- function(sample) {
  # Basic preprocessing
  sample <- SCTransform(sample, verbose = FALSE)
  sample <- RunPCA(sample, verbose = FALSE)
  sample <- RunUMAP(sample, dims = 1:30)
  sample <- FindNeighbors(sample, dims = 1:30)
  sample <- FindClusters(sample, resolution = 0.2)
  
  # Run parameter sweep to get best pK
  sweep.res <- paramSweep_v3(sample, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  best.pK <- as.numeric(as.vector(find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]))
  
  # Expected doublet rate (rule of thumb: ~1% per 1000 cells)
  nCells <- ncol(sample)
  exp.dbl.rate <- 0.01 * (nCells/1000)
  nExp <- round(exp.dbl.rate * nCells)
  
  # Homotypic adjustment
  homotypic.prop <- modelHomotypic(sample$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  # Run DoubletFinder
  sample <- doubletFinder_v3(sample,
                             PCs = 1:30,
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

# Apply to all samples
samples_clean <- lapply(seurat_list_no_liver, process_sample)


# Now merge for integration
all_samples <- merge(samples_clean[[1]], y = samples_clean[-1])



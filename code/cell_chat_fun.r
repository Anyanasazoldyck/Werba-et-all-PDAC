do_cell_chat <- function(object, group, assay = "RNA", min_cells = 10) {
  stopifnot(group %in% colnames(object@meta.data))
  
  # load DB
  if (!exists("CellChatDB.human")) data(CellChatDB.human, package = "CellChat")
  cellchat_DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")  # or "ECM-Receptor"/"Cell-Cell Contact"
  
  # build CellChat from Seurat
  cellchat <- createCellChat(object = object, group.by = group, assay = assay)
  cellchat@DB <- cellchat_DB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  
  Interactions <- subsetCommunication(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  list(cellchat_data_object = cellchat, interaction = Interactions)
}



# load package
library(CellChat)
library(NMF)
library(ggalluvial)
library(Seurat)


#create file ####
dir.create("cell_chat", recursive = T)

## prepare the 2 seurat object for each condition 
all_samples <- readRDS("D:/PDAC_SCRNA/objects/all_seurat_complete_infercnv.rds")

treated <- subset(all_samples, subset=Treatment=="treated")
naive <- subset(all_samples, subset=Treatment=="naive")


print(table(treated$cell.types))
print(table(naive$cell.types))


#group into MAJOR subtypes
major_map <- c(
  "TCell/NK"="Lymphoid", "Prolif_Lymphoid"="Lymphoid",
  "Myeloid"="Myeloid", "Prolif_Myeloid"="Myeloid",
  "Epithelial"="Epithelial", "Prolif_Epithelial"="Epithelial",
  "CAFs"="CAF", "Endothelial"="Endothelial"
)


treated$cell_major <- unname(major_map[as.character(treated$cell.types)])
treated$cell_major <- factor(treated$cell_major)

naive$cell_major <- unname(major_map[as.character(naive$cell.types)])
naive$cell_major <- factor(naive$cell_major)


treated_ch = do_cell_chat(object = treated,assay = "RNA",group = "cell_major",min_cells = 10)
naive_ch = do_cell_chat(object = naive,assay = "RNA",group = "cell_major",min_cells = 10)


#####################################################
#Visualize each condition #####
####################################################
#1- naive
cell_chate_naive = naive_ch$cellchat_data_object
groupSize = as.numeric(table(cell_chate_naive@idents)) 

# Plots a circos plot to show the number of interactions between cell groups
naive_network<-netVisual_circle(cell_chate_naive@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
png("cell_chat/naive_network.png", res = 300, height = 300*6, width = 300*6)
naive_network
dev.off()
# Plots a circos plot to show the strength of interactions between cell groups
naive_network2<-netVisual_circle(cell_chate_naive@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
png("cell_chat/naive_network_weight.png", res = 300, height = 300*6, width = 300*6)
naive_network2
dev.off()

png("cell_chat/naive_network_seperated.png", res = 300, height = 300*6, width = 300*12)
mat = cell_chate_naive@net$weight
par(mfrow = c(2,6), xpd=TRUE, mar = c(0.75,0.75,0.75,0.75))
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
naive_network3=par(mfrow=c(1,1))
dev.off()


#2- treated
cell_chate_treated = treated_ch$cellchat_data_object
groupSize = as.numeric(table(cell_chate_treated@idents)) 

# Plots a circos plot to show the number of interactions between cell groups
treated_network<-netVisual_circle(cell_chate_treated@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
png("cell_chat/treated_network.png", res = 300, height = 300*6, width = 300*6)
treated_network
dev.off()
# Plots a circos plot to show the strength of interactions between cell groups
treated_network2<-netVisual_circle(cell_chate_treated@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
png("cell_chat/treated_network_weight.png", res = 300, height = 300*6, width = 300*6)
treated_network2
dev.off()

png("cell_chat/treated_network_seperated.png", res = 300, height = 300*6, width = 300*12)
mat = cell_chate_treated@net$weight
par(mfrow = c(2,6), xpd=TRUE, mar = c(0.75,0.75,0.75,0.75))
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
treated_network3=par(mfrow=c(1,1))
dev.off()


#---------------------------------
# make a heatmap of interactions-------
#it seems that the chemotherapy decresed the cell-cell interaction post between myeloid/lymphoid and cancer cells.
#---------------------------------
p=netVisual_heatmap(cell_chate_naive, measure = "weight")
png("cell_chat/naive_hm.png", res = 300, height = 300*6, width = 300*6)
p
dev.off()


treated_heatmap=netVisual_heatmap(cell_chate_treated, measure = "weight")
png("cell_chat/treated_hm.png", res = 300, height = 300*6, width = 300*6)
treated_heatmap
dev.off()

#-----------------------------------------
#which pathways got activated /canceled by chemo-----
#--------------------------------------
naive <- cell_chate_naive@netP$pathways
treated <- cell_chate_treated@netP$pathways


diff_path <- setdiff(treated, naive)


#------------------------------------
#differental analysis of the interaction------
#-------------------------------------
cellchat.list = list(Control = cell_chate_naive, Treatement = cell_chate_treated) 

merged_cellchat = mergeCellChat(cellchat.list, add.names = names(cellchat.list))



# differential interaction heatmaps
g1 = netVisual_heatmap(merged_cellchat)
g2 = netVisual_heatmap(merged_cellchat , measure = "weight") 
g1 + g2
png("cell_chat/differential interactions_hm.png", res = 300, height = 300*6, width = 300*6)
g1 + g2
dev.off()

# plot comaprative in and out strengths
png("cell_chat/differential interactions_in_out.png", res = 300, height = 300*6, width = 300*6)

num.link = sapply(cellchat.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax = c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg = list()
for (i in 1:length(cellchat.list)) {
  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.list[[i]], title = names(cellchat.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()

# ranking them  by pathway
rankNet(merged_cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
png("cell_chat/differential interactions_pathway.png", res = 300, height = 300*6, width = 300*6)
rankNet(merged_cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
dev.off()

# Differential , how CAF interacts with each compartment in each condition.
png("cell_chat/differential interactions_diff.png", res = 300, height = 300*10, width = 300*6)
netVisual_bubble(merged_cellchat, sources.use = "CAF", targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45)
dev.off()



#-------------------------------------------
#-------what is popping up in treatment and how is it interacting -------
#--------------------------------------


# plot the sending and recieving scores



p=netAnalysis_signalingRole_heatmap(cell_chate_treated, pattern = "outgoing")
p1=netAnalysis_signalingRole_heatmap(cell_chate_treated, pattern = "incoming")
png("cell_chat/treated_system_level.png", res = 300, height = 300*5, width = 300*10)
p+p1
dev.off()
#it seems like CAF and Myeloid are talking through OSM pathway
netAnalysis_signalingRole_scatter(cell_chate_treated, signaling = c("OSM"))

# identify enriched pairs in osm
pairLR = extractEnrichedLR(cell_chate_treated, signaling = "OSM", geneLR.return = FALSE)
pairLR

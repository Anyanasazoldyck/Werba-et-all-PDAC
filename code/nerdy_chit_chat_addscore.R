setwd("D:/PDAC_SCRNA")


#libraries##############





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


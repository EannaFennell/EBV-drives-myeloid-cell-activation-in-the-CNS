library(Seurat)
library(ggplot2)
library(scRepertoire)
library(ggpubr)
library(viridis)

spleen_brain <- readRDS("D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint4.rds")
DimPlot(spleen_brain, reduction = "umap.rpca", label = TRUE)
table(spleen_brain@meta.data[["cloneType"]])

############
# Subset T cell subsets
############

Idents(spleen_brain) <- spleen_brain$predicted.id
spleen_brain_t_cells <- subset(spleen_brain, idents = c('Regulatory T cells',
                                                        "Tcm/Naive helper T cells",
                                                        "Tem/Trm cytotoxic T cells",
                                                        "Tem/Effector helper T cells",
                                                        "Tem/Temra cytotoxic T cells",
                                                        "MAIT cells",
                                                        "CD16- NK cells",
                                                        "CD16+ NK cells",
                                                        "Follicular helper T cells",
                                                        "CRTAM+ gamma-delta T cells",
                                                        "Tcm/Naive cytotoxic T cells",
                                                        'Cycling immune mix'))

DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE)

############
# Drop BCR data
############

drops <- c("CTgene","CTnt","CTaa","CTstrict","Frequency","cloneType")
spleen_brain_t_cells@meta.data <- spleen_brain_t_cells@meta.data[ , !(names(spleen_brain_t_cells@meta.data) %in% drops)]

############
# Add TCRs
############

# Spleen

spleen_tcrs <- read.csv("D:/Tbet_B_cell/scRNAseq_Exp2/Spleen_Plate/per_sample_outs/Spleen_Plate/vdj_t/filtered_contig_annotations.csv")

combined.TCR.spleen <- combineTCR(spleen_tcrs, 
                                  samples = "spleen")

head(combined.TCR.spleen[[1]])

spleen_brain_t_cells <- combineExpression(
  combined.TCR.spleen,
  spleen_brain_t_cells,
  cloneCall = "gene",
  cloneTypes = c(None = 0, Single = 1, Small = 5, Medium = 20, Large = 50,
                 Hyperexpanded = 200),
  proportion = FALSE
)

# Brain

brain_tcrs <- read.csv("D:/Tbet_B_cell/scRNAseq_Exp2/Brain_Plate/per_sample_outs/Brain_Plate/vdj_t/filtered_contig_annotations.csv")

combined.TCR.brain <- combineTCR(brain_tcrs, 
                                 samples = "brain")

head(combined.TCR.brain[[1]])

spleen_brain_t_cells <- combineExpression(
  combined.TCR.brain,
  spleen_brain_t_cells,
  cloneCall = "gene",
  cloneTypes = c(None = 0, Single = 1, Small = 5, Medium = 20, Large = 50,
                 Hyperexpanded = 200),
  proportion = FALSE
)

table(spleen_brain_t_cells@meta.data[["cloneType"]])

############
# Remove cells without TCR
############

Idents(spleen_brain_t_cells) <- "CTgene"

length(unique(spleen_brain_t_cells@meta.data$CTgene))
unique(spleen_brain_t_cells@meta.data$CTgene)[1]
length(unique(spleen_brain_t_cells@meta.data$CTgene))
keep <- unique(spleen_brain_t_cells@meta.data$CTgene)[c(2:2779)]

spleen_brain_t_cells <- subset(spleen_brain_t_cells, idents = keep)

Idents(spleen_brain_t_cells) <- spleen_brain_t_cells$predicted.id

DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE)

############
# Re-cluster (silence the TCR genes)
############

spleen_brain_t_cells <- PercentageFeatureSet(spleen_brain_t_cells, pattern = "^MT-", col.name = "percent.mt")
spleen_brain_t_cells <- SCTransform(spleen_brain_t_cells, vars.to.regress = "percent.mt", verbose = TRUE)

VariableFeatures(spleen_brain_t_cells)[1:20]
Seurat::VariableFeatures(spleen_brain_t_cells) <- quietVDJgenes(Seurat::VariableFeatures(spleen_brain_t_cells))
VariableFeatures(spleen_brain_t_cells)[1:20]

spleen_brain_t_cells <- FindNeighbors(spleen_brain_t_cells, reduction = "integrated.rpca", dims = 1:30)
spleen_brain_t_cells <- FindClusters(spleen_brain_t_cells, resolution = 0.5, cluster.name = "rpca_clusters")
spleen_brain_t_cells <- RunUMAP(spleen_brain_t_cells, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")

Idents(object = spleen_brain_t_cells) <- "organ"
Idents(object = spleen_brain_t_cells) <- "seurat_clusters"
Idents(object = spleen_brain_t_cells) <- "ebv_pbs"
Idents(object = spleen_brain_t_cells) <- "predicted.id"
Idents(object = spleen_brain_t_cells) <- "t_cell_states"
Idents(object = spleen_brain_t_cells) <- "cd4_cd8"

DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE, pt.size = 1.5) + xlab("UMAP 1") + ylab("UMAP 2")

DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = FALSE, pt.size = 1.5, 
        split.by = "ebv_pbs") + xlab("UMAP 1") + ylab("UMAP 2")


col.clusters <- c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88","black","pink")
col.clusters <- c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88")


DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE, pt.size = 1.5, 
        label.size = 5, split.by = "organ") +
  scale_shape_manual(values = 21) +
  scale_colour_manual(values = col.clusters) + 
  theme(legend.position = "none") + xlab("UMAP 1") + 
  ylab("UMAP 2")

p <- DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE, 
             pt.size = 1.5, 
             label.size = 4)

p +
  geom_point(
    data = p$data,
    aes(umaprpca_1, umaprpca_2, fill = ident),   # use the same grouping column name you used in group.by
    shape = 21,         # filled circle that supports stroke
    color = "black",    # border color
    stroke = 0.15,      # border thickness — tune this
    size = 2          # point size — tune to taste
  ) + scale_fill_manual(values = col.clusters)


spleen_brain_t_cells_for_anndata <- spleen_brain_t_cells

spleen_brain_t_cells_for_anndata@meta.data[] <- lapply(spleen_brain_t_cells_for_anndata@meta.data, function(col) {
  if (is.factor(col)) as.character(col) else col
})

seurat_to_anndata(spleen_brain_t_cells_for_anndata, h5path = "D:/Tbet_B_cell/scRNAseq_Exp2/anndata_files/t_cell_umap_test_3")

counts_matrix <- GetAssayData(spleen_brain_t_cells_for_anndata, assay = "RNA", slot = "counts")
write.csv(as.matrix(counts_matrix), "D:/Tbet_B_cell/scRNAseq_Exp2/anndata_files/raw_matrices/counts.csv")
write.csv(spleen_brain_t_cells_for_anndata@meta.data, "D:/Tbet_B_cell/scRNAseq_Exp2/anndata_files/raw_matrices/metadata.csv")
write.csv(spleen_brain_t_cells_for_anndata@reductions$umap.rpca@cell.embeddings, "D:/Tbet_B_cell/scRNAseq_Exp2/anndata_files/raw_matrices/umap_embeddings.csv")


# rds file checkpoint 1 - just t cells
saveRDS(spleen_brain_t_cells,"D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint6_just_t_cells.rds")
spleen_brain_t_cells <- readRDS("D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint6_just_t_cells.rds")

############
# Correct animal IDs
############

spleen_brain_t_cells@meta.data$ebv_pbs <- spleen_brain_t_cells@meta.data[["hash.ID"]]

spleen_brain_t_cells@meta.data$ebv_pbs[spleen_brain_t_cells@meta.data$ebv_pbs %in% c("Hash9","Hash12",
                                                                                     "Hash14","Hash17",
                                                                                     "Hash2","Hash21",
                                                                                     "Hash23","Hash4",
                                                                                     "Hash13","Hash8")] <- "EBV"

spleen_brain_t_cells@meta.data$ebv_pbs[spleen_brain_t_cells@meta.data$ebv_pbs %in% c("Hash1",
                                                                                     "Hash6","Hash10")] <- "PBS"

table(spleen_brain_t_cells@meta.data$ebv_pbs)


spleen_brain_t_cells@meta.data$mouse_id <- spleen_brain_t_cells@meta.data[["hash.ID"]]

spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash2")] <- "M22"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash4")] <- "M26"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash6")] <- "M33"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash8")] <- "M7"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash10")] <- "M8"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash12")] <- "M13"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash14")] <- "M16"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash17")] <- "M5"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash21")] <- "M6"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash23")] <- "M21"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash1")] <- "M11"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash9")] <- "M35"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash13")] <- "M3"

table(spleen_brain_t_cells@meta.data$mouse_id)


# rds file checkpoint
saveRDS(spleen_brain_t_cells,"D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint8_just_t_cells.rds")
spleen_brain_t_cells <- readRDS("D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint8_just_t_cells.rds")


############
# 
############

Idents(object = spleen_brain_t_cells) <- "cloneType"
Idents(object = spleen_brain_t_cells) <- "seurat_clusters"

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
colorblind_vector <- colorblind_vector[c(1,3,4,5,7)]

DimPlot(spleen_brain_t_cells, label = FALSE, pt.size = 1 , reduction = "umap.rpca", split.by = "ebv_pbs") +
  scale_color_manual(values = colorblind_vector) + 
  ylab("UMAP 2") + xlab("UMAP 1")


FeaturePlot(spleen_brain_t_cells, features = "CCR5", 
            order = T, max.cutoff = 'q97', 
            reduction = "umap.rpca", pt.size = 1) + xlab("UMAP 1") + ylab("UMAP 2") & 
  scale_color_viridis(option = "D") 

FeaturePlot(spleen_brain_t_cells, features = "CD8A", 
            order = T, max.cutoff = 'q97', 
            reduction = "umap.rpca", pt.size = 1) + xlab("UMAP 1") + ylab("UMAP 2") & 
  scale_color_viridis(option = "D") 

FeaturePlot(spleen_brain_t_cells, features = "CD69", 
            order = T, max.cutoff = 'q97', split.by = "organ",
            reduction = "umap.rpca", pt.size = 1) + xlab("UMAP 1") + ylab("UMAP 2") & 
  scale_color_viridis(option = "D") 




FeaturePlot(spleen_brain_t_cells, features = "IFNG", split.by = "ebv_pbs", 
            order = T, max.cutoff = 'q97', 
            reduction = "umap.rpca", pt.size = 1) & 
  scale_color_viridis(option = "D")

FeaturePlot(spleen_brain_t_cells, features = "IL1B", split.by = "organ", 
            order = T, max.cutoff = 'q97', 
            reduction = "umap.rpca", pt.size = 1) & 
  scale_color_viridis(option = "D")


chemotaxis.genes <- c("CCR7","CXCR4","CCR5","CXCR3","CCR6","CCR4",
                      "CCR8","CCR9","CCR10","CXCR5","CX3CR1","CCR2",
                      "CXCR6")



DotPlot(spleen_brain_t_cells, features = c("CD4","CD8A","FOXP3","IL2RA",
                                           "SELL","CCR7","CXCR5","ICOS",
                                           "GZMB","PRF1","IFNG","HLA-DRA",
                                           "CD69","ITGAE","MKI67","CCNA2",
                                           "CCNB1"))  + rotate_x_text(45) & 
  scale_color_viridis(option = "D")


spleen_brain_t_cells <- PrepSCTFindMarkers(spleen_brain_t_cells, assay = "SCT", verbose = TRUE)
all.markers <- FindAllMarkers(object = spleen_brain_t_cells)

all.markers.pos <- all.markers[all.markers$avg_log2FC > 0,]
subset <- all.markers.pos[all.markers.pos$cluster == 9,]


VlnPlot(spleen_brain_t_cells, features = c("CXCR6"), split.by = "organ")

write.csv(all.markers, "D:/Tbet_B_cell/scRNAseq_Exp2/dge/findallmarkers_tcells.csv")



annotation <- c("Tcm Helper","Tem Cytotoxic","Tcm Cytotoxic","Tem Cytotoxic","Tem Helper","Teff Cytotoxic","Teff Cytotoxic","Naive","Teff Helper","Tregs")
names(annotation) <- levels(spleen_brain_t_cells)
spleen_brain_t_cells <- RenameIdents(spleen_brain_t_cells, annotation)
spleen_brain_t_cells@meta.data$t_cell_states <- spleen_brain_t_cells@active.ident

Idents(object = spleen_brain_t_cells) <- "seurat_clusters"
annotation <- c("CD4","CD8","CD8","CD8","CD4","CD8","CD8","Naive","CD4","CD4")
names(annotation) <- levels(spleen_brain_t_cells)
spleen_brain_t_cells <- RenameIdents(spleen_brain_t_cells, annotation)
spleen_brain_t_cells@meta.data$cd4_cd8 <- spleen_brain_t_cells@active.ident


DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", label = TRUE, label.size = 6,
        pt.size = 1.5) + 
  scale_color_manual(values = c("#C11C84","#2BDD88","black")) + theme(legend.position = "none") + xlab("UMAP 1") + ylab("UMAP 2")

############
# Cell Abundances
############

library(tidyr)

Idents(object = spleen_brain_t_cells) <- "t_cell_states"

cell_abun <- as.data.frame(table(spleen_brain_t_cells@active.ident, spleen_brain_t_cells@meta.data$organ))

ggplot(cell_abun, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88")) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))


Idents(object = spleen_brain_t_cells) <- "cloneType"


spleen_brain_t_cells@active.ident <- factor(spleen_brain_t_cells@active.ident, levels = rev(c("Single (0 < X <= 1)",
                                                                                              "Small (1 < X <= 5)",
                                                                                              "Medium (5 < X <= 20)",
                                                                                              "Large (20 < X <= 50)",
                                                                                              "Hyperexpanded (50 < X <= 200)")))

cell_abun <- as.data.frame(table(spleen_brain_t_cells@active.ident, spleen_brain_t_cells@meta.data$ebv_pbs))

ggplot(cell_abun, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity", colour = "black") + scale_fill_manual(values = rev(colorblind_vector)) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))


############
# TCR analysis
############

# Correct gene calling

hyper_expansions <- names(table(spleen_brain_t_cells@meta.data[["CTgene"]][spleen_brain_t_cells@meta.data[["cloneType"]] == "Hyperexpanded (50 < X <= 200)"]))

large_expansions <- names(table(spleen_brain_t_cells@meta.data[["CTgene"]][spleen_brain_t_cells@meta.data[["cloneType"]] == "Large (20 < X <= 50)"]))

hyper_genes_in_large <- hyper_expansions[hyper_expansions %in% large_expansions]

table(spleen_brain_t_cells@meta.data[["cloneType"]][spleen_brain_t_cells@meta.data[["CTgene"]] == hyper_genes_in_large[1]])

table(spleen_brain_t_cells@meta.data[["cloneType"]][spleen_brain_t_cells@meta.data[["CTgene"]] == hyper_genes_in_large[2]])


spleen_brain_t_cells@meta.data[["cloneType"]][spleen_brain_t_cells@meta.data[["CTgene"]] == hyper_genes_in_large[1]] <- "Hyperexpanded (50 < X <= 200)"

spleen_brain_t_cells@meta.data[["cloneType"]][spleen_brain_t_cells@meta.data[["CTgene"]] == hyper_genes_in_large[2]] <- "Hyperexpanded (50 < X <= 200)"





# 

Idents(object = spleen_brain_t_cells) <- "cloneType"
Idents(object = spleen_brain_t_cells) <- "seurat_clusters"

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
colorblind_vector <- colorblind_vector[c(1,3,4,5,7)]

DimPlot(spleen_brain_t_cells, label = FALSE, pt.size = 1.5, 
        reduction = "umap.rpca") +
  scale_color_manual(values = rev(colorblind_vector)) + 
  ylab("UMAP 2") + xlab("UMAP 1")

############
# TCR plotting
############

Idents(object = spleen_brain_t_cells) <- "cloneType"
Idents(object = spleen_brain_t_cells) <- "t_cell_states"

spleen_brain_t_cells@meta.data$t_cell_state_organ <- paste0(spleen_brain_t_cells@meta.data$t_cell_states,"_",spleen_brain_t_cells@meta.data$organ)
spleen_brain_t_cells@meta.data$hash_organ <- paste0(spleen_brain_t_cells@meta.data$hash.ID,"_",spleen_brain_t_cells@meta.data$organ)


levels(spleen_brain_t_cells@active.ident)
spleen_brain_t_cells@active.ident <- factor(spleen_brain_t_cells@active.ident, levels = rev(c("Single (0 < X <= 1)",
                                                                                          "Small (1 < X <= 5)",
                                                                                          "Medium (5 < X <= 20)",
                                                                                          "Large (20 < X <= 50)",
                                                                                          "Hyperexpanded (50 < X <= 200)")))

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
colorblind_vector <- colorblind_vector[c(1,3,4,5,7)]

cell_abun <- as.data.frame(table(spleen_brain_t_cells@active.ident, spleen_brain_t_cells@meta.data$hash_organ))

cell_abun$Var2 <- as.factor(cell_abun$Var2)
cell_abun$Var2 <- factor(cell_abun$Var2, levels = c(paste0("Hash",c(1,9,13,10,12,14,17,2,21,23,4,6,8),"_brain"),
                                                    paste0("Hash",c(1,9,13,10,12,14,17,2,21,23,4,6,8),"_spleen")))

cell_abun <- as.data.frame(table(spleen_brain_t_cells@active.ident[spleen_brain_t_cells@meta.data$ebv_pbs == "EBV"], spleen_brain_t_cells@meta.data$t_cell_state_organ[spleen_brain_t_cells@meta.data$ebv_pbs == "EBV"]))


ggplot(cell_abun, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity", color = "black") + scale_fill_manual(values = rev(colorblind_vector)) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black")) + rotate_x_text(45)

############
# TCR brain vs. spleen
############

condition <- "EBV"

percent_expansions_brain <- as.data.frame(table(spleen_brain_t_cells@meta.data$cloneType[spleen_brain_t_cells@meta.data$organ == "brain" & spleen_brain_t_cells@meta.data$ebv_pbs == condition]) / sum(table(spleen_brain_t_cells@meta.data$cloneType[spleen_brain_t_cells@meta.data$organ == "brain" & spleen_brain_t_cells@meta.data$ebv_pbs == condition])) * 100)
percent_expansions_brain$cells <- rep("Brain", times = nrow(percent_expansions_brain))

percent_expansions_spleen <- as.data.frame(table(spleen_brain_t_cells@meta.data$cloneType[spleen_brain_t_cells@meta.data$organ == "spleen" & spleen_brain_t_cells@meta.data$ebv_pbs == condition]) / sum(table(spleen_brain_t_cells@meta.data$cloneType[spleen_brain_t_cells@meta.data$organ == "spleen" & spleen_brain_t_cells@meta.data$ebv_pbs == condition])) * 100)
percent_expansions_spleen$cells <- rep("Spleen", times = nrow(percent_expansions_spleen))

percent_expansions <- rbind(percent_expansions_brain,percent_expansions_spleen)

percent_expansions$Var1 <- as.factor(percent_expansions$Var1)
percent_expansions$Var1 <- factor(percent_expansions$Var1, levels = rev(c("Single (0 < X <= 1)","Small (1 < X <= 5)","Medium (5 < X <= 20)",
                                                                          "Large (20 < X <= 50)","Hyperexpanded (50 < X <= 200)")))



ggbarplot(percent_expansions, x = "cells", y = "Freq", fill = "Var1", palette = rev(colorblind_vector),
          add = c("mean_se"), 
          add.params = list(size = 1, color = "black")) + 
  scale_y_continuous(expand = c(0, 0)) + ylab("% of cells") + xlab("") + ggdist::theme_ggdist() + theme(legend.position = "right",
                                                                                                        axis.text.x = element_text(size = 18, colour = "black"),
                                                                                                        axis.text.y = element_text(size = 18, colour = "black"),
                                                                                                        axis.title.y = element_text(size = 20, colour = "black"),
                                                                                                        axis.ticks.length=unit(.25, "cm"),
                                                                                                        axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                        axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                        axis.ticks = element_line(color="black")) + rotate_x_text(45) + labs(fill='') 




############
# TCR data
############

bcr.no <- 1
table(spleen_brain_t_cells@meta.data$organ[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]]) # / table(spleen_brain_b_cells@meta.data$organ)
(table(spleen_brain_t_cells@meta.data$organ[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]]) / table(spleen_brain_b_cells@meta.data$organ)) * 100
table(spleen_brain_t_cells@meta.data$ebv_pbs[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]])
table(spleen_brain_t_cells@meta.data$hash.ID[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]])
table(spleen_brain_t_cells@meta.data$cloneType[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]])
table(spleen_brain_t_cells@meta.data$CTgene[spleen_brain_t_cells@meta.data$CTgene == top_genes[c(bcr.no)]])



############
# Shared TCRs
############

no.of.genes <- 20

gene_tcr_table <- rev(sort(table(spleen_brain_t_cells@meta.data[["CTgene"]])))
top_genes <- names(gene_tcr_table[c(1:no.of.genes)])
top_genes_counts <- as.numeric(gene_tcr_table[c(1:no.of.genes)])
gene_tcr_table[c(1:no.of.genes)]

df.for.heatmap <- spleen_brain_t_cells@meta.data[spleen_brain_t_cells@meta.data$CTgene %in% top_genes,c("organ","hash.ID","ebv_pbs","CTgene")]

#  df.for.heatmap <- df.for.heatmap[df.for.heatmap$ebv_pbs == "ebv"]

df.for.heatmap$hash_organ <- paste0(df.for.heatmap$hash.ID, "_", df.for.heatmap$organ)

library(dplyr)
library(tidyr)

heatmap_data <- df.for.heatmap %>%
  group_by(hash_organ, CTgene) %>%
  summarise(Count = n(), .groups = "drop")

heatmap_data_wide <- heatmap_data %>%
  pivot_wider(names_from = hash_organ, values_from = Count, values_fill = 0)





heatmap_data_wide <- as.data.frame(heatmap_data_wide)
rownames(heatmap_data_wide) <- heatmap_data_wide[,1]
heatmap_data_wide <- heatmap_data_wide[,c(2:27)]

heatmap_data_wide <- heatmap_data_wide[top_genes, , drop = FALSE]

hold.names <- rownames(heatmap_data_wide)


row_sums <- rowSums(heatmap_data_wide)


# Total cells in brain and spleen


hash_organ_total_counts <- table(spleen_brain_t_cells@meta.data$hash_organ)

hash_organ_total_counts_ebv <- hash_organ_total_counts[names(hash_organ_total_counts) %in% c("Hash9_brain","Hash9_spleen",
                                                                                             "Hash12_brain","Hash12_spleen",
                                                                                             "Hash14_brain","Hash14_spleen",
                                                                                             "Hash17_brain","Hash17_spleen",
                                                                                             "Hash2_brain","Hash2_spleen",
                                                                                             "Hash21_brain","Hash21_spleen",
                                                                                             "Hash23_brain","Hash23_spleen",
                                                                                             "Hash4_brain","Hash4_spleen",
                                                                                             "Hash13_brain","Hash13_spleen",
                                                                                             "Hash8_brain","Hash8_spleen")]


hash_organ_total_counts_pbs <- hash_organ_total_counts[names(hash_organ_total_counts) %in% c("Hash1_brain","Hash1_spleen",
                                                                                             "Hash6_brain","Hash6_spleen",
                                                                                             "Hash10_brain","Hash10_spleen")]

# Separate PBS and EBV

heatmap_data_wide_ebv <- heatmap_data_wide[,c("Hash9_brain","Hash9_spleen",
                                              "Hash12_brain","Hash12_spleen",
                                              "Hash14_brain","Hash14_spleen",
                                              "Hash17_brain","Hash17_spleen",
                                              "Hash2_brain","Hash2_spleen",
                                              "Hash21_brain","Hash21_spleen",
                                              "Hash23_brain","Hash23_spleen",
                                              "Hash4_brain","Hash4_spleen",
                                              "Hash13_brain","Hash13_spleen",
                                              "Hash8_brain","Hash8_spleen")]

heatmap_data_wide_pbs <- heatmap_data_wide[,c("Hash1_brain","Hash1_spleen",
                                              "Hash6_brain","Hash6_spleen",
                                              "Hash10_brain","Hash10_spleen")]


row_sums_ebv <- rowSums(heatmap_data_wide_ebv)
row_sums_pbs <- rowSums(heatmap_data_wide_pbs)

hold.names.col.ebv <- colnames(heatmap_data_wide_ebv)
hold.names.col.pbs <- colnames(heatmap_data_wide_pbs)

# Normalise
#heatmap_data_wide_ebv <- apply(heatmap_data_wide_ebv, MARGIN = 2, scale)
#heatmap_data_wide_pbs <- apply(heatmap_data_wide_pbs, MARGIN = 2, scale)
#heatmap_data_wide_ebv <- t(apply(heatmap_data_wide_ebv, MARGIN = 1, scale))
#heatmap_data_wide_pbs <- t(apply(heatmap_data_wide_pbs, MARGIN = 1, scale))
#heatmap_data_wide <- log2(heatmap_data_wide)
#heatmap_data_wide[heatmap_data_wide == "-Inf"] = 0

#heatmap_data_wide_ebv[heatmap_data_wide_ebv == "NaN"] = 0
#heatmap_data_wide_pbs[heatmap_data_wide_pbs == "NaN"] = 0

for(i in 1:no.of.genes){
  heatmap_data_wide_ebv[i,] <- (heatmap_data_wide_ebv[i,] / hash_organ_total_counts_ebv) * 100
}


for(i in 1:no.of.genes){
  heatmap_data_wide_pbs[i,] <- (heatmap_data_wide_pbs[i,] / hash_organ_total_counts_pbs) * 100
}



colnames(heatmap_data_wide_ebv) <- hold.names.col.ebv # for MARGIN = 1
colnames(heatmap_data_wide_pbs) <- hold.names.col.pbs

library(viridis)
library(circlize)
library(grid)
viridis_colors <- viridis(100)

#custom_groups <- as.numeric(c("1","1","9","9","13",
#                   "13","10","10","12","14",
#                   "14","17","2","2","21",
#                   "21","23","23","4","4",
#                   "6","6","8","8"))

custom_groups <- as.numeric(c("9","9","12","12","14",
                              "14","17","17","2","2","21",
                              "21","23","23","4","4",
                              "13","13","8","8"))
custom_groups <- as.factor(custom_groups)
#custom_groups <- factor(custom_groups, levels = c(1,9,13,10,12,14,17,2,21,23,4,6,8))
custom_groups <- factor(custom_groups, levels = c(9,12,14,17,2,21,23,4,13,8))


custom_groups_pbs <- as.numeric(c("1","1","6","6","10","10"))
custom_groups_pbs <- as.factor(custom_groups_pbs)
custom_groups_pbs <- factor(custom_groups_pbs, levels = c(1,6,10))



#row.ordering <- c("Hash1_brain","Hash1_spleen","Hash9_brain","Hash9_spleen","Hash13_brain","Hash13_spleen","Hash10_brain","Hash10_spleen",
#                  "Hash12_spleen","Hash14_brain","Hash14_spleen","Hash17_spleen","Hash2_brain","Hash2_spleen","Hash21_brain","Hash21_spleen",
#                  "Hash23_brain","Hash23_spleen","Hash4_brain","Hash4_spleen","Hash6_brain","Hash6_spleen","Hash8_brain","Hash8_spleen")

row.ordering <- c("Hash9_brain","Hash9_spleen","Hash12_brain",
                  "Hash12_spleen","Hash14_brain","Hash14_spleen","Hash17_brain","Hash17_spleen","Hash2_brain","Hash2_spleen","Hash21_brain","Hash21_spleen",
                  "Hash23_brain","Hash23_spleen","Hash4_brain","Hash4_spleen","Hash13_brain","Hash13_spleen","Hash8_brain","Hash8_spleen")

row.ordering.pbs <- c("Hash1_brain","Hash1_spleen","Hash6_brain","Hash6_spleen","Hash10_brain","Hash10_spleen")


heatmap_data_wide_ebv <- heatmap_data_wide_ebv[, row.ordering, drop = FALSE]
rownames(heatmap_data_wide_ebv) <- paste0("TCR",c(1:no.of.genes))

heatmap_data_wide_pbs <- heatmap_data_wide_pbs[, row.ordering.pbs, drop = FALSE]
rownames(heatmap_data_wide_pbs) <- paste0("TCR",c(1:no.of.genes))


# PBS = 1, 6, 10

library(ComplexHeatmap)

col_anno <- columnAnnotation(
  Counts = anno_barplot(row_sums_ebv, 
                        gp = gpar(fill = "black"),  # Bar color
                        border = FALSE)  # Remove borders
)

#row_annotation <- data.frame(
#  group = c(rep("PBS",times = 6), rep("EBV",times = 18))
#)

# Convert the 'group' column to a factor
#row_annotation$group <- factor(row_annotation$group, levels = c("PBS", "EBV"))

#group_colors <- c("PBS" = "#797979", "EBV" = "#5194C7")

#row_anno <- rowAnnotation(Group = row_annotation$group,
#                          col = list(Group = group_colors))

heatmap_data_wide_ebv <- log2(heatmap_data_wide_ebv + 1)

library(tibble)

#heatmap_data_wide_ebv <- add_column(heatmap_data_wide_ebv, Hash12_brain = rep(0, nrow(heatmap_data_wide_ebv)), 
#                                    .before = "Hash12_spleen")

#heatmap_data_wide_ebv <- add_column(heatmap_data_wide_ebv, Hash17_brain = rep(0, nrow(heatmap_data_wide_ebv)), 
#                                    .before = "Hash17_spleen")

ComplexHeatmap::Heatmap(t(heatmap_data_wide_ebv), 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        row_split = custom_groups,
                        
                        #left_annotation = row_anno,
                        col = colorRamp2(seq(min(heatmap_data_wide_ebv), 5, length.out = 100), viridis_colors),
                        rect_gp = gpar(col = "black", lwd = 1),
                        top_annotation = col_anno,
                        name = "log2(% TCR of sample cells)")





col_anno_pbs <- columnAnnotation(
  Counts = anno_barplot(row_sums_pbs, 
                        ylim = c(0, 150),
                        gp = gpar(fill = "black"),  # Bar color
                        border = FALSE)  # Remove borders
)

heatmap_data_wide_pbs <- log2(heatmap_data_wide_pbs + 1)

ComplexHeatmap::Heatmap(t(heatmap_data_wide_pbs), 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        row_split = custom_groups_pbs,
                        
                        #left_annotation = row_anno,
                        col = colorRamp2(seq(min(heatmap_data_wide_pbs), 6, length.out = 100), viridis_colors),
                        rect_gp = gpar(col = "black", lwd = 1),
                        top_annotation = col_anno_pbs,
                        name = "log2(% BCR of sample cells)")

############
# Clonotype plotting on UMAP
############

Idents(object = spleen_brain_t_cells) <- "CTgene"

no.tcr <- 1
highlight_cluster <- names(gene_tcr_table[c(1:no.of.genes)])[no.tcr]

cluster_ids <- spleen_brain_t_cells@active.ident
cluster_colors <- ifelse(cluster_ids == highlight_cluster, colorblind_vector[2], "grey")

DimPlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
        reduction = "umap.rpca", cols = cluster_colors, 
        split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") + ggtitle(paste0("TCR",no.tcr," - ",highlight_cluster))


############
# TCR 1 DGE
############

library(EnhancedVolcano)

Idents(object = spleen_brain_t_cells) <- "CTgene"

tcr1.dge <- FindMarkers(spleen_brain_t_cells, ident.1 = names(gene_tcr_table[c(1:no.of.genes)])[1])

Idents(object = spleen_brain_t_cells) <- "cd4_cd8"
cd4_t_cells <- subset(spleen_brain_t_cells, idents = "CD4")
Idents(object = cd4_t_cells) <- "CTgene"

tcr1.dge.cd4 <- FindMarkers(cd4_t_cells, 
                            ident.1 = names(gene_tcr_table[c(1:no.of.genes)])[1], 
                            recorrect_umi = FALSE)

VlnPlot(spleen_brain_t_cells, features = c("S1PR1"), split.by = "organ", idents = c(17,30))

VlnPlot(spleen_brain_t_cells, features = c("S1PR1"))

keyvals <- ifelse(
  tcr1.dge.cd4$avg_log2FC < -2, '#DBE7F6',
  ifelse(tcr1.dge.cd4$avg_log2FC > 2, '#541C25',
         '#C4C4C4'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#541C25'] <- 'high'
names(keyvals)[keyvals == '#C4C4C4'] <- 'mid'
names(keyvals)[keyvals == '#DBE7F6'] <- 'low'

EnhancedVolcano(tcr1.dge.cd4,
                lab = rownames(tcr1.dge.cd4),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #cutoffLineType = 'blank',
                selectLab = c("CD27","SELL","CTSH","IL26",
                              "HRH4","RORC","CCL20","GNLY","HLA-DRB1",
                              "CEBPD"),
                cutoffLineWidth = 0.7,
                pointSize = 1.2,
                pCutoff = 0.05,
                FCcutoff = 2,
                labFace = 'italic',
                colCustom = keyvals,
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + xlab("Log2FC") + ylab("-Log10(Padj)") + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                  legend.position = "none",
                                                                                                                  axis.text.x = element_text(size = 16, colour = "black"),
                                                                                                                  axis.text.y = element_text(size = 16, colour = "black"),
                                                                                                                  axis.title.y = element_text(size = 18, colour = "black"),
                                                                                                                  axis.title.x = element_text(size = 18, colour = "black"),
                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.ticks = element_line(color="black"))


############
# TCR Phenotype
############

Idents(object = spleen_brain_t_cells) <- "seurat_clusters"
spleen_brain_t_cells@meta.data$original.manual.clusters <- spleen_brain_t_cells@active.ident

spleen_brain_t_cells <- FindClusters(spleen_brain_t_cells, resolution = 5, cluster.name = "rpca_clusters")

DimPlot(spleen_brain_t_cells, label = TRUE, pt.size = 1, 
        reduction = "umap.rpca") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1")

spleen_brain_t_cells <- PrepSCTFindMarkers(spleen_brain_t_cells)

tcr1.cluster17.dge <- FindMarkers(spleen_brain_t_cells, ident.1 = 17)

tcr1.cluster17.vs.30.dge <- FindMarkers(spleen_brain_t_cells, ident.1 = 17,
                                        ident.2 = 30)

# rds file checkpoint 2 - just t cells
saveRDS(spleen_brain_t_cells,"D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint7_just_t_cells.rds")
spleen_brain_t_cells <- readRDS("D:/Tbet_B_cell/scRNAseq_Exp2/rds_files/spleen_brain_checkpoint7_just_t_cells.rds")



spleen_brain_t_cells@meta.data$ebv_pbs <- spleen_brain_t_cells@meta.data[["hash.ID"]]

spleen_brain_t_cells@meta.data$ebv_pbs[spleen_brain_t_cells@meta.data$ebv_pbs %in% c("Hash9","Hash12",
                                                                                     "Hash14","Hash17",
                                                                                     "Hash2","Hash21",
                                                                                     "Hash23","Hash4",
                                                                                     "Hash13","Hash8")] <- "EBV"

spleen_brain_t_cells@meta.data$ebv_pbs[spleen_brain_t_cells@meta.data$ebv_pbs %in% c("Hash1",
                                                                                     "Hash6","Hash10")] <- "PBS"

table(spleen_brain_t_cells@meta.data$ebv_pbs)


spleen_brain_t_cells@meta.data$mouse_id <- spleen_brain_t_cells@meta.data[["hash.ID"]]

spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash2")] <- "M22"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash4")] <- "M26"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash6")] <- "M33"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash8")] <- "M7"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash10")] <- "M8"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash12")] <- "M13"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash14")] <- "M16"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash17")] <- "M5"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash21")] <- "M6"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash23")] <- "M21"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash1")] <- "M11"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash9")] <- "M35"
spleen_brain_t_cells@meta.data$mouse_id[spleen_brain_t_cells@meta.data$mouse_id %in% c("Hash13")] <- "M3"

table(spleen_brain_t_cells@meta.data$mouse_id)



keyvals <- ifelse(
  tcr1.cluster17.vs.30.dge$avg_log2FC < -1.5, '#DBE7F6',
  ifelse(tcr1.cluster17.vs.30.dge$avg_log2FC > 1.5, '#541C25',
         '#C4C4C4'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#541C25'] <- 'high'
names(keyvals)[keyvals == '#C4C4C4'] <- 'mid'
names(keyvals)[keyvals == '#DBE7F6'] <- 'low'

EnhancedVolcano(tcr1.cluster17.vs.30.dge,
                lab = rownames(tcr1.cluster17.vs.30.dge),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #cutoffLineType = 'blank',
                #selectLab = c("CD27","SELL","CTSH","IL26",
                #              "HRH4","RORC","CCL20","GNLY","HLA-DRB1",
                #              "CEBPD"),
                cutoffLineWidth = 0.7,
                pointSize = 1.2,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                labFace = 'italic',
                colCustom = keyvals,
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + xlab("Log2FC") + ylab("-Log10(Padj)") + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                  legend.position = "none",
                                                                                                                  axis.text.x = element_text(size = 16, colour = "black"),
                                                                                                                  axis.text.y = element_text(size = 16, colour = "black"),
                                                                                                                  axis.title.y = element_text(size = 18, colour = "black"),
                                                                                                                  axis.title.x = element_text(size = 18, colour = "black"),
                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.ticks = element_line(color="black"))



############
# Th17 signature
############

th17.signature <- c("RORC","RORA","CCL20","CCR6","CTSH","IL26","IL23R",
                    "IFNG","TNF","IL17A","CXCR3")

spleen_brain_t_cells <- Seurat::AddModuleScore(spleen_brain_t_cells,
                                               features = list(th17.signature), name = "TH17_")

FeaturePlot(spleen_brain_t_cells, label = TRUE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "TH17_1", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")
#+ ggtitle("TH17 Gene Signature")

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 2, order = T,
reduction = "umap.rpca", max.cutoff = 'q97', 
features = "STAT4", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")

library(tidyr)

Idents(object = spleen_brain_t_cells) <- "seurat_clusters"

cell_abun <- as.data.frame(table(spleen_brain_t_cells@active.ident, spleen_brain_t_cells@meta.data$organ))

#cell_abun <- cell_abun[cell_abun$Var1 %in% c(17,30),]
#cell_abun$Var2 <- as.factor(cell_abun$Var2)
#cell_abun$Var2 <- factor(cell_abun$Var2, levels = c(paste0("Hash",c(1,9,13,10,12,14,17,2,21,23,4,6,8))))


ggplot(cell_abun, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black")) + rotate_x_text(45)

############
# TCR1 in spleen vs. brain
############

Idents(object = spleen_brain_t_cells) <- "CTgene"
tcr1_seurat <- subset(spleen_brain_t_cells, idents = names(gene_tcr_table[c(1:no.of.genes)])[1])
Idents(object = tcr1_seurat) <- "organ"

tcr1.brain.vs.spleen.dge <- FindMarkers(tcr1_seurat, ident.1 = "brain",
                                        ident.2 = "spleen", 
                                        recorrect_umi = FALSE)

keyvals <- ifelse(
  tcr1.brain.vs.spleen.dge$avg_log2FC < -1, '#DBE7F6',
  ifelse(tcr1.brain.vs.spleen.dge$avg_log2FC > 1, '#541C25',
         '#C4C4C4'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#541C25'] <- 'high'
names(keyvals)[keyvals == '#C4C4C4'] <- 'mid'
names(keyvals)[keyvals == '#DBE7F6'] <- 'low'

EnhancedVolcano(tcr1.brain.vs.spleen.dge,
                lab = rownames(tcr1.brain.vs.spleen.dge),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #cutoffLineType = 'blank',
                #selectLab = c("CD27","SELL","CTSH","IL26",
                #              "HRH4","RORC","CCL20","GNLY","HLA-DRB1",
                #              "CEBPD"),
                cutoffLineWidth = 0.7,
                pointSize = 1.2,
                pCutoff = 0.05,
                FCcutoff = 1,
                labFace = 'italic',
                colCustom = keyvals,
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + xlab("Log2FC") + ylab("-Log10(Padj)") + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                  legend.position = "none",
                                                                                                                  axis.text.x = element_text(size = 16, colour = "black"),
                                                                                                                  axis.text.y = element_text(size = 16, colour = "black"),
                                                                                                                  axis.title.y = element_text(size = 18, colour = "black"),
                                                                                                                  axis.title.x = element_text(size = 18, colour = "black"),
                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.ticks = element_line(color="black"))

############
# Th17 spleen vs. brain
############

Idents(object = spleen_brain_t_cells) <- "seurat_clusters"
th17_seurat <- subset(spleen_brain_t_cells, idents = c(17,30))
Idents(object = th17_seurat) <- "organ"

th17.brain.vs.spleen.dge <- FindMarkers(th17_seurat, ident.1 = "brain",
                                        ident.2 = "spleen", 
                                        recorrect_umi = FALSE)

keyvals <- ifelse(
  th17.brain.vs.spleen.dge$avg_log2FC < -1, '#DBE7F6',
  ifelse(th17.brain.vs.spleen.dge$avg_log2FC > 1, '#541C25',
         '#C4C4C4'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#541C25'] <- 'high'
names(keyvals)[keyvals == '#C4C4C4'] <- 'mid'
names(keyvals)[keyvals == '#DBE7F6'] <- 'low'

EnhancedVolcano(th17.brain.vs.spleen.dge,
                lab = rownames(th17.brain.vs.spleen.dge),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #cutoffLineType = 'blank',
                #selectLab = c("CD27","SELL","CTSH","IL26",
                #              "HRH4","RORC","CCL20","GNLY","HLA-DRB1",
                #              "CEBPD"),
                cutoffLineWidth = 0.7,
                pointSize = 1.2,
                pCutoff = 0.05,
                FCcutoff = 1,
                labFace = 'italic',
                colCustom = keyvals,
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 4.0,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + xlab("Log2FC") + ylab("-Log10(Padj)") + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                  legend.position = "none",
                                                                                                                  axis.text.x = element_text(size = 16, colour = "black"),
                                                                                                                  axis.text.y = element_text(size = 16, colour = "black"),
                                                                                                                  axis.title.y = element_text(size = 18, colour = "black"),
                                                                                                                  axis.title.x = element_text(size = 18, colour = "black"),
                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.ticks = element_line(color="black"))


############
# Chemotaxis Analysis
############

chemotaxis.genes <- c("CCR7","CXCR4","CCR5","CXCR3","CCR6","CCR4",
                      "CCR8","CCR9","CCR10","CXCR5","CX3CR1","CCR2",
                      "CXCR6","CSF2")

Idents(object = spleen_brain_t_cells) <- "t_cell_states"
tem_help_seurat <- subset(spleen_brain_t_cells, idents = c("Tem Helper"))
Idents(object = tem_help_seurat) <- "organ"

them.helper.spleen.vs.brain.dge <- FindMarkers(tem_help_seurat, ident.1 = "brain",
                                        ident.2 = "spleen",
                                        recorrect_umi = FALSE)

them.helper.spleen.vs.brain.dge <- them.helper.spleen.vs.brain.dge[chemotaxis.genes,]

keyvals <- ifelse(
  them.helper.spleen.vs.brain.dge$avg_log2FC < -1, '#DBE7F6',
  ifelse(them.helper.spleen.vs.brain.dge$avg_log2FC > 1, '#541C25',
         '#C4C4C4'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#541C25'] <- 'high'
names(keyvals)[keyvals == '#C4C4C4'] <- 'mid'
names(keyvals)[keyvals == '#DBE7F6'] <- 'low'

EnhancedVolcano(them.helper.spleen.vs.brain.dge,
                lab = rownames(them.helper.spleen.vs.brain.dge),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #cutoffLineType = 'blank',
                selectLab = chemotaxis.genes,
                cutoffLineWidth = 0.7,
                pointSize = 3,
                pCutoff = 0.05,
                FCcutoff = 1,
                labFace = 'italic',
                colCustom = keyvals,
                #boxedLabels = TRUE,
                #drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "",
                labSize = 6.0,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + xlab("Log2FC") + ylab("-Log10(Padj)") + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                  legend.position = "none",
                                                                                                                  axis.text.x = element_text(size = 16, colour = "black"),
                                                                                                                  axis.text.y = element_text(size = 16, colour = "black"),
                                                                                                                  axis.title.y = element_text(size = 18, colour = "black"),
                                                                                                                  axis.title.x = element_text(size = 18, colour = "black"),
                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                  axis.ticks = element_line(color="black"))



library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)

chemotaxis.genes <- c("CCR7","CXCR4","CCR5","CXCR3","CCR6","CCR4",
                      "CCR8","CCR9","CCR10","CXCR5","CX3CR1","CCR2",
                      "CXCR6","CSF2")

Idents(object = spleen_brain_t_cells) <- "ebv_pbs"
ebv_seurat <- subset(spleen_brain_t_cells, idents = c("PBS"))
Idents(object = ebv_seurat) <- "t_cell_state_organ"

avg_exp_chemotaxis <- as.data.frame(AverageExpression(ebv_seurat, features = chemotaxis.genes, assays = c("SCT"))$SCT)

viridis_colors <- viridis(100)

ComplexHeatmap::Heatmap(t(avg_exp_chemotaxis), 
                        cluster_rows = TRUE, 
                        cluster_columns = TRUE,
                        #left_annotation = row_anno,
                        col = colorRamp2(seq(0, 4, length.out = 100), viridis_colors),
                        rect_gp = gpar(col = "black", lwd = 1),
                        #top_annotation = col_anno_pbs,
                        name = "Average Normalised Expression")

##################
# Extract top TCRs to a excel file
##################

top_genes <- names(rev(sort(table(spleen_brain_t_cells@meta.data[["CTgene"]]))))
top_genes[c(1:100)]

top_clones <- rev(sort(table(spleen_brain_t_cells@meta.data[["CTgene"]])))
amount_clones <- rev(as.numeric(table(spleen_brain_t_cells@meta.data[["CTgene"]])))

datalist = list()
count = 0

for(i in 1:100){
  
  num_genes <- length(unique(na.omit(spleen_brain_t_cells@meta.data[["CTgene"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]])))
  
  for(j in 1:num_genes){
    
    num_rearrange <- length(unique(na.omit(spleen_brain_t_cells@meta.data[["CTgene"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]])))
    
    for(k in 1:num_rearrange){
      
      count = count + 1
      
      amount_organ <- table(spleen_brain_t_cells@meta.data[["organ"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]])
      
      datalist[[count]] <- c(i, as.numeric(top_clones[i]),
                             unique(na.omit(spleen_brain_t_cells@meta.data[["CTgene"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]]))[k],
                             unique(na.omit(spleen_brain_t_cells@meta.data[["CTaa"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]]))[j],
                             if (length(amount_organ[names(amount_organ) == "brain"]) != 0) {
                               amount_organ[names(amount_organ) == "brain"]
                             } else {
                               0
                             },
                             if (length(amount_organ[names(amount_organ) == "spleen"]) != 0) {
                               amount_organ[names(amount_organ) == "spleen"]
                             } else {
                               0
                             },
                             unique(na.omit(spleen_brain_t_cells@meta.data[["ebv_pbs"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]]))[1],
                             names(which.max(table(as.character(spleen_brain_t_cells@meta.data[["cd4_cd8"]][spleen_brain_t_cells@meta.data[["CTgene"]] == top_genes[i]]))))
      )
      
    }
    
  }
  
}

big_Data <- as.data.frame(do.call(rbind, datalist))
colnames(big_Data) <- c("Clone Number","Size","Gene","Amino Acids","# in brain","# in spleen","EBV/PBS","CD4/CD8")

openxlsx::write.xlsx(big_Data, "D:/Tbet_B_cell/scRNAseq_Exp2/tcr_sequences.xlsx")


openxlsx::write.xlsx(top_genes[c(1:50)], "D:/Tbet_B_cell/scRNAseq_Exp2/bcr_sequences.xlsx")






FeaturePlot(spleen_brain_t_cells, features = "IFNG", 
            order = T, max.cutoff = 'q97', 
            reduction = "umap.rpca", pt.size = 1, split.by = "ebv_pbs") + xlab("UMAP 1") + ylab("UMAP 2") & 
  scale_color_viridis(option = "D") 



ebv_pbs <- FindMarkers(spleen_brain_t_cells, ident.1 = "EBV", ident.2 = "PBS", features = "CD40LG")



Idents(spleen_brain_t_cells) <- "ebv_pbs"
Idents(spleen_brain_t_cells) <- "t_cell_states"

spleen_brain_t_cells_teff_help <- subset(spleen_brain_t_cells, idents = "Teff Helper")

Idents(spleen_brain_t_cells_teff_help) <- "ebv_pbs"

Idents(spleen_brain_t_cells_teff_help) <- factor(Idents(spleen_brain_t_cells_teff_help), levels = c("PBS","EBV"))

VlnPlot(spleen_brain_t_cells_teff_help, features = "CD40LG", 
        pt.size = 0) + xlab("") + theme(legend.position = "none")

spleen_brain_t_cells_teff_help <- PrepSCTFindMarkers(spleen_brain_t_cells_teff_help)

ebv_pbs <- FindMarkers(spleen_brain_t_cells_teff_help, ident.1 = "EBV", ident.2 = "PBS", features = "CD40LG")


##################
# HLA-DR, GrA, GrB
##################


FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "HLA-DRA", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "HLA-DRA", split.by = "ebv_pbs") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")



FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "GZMB", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "GZMB", split.by = "ebv_pbs") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")





FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "GZMA", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "IFNG", split.by = "ebv_pbs") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")





tp_cell_ids <- WhichCells(spleen_brain_t_cells, 
                           expression = `HLA-DRA` > 0 & GZMB > 0 & GZMA > 0)

pop.name <- "HLA-DRA+GMBA+GZMB+"

spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB <- rep("NOI", times = nrow(spleen_brain_t_cells@meta.data))
spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB[rownames(spleen_brain_t_cells@meta.data) %in% tp_cell_ids] <- pop.name

spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB <- as.factor(spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB)
spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB <- factor(spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB, levels = c("NOI",pop.name))

Idents(object = spleen_brain_t_cells) <- "TP_HLADR_GMZAB"
Idents(object = spleen_brain_t_cells) <- "cd4_cd8"

DimPlot(spleen_brain_t_cells, reduction = "umap.rpca", pt.size = 1, split.by = "organ") + scale_colour_manual(values = c("grey","red")) + 
  ylab("UMAP 2") + xlab("UMAP 1")

DimPlot(spleen_brain_t_cells, reduction = "tsne.rpca", pt.size = 1, order = TRUE) + 
  scale_colour_manual(values = c("grey","red")) + 
  ylab("tSNE 2") + xlab("tSNE 1")

table(spleen_brain_t_cells@meta.data$TP_HLADR_GMZAB[spleen_brain_t_cells@meta.data$ebv_pbs == "EBV" & spleen_brain_t_cells@meta.data$organ == "brain" & spleen_brain_t_cells@meta.data$cd4_cd8 == "CD8"])

autoAggression <- readxl::read_xlsx("D:/Tbet_B_cell/scRNAseq_Exp2/autoAggression.xlsx")

autoAggression <- autoAggression$Genes

spleen_brain_t_cells <- UCell::AddModuleScore_UCell(spleen_brain_t_cells, features = list(AutoAggression = autoAggression))

spleen_brain_t_cells@meta.data$AutoAggression_UCell

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 1, 
            reduction = "umap.rpca",  
            features = "AutoAggression_UCell", split.by = "organ") + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")

##################
# Macrophage activating T cell cytokines - only brain of EBV-infected mice (Fig. 5c)
##################

Idents(object = spleen_brain_t_cells) <- "organ"
Idents(object = spleen_brain_t_cells) <- "seurat_clusters"
Idents(object = spleen_brain_t_cells) <- "ebv_pbs"
Idents(object = spleen_brain_t_cells) <- "predicted.id"
Idents(object = spleen_brain_t_cells) <- "t_cell_states"
Idents(object = spleen_brain_t_cells) <- "cd4_cd8"

ebv_cells <- WhichCells(object = spleen_brain_t_cells,
                        idents = "EBV")

brain_cells <- WhichCells(object = spleen_brain_t_cells,
                        idents = "brain")

ebv_brain_cells <- intersect(ebv_cells, brain_cells)
ebv_brain_cells

FeaturePlot(spleen_brain_t_cells, label = FALSE, pt.size = 2, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "IFNG", split.by = "organ", 
            order = T, cells = ebv_brain_cells) + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")
#+ ggtitle("TH17 Gene Signature")



VlnPlot(spleen_brain_t_cells, features = "IFNG", 
        split.by = "organ")


##################
# 
##################

Idents(object = spleen_brain_t_cells) <- "t_cell_states"
Idents(object = spleen_brain_t_cells) <- "ebv_pbs"

ebv_spleen_brain_t_cells <- subset(spleen_brain_t_cells,  idents = 'EBV')

Idents(object = ebv_spleen_brain_t_cells) <- "t_cell_states"

cd8_em_ebv_spleen_brain_t_cells <- subset(ebv_spleen_brain_t_cells,  idents = 'Tem Cytotoxic')

Idents(object = cd8_em_ebv_spleen_brain_t_cells) <- "organ"

cd8_em_ebv_spleen_brain_t_cells <- PrepSCTFindMarkers(cd8_em_ebv_spleen_brain_t_cells, assay = "SCT", verbose = TRUE)
CD8_EM <- FindMarkers(cd8_em_ebv_spleen_brain_t_cells, ident.1 = "brain", ident.2 = "spleen", recorrect_umi=FALSE, assay = "SCT")

writexl::write_xlsx(cbind(CD8_EM, rownames(CD8_EM)), path = "D:/MS_T_cell/Figures/Fig5/cd8_em_brain_vs_spleen.xlsx")

#inflam_vs_spp1 <- readxl::read_xlsx("D:/MS_T_cell/Figures/Fig3/inflam_vs_spp1.xlsx")

keyvals <- ifelse(
  CD8_EM$avg_log2FC < -1 & CD8_EM$p_val_adj < 10e-10, 'grey35',
  ifelse(CD8_EM$avg_log2FC > 1 & CD8_EM$p_val_adj < 10e-10, '#E5337F',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey35'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == '#E5337F'] <- 'low'

EnhancedVolcano::EnhancedVolcano(CD8_EM,
                                 lab = rownames(CD8_EM),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 FCcutoff = 1,
                                 pCutoff = 10e-10,
                                 selectLab = c("CXCR4", "CD69", "RPS29",
                                               "HLA-A","RGS1","PRF1",
                                               "CRIP1","DUSP2","RGS2","LAG3",
                                               "PDCD1","IRF1","IRF7","CXCR3",
                                               "CREM","ITGAL","IFNG","CCL4",
                                               "CCR5"),
                                 colCustom = keyvals,
                                 drawConnectors = TRUE,
                                 colConnectors = 'black',
                                 typeConnectors = 'closed',
                                 lengthConnectors = unit(0.02,'npc'),
                                 labSize = 5,
                                 colAlpha = 1,
                                 title = "",
                                 titleLabSize = 0,
                                 subtitle = "",
                                 subtitleLabSize = 0,
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE,
                                 legendPosition = "right",
                                 max.overlaps = Inf) +
  theme(axis.ticks.length=unit(.3, "cm"),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none")

##################
# 
##################

Idents(object = spleen_brain_t_cells) <- "t_cell_states"
Idents(object = spleen_brain_t_cells) <- "ebv_pbs"

tregs <- subset(spleen_brain_t_cells,  idents = 'Tregs')

DimPlot(tregs, reduction = "umap.rpca", label = TRUE)

FeaturePlot(tregs, label = FALSE, pt.size = 3, 
            reduction = "umap.rpca", max.cutoff = 'q97', 
            features = "IL10", split.by = "organ_condition", 
            order = T) + theme(legend.position = "none") +
  ylab("UMAP 2") + xlab("UMAP 1") & 
  scale_color_viridis(option = "D")
#+ ggtitle("TH17 Gene Signature")

tregs$organ_condition <- paste0(tregs$organ,"_",tregs$ebv_pbs)
Idents(object = tregs) <- "organ_condition"

organ_condition_dge <- FindAllMarkers(tregs, recorrect_umi=FALSE)

brain_ebv_dge <- organ_condition_dge[organ_condition_dge$cluster == "spleen_EBV",]

brain_ebv_dge_pos <- brain_ebv_dge[brain_ebv_dge$avg_log2FC > 1,]


organ_condition_specific <- FindMarkers(tregs, ident.1 = "spleen_EBV",
                                        ident.2 = "spleen_PBS", recorrect_umi=FALSE)


##################
# Dotplot
##################


pdf("D:/MS_T_cell/Figures/Ex9 - T cell scRNAseq extras/dotplot_ex_9a.pdf", width = 7, height = 4)

DotPlot(spleen_brain_t_cells, features = c("CD4","CD8A","FOXP3",
                                                 "IL2RA","SELL","CCR7",
                                                 "ICOS","GZMB","PRF1",
                                                 "IFNG","HLA-DRA","CD69",
                                                 "ITGAE","MKI67","CCNA2",
                                                 "CCNB1"), cols = "Blues") + 
  rotate_x_text(90) + scale_color_distiller(direction = 0) + geom_point(
    aes(size = pct.exp),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) + xlab("") + ylab("")

dev.off()




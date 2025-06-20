knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, knitr.table.format = "html") 
# input for this report: sce 
library(Seurat)
library(HDF5Array)
library(ezRun)
library(RColorBrewer)
library(SingleCellExperiment)
library(pheatmap)
library(DropletUtils)
library(scuttle)
library(BiocParallel)
library(data.table)
library(ggplot2)
library(ezRun)
library(Seurat)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(Gviz)
library(GenomicFeatures)
library(grDevices)
library(S4Vectors)
library(Seurat)
library(STACAS)
library(SignatuR)
library(patchwork)
library(paletteer)
library(Azimuth)
library(SeuratDisk)
library(scRepertoire)
library(SeuratData)
library(zellkonverter)

spleen.merged <- readRDS(file = paste0(project_wd, "/spleen.merged.rds"))
mycols <- as.character(pals::trubetskoy(12))
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

###### Figure 3b

DimPlot(spleen.merged, 
        reduction = "umap.rpca", 
        label = TRUE, 
        repel = TRUE, 
        label.size = 4) + 
  scale_color_manual(values = mycols) +
  ggtitle("B cell subtypes")

###### Figure 3c & d

FeaturePlot(spleen.merged, 
            features = c("TBX21", "CXCR3"), 
            reduction = "umap.rpca", 
            order = TRUE, 
            max.cutoff = 'q95') & 
  scale_color_viridis(option = "D")

###### Figure 3e

DimPlot(spleen.merged, 
        group.by = "cloneSize", 
        reduction = "umap.rpca", 
        order = TRUE) +
  scale_color_manual(values = rev(colorblind_vector)) +
  ggtitle("BCR expression (cloneSize)")

###### Figure 3f

DefaultAssay(spleen.merged)

percent_expansions_TbetCXCR3 <- as.data.frame(table(spleen.merged@meta.data$cloneSize[spleen.merged@meta.data$DN2_cells == "Tbet+CXCR3+" & spleen.merged@meta.data$predicted.id %in% c("Memory B cells","Cycling immune mix","Plasma cells","Plasmablasts")]) / sum(table(spleen.merged@meta.data$cloneSize[spleen.merged@meta.data$DN2_cells == "Tbet+CXCR3+" & spleen.merged@meta.data$predicted.id %in% c("Memory B cells","Cycling immune mix","Plasma cells","Plasmablasts")])) * 100)
percent_expansions_TbetCXCR3$cells <- rep("Tbet+CXCR3+", times = nrow(percent_expansions_TbetCXCR3))

percent_expansions_NOI <- as.data.frame(table(spleen.merged@meta.data$cloneSize[spleen.merged@meta.data$DN2_cells == "NOI" & spleen.merged@meta.data$predicted.id %in% c("Memory B cells","Cycling immune mix","Plasma cells","Plasmablasts")]) / sum(table(spleen.merged@meta.data$cloneSize[spleen.merged@meta.data$DN2_cells == "NOI" & spleen.merged@meta.data$predicted.id %in% c("Memory B cells","Cycling immune mix","Plasma cells","Plasmablasts")])) * 100)
percent_expansions_NOI$cells <- rep("NOI", times = nrow(percent_expansions_NOI))

percent_expansions <- rbind(percent_expansions_TbetCXCR3,percent_expansions_NOI)

percent_expansions$Var1 <- as.factor(percent_expansions$Var1)
percent_expansions$Var1 <- factor(percent_expansions$Var1, levels = rev(c("Single (0 < X <= 1)",
                                                                          "Small (1 < X <= 5)",
                                                                          "Medium (5 < X <= 20)",
                                                                          "Large (20 < X <= 50)",
                                                                          "Hyperexpanded (50 < X <= 200)",
                                                                          "Hyperexpanded (50 < X <= 649)",
                                                                          "Hyperexpanded (50 < X <= 2078)")))



ggbarplot(percent_expansions, x = "cells", y = "Freq", fill = "Var1", palette = rev(viridis::magma(7)),
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


###### Figure 3g

FeaturePlot(spleen.merged, 
            features = 'EBNA', 
            reduction = "umap.rpca", 
            order = TRUE, 
            max.cutoff = 'q95', 
            ncol = 1) & 
  scale_color_viridis(option = "D")

###### Figure 3h

sce <- readRDS("/sce_trajectory_computed.rds")

# Plotting the trajectories
plot(reducedDims(sce)$UMAP, col = mycols[Idents(spleen.merged)], pch=16)
lines(SlingshotDataSet(sce), col = 'black', lwd = 2)

###### Figure 3i

# Define color map to match DimPlot
cluster_color_map <- c(
  "Memory B cells" = "#3cb44b",
  "Naive B cells" = "#ffe119",
  "Cycling B cells" = "#e6194b",
  "Plasma cells" = "#4363d8",
  "Plasmablasts" = "#f58231"
)

# Get cluster identities
cluster_ids <- Idents(spleen.merged)

# Genes of interest
genes_of_interest <- c("BZLF1", "EBNA", "LMP1", "LMP2A", "TBX21", "CXCR3")

# Predict smooth gene expression
yhatSmooth <- predictSmooth(sce, gene = genes_of_interest, nPoints = 200, tidy = FALSE)

# Get pseudotime values
pseudotime <- slingPseudotime(sce)

# If there are multiple lineages, use the first one
if (is.matrix(pseudotime)) {
  pseudotime <- pseudotime[,1]
}

# Order cells by pseudotime
ordered_cells <- names(sort(pseudotime))

# Get cluster IDs in the order of cells along the pseudotime
orderedClusters <- cluster_ids[ordered_cells]

# Create annotation for top of heatmap
annotation_col <- data.frame(Cluster = orderedClusters[seq(1, length(orderedClusters), length.out = ncol(yhatSmooth))])
rownames(annotation_col) <- colnames(yhatSmooth)

# Create the heatmap
heatSmooth <- pheatmap(t(scale(t(yhatSmooth))), 
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = FALSE,
                       main = "Gene Expression along Slingshot Trajectory",
                       fontsize = 10,
                       scale = "row",
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                       annotation_col = annotation_col,
                       annotation_colors = list(Cluster = cluster_color_map),
                       annotation_legend = TRUE)

# Print the heatmap
print(heatSmooth)

###### Figure 3j

cns_markers <- c("ALCAM", "ITGAL")

plot_data <- FetchData(spleen.merged, vars = c(cns_markers, "dp_status"))
plot_data_long <- pivot_longer(plot_data, cols = all_of(cns_markers), names_to = "gene", values_to = "expression")

ggplot(plot_data_long, aes(x = gene, y = expression, fill = dp_status)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_fill_manual(values = c("Non_DP_cells" = "#E69F00", "DP_cells" = "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of CNS Migration Markers in DP vs Non-DP Cells",
       x = "Gene",
       y = "Expression Level",
       fill = "Cell Type") +
  stat_compare_means(method = "t.test", 
                     label = "p.format",
                     label.x = 1.5) +
  theme_classic()

###### Figure 3k

VlnPlot(spleen.merged, 
        features = "chemotaxis_genes_UCell", 
        group.by = "dp_status", 
        pt.size = 0,
        cols = c("Non_DP_cells" = "#E69F00", "DP_cells" = "#56B4E9")) +
  stat_compare_means(method = "t.test", 
                     label = "p.format",
                     label.x = 1.5) +
  labs(title = "Chemotaxis Genes UCell Scores",
       x = "Cell Type",
       y = "UCell Score") +
  theme_classic()

###### Figure 3l

# List of chemotaxis genes
chemotaxis_genes <- c("CCL3", "CCL4", "CCL5", "CCL19", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCR3", "CCR7")

# Perform differential expression analysis
Idents(spleen.merged) <- "dp_status"
de_results <- FindMarkers(spleen.merged, 
                          ident.1 = "DP_cells", 
                          ident.2 = "Non_DP_cells", 
                          features = intersect(chemotaxis_genes, Features(spleen.merged)))

# Extract log fold changes and p-values
heatmap_data <- de_results[, c("avg_log2FC", "p_val_adj")]
heatmap_data$gene <- rownames(heatmap_data)
heatmap_data <- heatmap_data[order(heatmap_data$avg_log2FC, decreasing = TRUE), ]

# Create color scale
color_scale <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Create heatmap
pheatmap(heatmap_data[, "avg_log2FC", drop = FALSE],
         color = color_scale,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         labels_row = heatmap_data$gene,
         main = "Differential Expression of Chemotaxis Genes (DP_cells vs Non_DP_cells)",
         angle_col = 0,
         display_numbers = matrix(ifelse(heatmap_data$p_val_adj < 0.05, "*", ""), ncol = 1),
         number_color = "black",
         fontsize_number = 17)

###### Extended Data Fig. 4a

DoublePos.cells <- subset(x = spleen.merged, subset = CXCR3 > 0 & TBX21 > 0)

DimPlot(spleen.merged, sizes.highlight = 0.2, cells.highlight = Cells(DoublePos.cells), reduction = "umap.rpca") + ggtitle("CXCR3+ / TBX21+ cells") + 
  scale_color_manual(labels = c("Other", "Double Positive cells"), values = c("grey", "red")) 

###### Extended Data Fig. 4b

dittoBarPlot(spleen.merged, var =  "predicted.id", 
             group.by = "dp_status", retain.factor.levels = T) + 
  scale_fill_manual(values = mycols) 

###### Extended Data Fig. 4b

cell_abun <- as.data.frame(table(spleen.merged@meta.data$predicted.id, spleen.merged@meta.data$DN2_cells))

cell_abun$Var1 <- as.factor(cell_abun$Var1)
cell_abun$Var1 <- factor(cell_abun$Var1, levels = c("Naive B cells",
                                                    "Memory B cells",
                                                    "Cycling immune mix",
                                                    "Plasmablasts",
                                                    "Plasma cells"))

ggplot(cell_abun, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity", colour = "black", linewidth = 1) + scale_fill_manual(values = c("#EBE747","#1BC4D1","#43DB95","#E276BF","#E65658")) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))

###### Extended Data Fig. 4c

pre_char <- sapply(strsplit(spleen.merged@meta.data[["CTgene"]], "_"), `[`, 1)
class_switch <- sapply(strsplit(pre_char, "[.]"), `[`, 4)

spleen.merged@meta.data$class <- class_switch 

spleen.merged@meta.data$isotype <- spleen.merged@meta.data$class
spleen.merged@meta.data$isotype[spleen.merged@meta.data$isotype %in% c("IGHA1","IGHA2")] <- "IgA"
spleen.merged@meta.data$isotype[spleen.merged@meta.data$isotype %in% c("IGHG1","IGHG2",
                                                                     "IGHG3","IGHG4")] <- "IgG"
spleen.merged@meta.data$isotype[spleen.merged@meta.data$isotype %in% c("IGHM","IGHM;IGHV2-5")] <- "IgM"

isotope_cleaning <- spleen.merged@meta.data[!(spleen.merged@meta.data$isotype == "NA"),]
isotope_cleaning <- isotope_cleaning[!is.na(isotope_cleaning$isotype),]


# Calculate the proportion of cells in each isotype category for each cluster
isotype_prop <- isotope_cleaning %>%
  group_by(DN2_cells, isotype) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(DN2_cells) %>%
  mutate(prop = n / sum(n))

ggplot(isotype_prop, aes(fill=isotype, y=prop, x=DN2_cells)) + 
  geom_bar(position="fill", stat="identity", colour = "black", linewidth = 1) + scale_fill_manual(values = rev(c("grey90","grey50","grey20","black"))) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))

###### Extended Data Fig. 4d

genes <- c('IGHG3', 'IGHG1')

FeaturePlot(spleen.merged, 
            features = genes, 
            reduction = 'umap.rpca', 
            order = TRUE, 
            max.cutoff = 'q97', 
            ncol = 2) & 
  scale_color_viridis_c(option = 'D')

###### Extended Data Fig. 4e

genes <- c('BZLF1', 'LMP1', 'LMP2A')

FeaturePlot(spleen.merged, 
            features = genes, 
            reduction = 'umap.rpca', 
            order = TRUE, 
            max.cutoff = 'q97',
            ncol = 3) & 
  scale_color_viridis_c(option = 'D')

###### Extended Data Fig. 4f

clonalOverlap(spleen.merged, group.by = 'predicted.id',
              cloneCall = "strict", 
              method = "morisita")+ RotatedAxis()

###### Extended Data Fig. 4g-h

Idents(spleen.merged) <- "dp_status"

spleen.merged@meta.data$dp_status

spleen.merged <- subset(x = spleen.merged, idents = "DP_cells")

ebv_cells <- WhichCells(spleen.merged, 
                        expression = EBNA > 0)

pop.name <- "EBV+"

spleen.merged@meta.data$ebv_pos <- rep("EBV-", times = nrow(spleen.merged@meta.data))
spleen.merged@meta.data$ebv_pos[rownames(spleen.merged@meta.data) %in% ebv_cells] <- pop.name

spleen.merged@meta.data$ebv_pos <- as.factor(spleen.merged@meta.data$ebv_pos)
spleen.merged@meta.data$ebv_pos <- factor(spleen.merged@meta.data$ebv_pos, levels = c("EBV-",pop.name))

Idents(object = spleen.merged) <- "ebv_pos"

spleen.merged <- AddModuleScore(spleen.merged, features = list(c("MKI67","PCNA","TOP2A",
                                                                             "CCNB1","CCNB2",
                                                                             "CDK1")), name = "PROLIF")

VlnPlot(spleen.merged, features = "MKI67", pt.size = 0) +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95) + 
  theme(legend.position = "none") + ylab("MKI67 Expression") + xlab("") + ggtitle("") + ggdist::theme_ggdist() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
        axis.ticks = element_line(color="black"))  + rotate_x_text(0) 


VlnPlot(spleen.merged, features = "PROLIF1", pt.size = 0) +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95) + 
  theme(legend.position = "none") + ylab("Proliferation Gene Module Expression") + xlab("") + ggtitle("") + ggdist::theme_ggdist() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
        axis.ticks = element_line(color="black"))  + rotate_x_text(0) 

###### Extended Data Fig. 4i-j

table(spleen.merged@meta.data$cloneSize)

spleen.merged@meta.data$cloneSizeBin <- as.character(spleen.merged@meta.data$cloneSize)
spleen.merged@meta.data$cloneSizeBin[spleen.merged@meta.data$cloneSizeBin %in% c("Single (0 < X <= 1)",
                                                                                             "Small (1 < X <= 5)",
                                                                                             "Medium (5 < X <= 20)",
                                                                                             "Large (20 < X <= 50)")] <- "Small"


spleen.merged@meta.data$cloneSizeBin[spleen.merged@meta.data$cloneSizeBin %in% c("Hyperexpanded (50 < X <= 200)",
                                                                                             "Hyperexpanded (50 < X <= 649)",
                                                                                             "Hyperexpanded (50 < X <= 2078)")] <- "Large"


Idents(object = spleen.merged) <- "cloneSizeBin"

VlnPlot(spleen.merged, features = "MKI67", pt.size = 0, idents = c("Small","Large")) +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95) + 
  theme(legend.position = "none") + ylab("MKI67 Expression") + xlab("") + ggtitle("") + ggdist::theme_ggdist() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
        axis.ticks = element_line(color="black"))  + rotate_x_text(0) 

VlnPlot(spleen.merged, features = "PROLIF1", pt.size = 0, idents = c("Small","Large")) +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95) + 
  theme(legend.position = "none") + ylab("Proliferation Gene Module Expression") + xlab("") + ggtitle("") + ggdist::theme_ggdist() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
        axis.ticks = element_line(color="black"))  + rotate_x_text(0) 






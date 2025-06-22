library(Seurat)
library(ggplot2)
library(scRepertoire)
library(ggpubr)
library(viridis)

# Done: Ex. 5 a,b,c ; Ex 6. a,c,d,e,f,h
# To do: Ex. 6 b

spleen_brain_b_cells <- readRDS("/spleen_brain_B_cells_seurat_object.rds")

###### Ex Figure 5a

ebv_cells <- rownames(spleen_brain_b_cells@meta.data[spleen_brain_b_cells@meta.data$ebv_pbs == "EBV",])

col.clusters <- c("#EBE747","#1BC4D1","#43DB95","#E276BF","#E65658")

DimPlot(spleen_brain_b_cells, reduction = "tsne.rpca", label = FALSE, pt.size = 1.5, 
        label.size = 5, split.by = "organ", cells = ebv_cells) + 
  scale_color_manual(values = col.clusters) + xlab("tSNE 1") + ylab("tSNE 2")

###### Ex Figure 5b

Idents(object = spleen_brain_b_cells) <- "predicted.id"

ebv_cells <- rownames(spleen_brain_b_cells@meta.data[spleen_brain_b_cells@meta.data$ebv_pbs == "EBV",])

spleen_brain_b_cells@meta.data$cloneType <- as.factor(spleen_brain_b_cells@meta.data$cloneType)
spleen_brain_b_cells@meta.data$cloneType <- factor(spleen_brain_b_cells@meta.data$cloneType, levels = c("Single (0 < X <= 1)",
                                                                                                        "Small (1 < X <= 5)",
                                                                                                        "Medium (5 < X <= 20)",
                                                                                                        "Large (20 < X <= 50)",
                                                                                                        "Hyperexpanded (50 < X <= 200)"))

Idents(object = spleen_brain_b_cells) <- "cloneType"

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
colorblind_vector <- colorblind_vector[c(1:3,6:7)]

DimPlot(spleen_brain_b_cells, reduction = "tsne.rpca", pt.size = 1.5, 
        order = TRUE, split.by = "organ", cells = ebv_cells) +
  scale_color_manual(values = colorblind_vector) + 
  ylab("tSNE 2") + xlab("tSNE 1")

###### Ex Figure 5c & Ex Figure 6f

no.of.genes <- 20

gene_bcr_table <- rev(sort(table(spleen_brain_b_cells@meta.data[["CTgene"]])))
top_genes <- names(gene_bcr_table[c(1:no.of.genes)])
top_genes_counts <- as.numeric(gene_bcr_table[c(1:no.of.genes)])
gene_bcr_table[c(1:no.of.genes)]

df.for.heatmap <- spleen_brain_b_cells@meta.data[spleen_brain_b_cells@meta.data$CTgene %in% top_genes,c("organ","hash.ID","ebv_pbs","CTgene")]

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
heatmap_data_wide <- heatmap_data_wide[,c(2:25)]

heatmap_data_wide <- heatmap_data_wide[top_genes, , drop = FALSE]

hold.names <- rownames(heatmap_data_wide)

row_sums <- rowSums(heatmap_data_wide)

spleen_brain_b_cells@meta.data$hash_organ <- paste0(spleen_brain_b_cells@meta.data$hash.ID, "_", spleen_brain_b_cells@meta.data$organ)

hash_organ_total_counts <- table(spleen_brain_b_cells@meta.data$hash_organ)

hash_organ_total_counts_ebv <- hash_organ_total_counts[names(hash_organ_total_counts) %in% c("Hash9_brain","Hash9_spleen",
                                                                                             "Hash12_spleen",
                                                                                             "Hash14_brain","Hash14_spleen",
                                                                                             "Hash17_spleen",
                                                                                             "Hash2_brain","Hash2_spleen",
                                                                                             "Hash21_brain","Hash21_spleen",
                                                                                             "Hash23_brain","Hash23_spleen",
                                                                                             "Hash4_brain","Hash4_spleen",
                                                                                             "Hash13_brain","Hash13_spleen",
                                                                                             "Hash8_brain","Hash8_spleen")]


hash_organ_total_counts_pbs <- hash_organ_total_counts[names(hash_organ_total_counts) %in% c("Hash1_brain","Hash1_spleen",
                                                                                             "Hash6_brain","Hash6_spleen",
                                                                                             "Hash10_brain","Hash10_spleen")]

heatmap_data_wide_ebv <- heatmap_data_wide[,c("Hash9_brain","Hash9_spleen",
                                              "Hash12_spleen",
                                              "Hash14_brain","Hash14_spleen",
                                              "Hash17_spleen",
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

custom_groups <- as.numeric(c("9","9","12","12","14",
                              "14","17","17","2","2","21",
                              "21","23","23","4","4",
                              "13","13","8","8"))
custom_groups <- as.factor(custom_groups)

custom_groups <- factor(custom_groups, levels = c(9,12,14,17,2,21,23,4,13,8))

custom_groups_pbs <- as.numeric(c("1","1","6","6","10","10"))
custom_groups_pbs <- as.factor(custom_groups_pbs)
custom_groups_pbs <- factor(custom_groups_pbs, levels = c(1,6,10))

row.ordering <- c("Hash9_brain","Hash9_spleen",
                  "Hash12_spleen","Hash14_brain","Hash14_spleen","Hash17_spleen","Hash2_brain","Hash2_spleen","Hash21_brain","Hash21_spleen",
                  "Hash23_brain","Hash23_spleen","Hash4_brain","Hash4_spleen","Hash13_brain","Hash13_spleen","Hash8_brain","Hash8_spleen")

row.ordering.pbs <- c("Hash1_brain","Hash1_spleen","Hash6_brain","Hash6_spleen","Hash10_brain","Hash10_spleen")


heatmap_data_wide_ebv <- heatmap_data_wide_ebv[, row.ordering, drop = FALSE]
rownames(heatmap_data_wide_ebv) <- paste0("BCR",c(1:no.of.genes))

heatmap_data_wide_pbs <- heatmap_data_wide_pbs[, row.ordering.pbs, drop = FALSE]
rownames(heatmap_data_wide_pbs) <- paste0("BCR",c(1:no.of.genes))

library(ComplexHeatmap)

col_anno <- columnAnnotation(
  Counts = anno_barplot(row_sums_ebv, 
                        gp = gpar(fill = "black"),  # Bar color
                        border = FALSE)  # Remove borders
)

heatmap_data_wide_ebv_backup <- heatmap_data_wide_ebv
heatmap_data_wide_ebv <- heatmap_data_wide_ebv_backup


heatmap_data_wide_ebv <- log2(heatmap_data_wide_ebv + 1)


library(tibble)

heatmap_data_wide_ebv <- add_column(as.tibble(heatmap_data_wide_ebv), Hash12_brain = rep(0, nrow(heatmap_data_wide_ebv)), 
                                    .before = "Hash12_spleen")

heatmap_data_wide_ebv <- add_column(as.tibble(heatmap_data_wide_ebv), Hash17_brain = rep(0, nrow(heatmap_data_wide_ebv)), 
                                    .before = "Hash17_spleen")

ComplexHeatmap::Heatmap(t(heatmap_data_wide_ebv), 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        row_split = custom_groups,
                        
                        #left_annotation = row_anno,
                        col = colorRamp2(seq(0, 6, length.out = 100), viridis_colors),
                        rect_gp = gpar(col = "black", lwd = 1),
                        top_annotation = col_anno,
                        name = "log2(% BCR of sample cells)")





col_anno_pbs <- columnAnnotation(
  Counts = anno_barplot(row_sums_pbs, 
                        ylim = c(0, 400),
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



###### Ex Figure 6a

spleen_brain_b_cells@meta.data$predicted.id <- as.factor(spleen_brain_b_cells@meta.data$predicted.id)
spleen_brain_b_cells@meta.data$predicted.id <- factor(spleen_brain_b_cells@meta.data$predicted.id, levels = c("Naive B cells",
                                                                                                              "Memory B cells",
                                                                                                              "Cycling immune mix",
                                                                                                              "Plasmablasts",
                                                                                                              "Plasma cells"))
Idents(object = spleen_brain_b_cells) <- "predicted.id"

meta.data <- spleen_brain_b_cells@meta.data
meta.data <- meta.data[meta.data$ebv_pbs == "EBV",]

cell_abun <- as.data.frame(table(meta.data$predicted.id, meta.data$organ))

ggplot(cell_abun, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity", color = "black") + scale_fill_manual(values = c("#EBE747","#1BC4D1","#43DB95","#E276BF","#E65658")) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))

###### Ex Figure 6b

DefaultAssay(spleen_brain_b_cells) <- "RNA"

spleen_brain_b_cells <- JoinLayers(spleen_brain_b_cells)

# List genes
ebv.genes <- c("exon-V01555.2:58..272-1","exon-V01555.2:540..788-1",
               "exon-V01555.2:871..951-1","exon-V01555.2:1026..1196-1",
               "exon-V01555.2:1280..1495-1","exon-V01555.2:5408..5856-1",
               "exon-V01555.2:11336..11480-1","exon-V01555.2:11626..11657-1",
               "exon-V01555.2:47761..47793-1","exon-V01555.2:47878..47999-1",
               "exon-V01555.2:48386..48444-1","exon-V01555.2:49852..50032-1",
               "exon-V01555.2:67477..67649-1","exon-V01555.2:92238..92581-1",
               "exon-V01555.2:92670..95248-1","exon-V01555.2:98364..98730-1",
               "exon-V01555.2:98805..99050-1","exon-V01555.2:102655..103194-1",
               "exon-V01555.2:166498..166916-1","exon-V01555.2:168163..168965-1",
               "exon-V01555.2:169207..169474-1")

# Load table
list.merge.ebv <- readxl::read_xlsx('.../EBV_EXONS.xlsx')

# match genes with table
ebv.genes.match <- c()
ebv.genes.match$exons <- ebv.genes
ebv.genes.match$genes <- list.merge.ebv$Protein
ebv.genes.match <- as.data.frame(ebv.genes.match)

# Collapse into list by genes
ebv.genes.list <- aggregate(exons ~ genes, ebv.genes.match, paste, collapse = ", ")

# Isolate count matrix
counts_matrix <- GetAssayData(spleen_brain_b_cells, assay = "RNA")
counts_df <- as.data.frame(counts_matrix)

## Convert exons to gene by adding up the counts
for (i in 1:nrow(ebv.genes.list)) {
  gene <- ebv.genes.list$genes[i]
  exons <- unlist(strsplit(ebv.genes.list$exons[i], ", "))
  
  # Ensure exons are numeric and index within counts_df
  exon_indices <- match(exons, rownames(counts_df))
  exon_indices <- exon_indices[!is.na(exon_indices)] # Remove NA values
  
  # Subset to numeric columns if necessary
  numeric_columns <- sapply(counts_df, is.numeric)
  numeric_counts_df <- counts_df[, numeric_columns, drop = FALSE]
  
  # Summing the selected rows for numeric columns
  if (length(exon_indices) > 0) {
    gene.vector <- colSums(numeric_counts_df[exon_indices, , drop = FALSE])
    
    # Creating a data frame row from gene.vector
    gene_df <- as.data.frame(t(gene.vector), stringsAsFactors = FALSE)
    colnames(gene_df) <- names(gene.vector)
    rownames(gene_df) <- gene
    
    # Adding the summed row
    counts_df <- rbind(counts_df, gene_df)
  }
  
  # Removing the rows by exon indices
  counts_df <- counts_df[-exon_indices, ]
}

# Create new Seurat object
spleen_brain_b_cells.merged.exons <- CreateSeuratObject(counts = counts_df,  meta.data = spleen_brain_b_cells@meta.data)

# Append previous umap rpca reductions
spleen_brain_b_cells.merged.exons@reductions$umap.rpca <- spleen_brain_b_cells@reductions$tsne.rpca

# Normalize / Scale / Find Variable features / PCA
spleen_brain_b_cells.merged.exons <- SCTransform(spleen_brain_b_cells.merged.exons)

spleen_brain_b_cells.merged.exons <- AddModuleScore(spleen_brain_b_cells.merged.exons, 
                                                    features = list(c("EBNA","LMP1","LMP2A")),
                                                    name = "EBV_transcripts")


spleen_brain_b_cells.merged.exons@meta.data$organ <- as.factor(spleen_brain_b_cells.merged.exons@meta.data$organ)
spleen_brain_b_cells.merged.exons@meta.data$organ <- factor(spleen_brain_b_cells.merged.exons@meta.data$organ, levels = c("spleen",
                                                                                                                          "brain"))


FeaturePlot(spleen_brain_b_cells.merged.exons, features = "EBV_transcripts1", 
            order = T, max.cutoff = 'q97', split.by = "organ", 
            reduction = "umap.rpca", ncol = 3, pt.size = 1.5, cells = ebv_cells) & 
  scale_color_viridis(option = "D") & theme(legend.position = "right")




###### Ex Figure 6c

cell_abun <- as.data.frame(table(spleen_brain_b_cells@active.ident, spleen_brain_b_cells@meta.data$organ, spleen_brain_b_cells@meta.data$predicted.id, spleen_brain_b_cells@meta.data$ebv_pbs))
cell_abun <- cell_abun[cell_abun$Var4 %in% c("EBV"),]
cell_abun <- cell_abun[cell_abun$Var3 %in% c("Memory B cells",
                                             "Cycling immune mix",
                                             "Plasmablasts",
                                             "Plasma cells"),]

cell_abun <- cell_abun[,-c(3,4)]

library(dplyr)

df_combined <- cell_abun %>%
  group_by(Var1, Var2) %>%
  summarise(total_amount = sum(Freq))

df_combined$Var1 <- factor(df_combined$Var1, levels = rev(levels(df_combined$Var1)))

ggplot(df_combined, aes(fill=Var1, y=total_amount, x=Var2)) + 
  geom_bar(position="fill", stat="identity", color = "black") + scale_fill_manual(values = rev(c(colorblind_vector))) + xlab("") + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Frequency") + ggdist::theme_ggdist() + theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))

###### Ex Figure 6d

pbs_cells <- rownames(spleen_brain_b_cells@meta.data[spleen_brain_b_cells@meta.data$ebv_pbs == "PBS",])

col.clusters <- c("#EBE747","#1BC4D1","#43DB95","#E276BF","#E65658")

DimPlot(spleen_brain_b_cells, reduction = "tsne.rpca", label = FALSE, pt.size = 1.5, 
        label.size = 5, split.by = "organ", cells = pbs_cells) + 
  scale_color_manual(values = col.clusters) + xlab("tSNE 1") + ylab("tSNE 2")

###### Ex Figure 6e

Idents(object = spleen_brain_b_cells) <- "predicted.id"

spleen_brain_b_cells@meta.data$cloneType <- as.factor(spleen_brain_b_cells@meta.data$cloneType)
spleen_brain_b_cells@meta.data$cloneType <- factor(spleen_brain_b_cells@meta.data$cloneType, levels = c("Single (0 < X <= 1)",
                                                                                                        "Small (1 < X <= 5)",
                                                                                                        "Medium (5 < X <= 20)",
                                                                                                        "Large (20 < X <= 50)",
                                                                                                        "Hyperexpanded (50 < X <= 200)"))

Idents(object = spleen_brain_b_cells) <- "cloneType"

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
colorblind_vector <- colorblind_vector[c(1:3,6:7)]

pbs_cells <- rownames(spleen_brain_b_cells@meta.data[spleen_brain_b_cells@meta.data$ebv_pbs == "PBS",])

DimPlot(spleen_brain_b_cells, reduction = "tsne.rpca", pt.size = 1.5, 
        order = TRUE, split.by = "organ", cells = pbs_cells) +
  scale_color_manual(values = colorblind_vector) + 
  ylab("tSNE 2") + xlab("tSNE 1")

###### Ex Figure 6g

Idents(object = spleen_brain_b_cells) <- "ebv_pbs"

spleen_brain_b_cells_ebv <- subset(spleen_brain_b_cells, idents = "EBV")

Idents(spleen_brain_b_cells_ebv) <- "predicted.id"

VlnPlot(
  spleen_brain_b_cells_ebv,
  features = c("TNFSF9"),
  pt.size = 0, split.by = "organ", idents = "Memory B cells",
) +
  RestoreLegend() + xlab("") + ylab("Expression") + rotate_x_text(0)

###### Ex Figure 6h

heavy_gene_cleaning <- sapply(strsplit(spleen_brain_b_cells@meta.data[["CTgene"]], "[.]"), `[`, 1)
heavy_gene_cleaning_2 <- sapply(strsplit(heavy_gene_cleaning, "_"), `[`, 1)

heavy_gene <- data.frame(heavy_gene = heavy_gene_cleaning_2,
                         organ = spleen_brain_b_cells@meta.data$organ,
                         ebv_pbs = spleen_brain_b_cells@meta.data$ebv_pbs,
                         sample = spleen_brain_b_cells@meta.data$hash.ID)

heavy_gene <- heavy_gene[heavy_gene$ebv_pbs == "EBV",]
heavy_gene <- heavy_gene[!(heavy_gene$heavy_gene == "NA"),]

cell_abun <- as.data.frame(table(heavy_gene$heavy_gene, heavy_gene$organ))

library(tidyr)
library(dplyr)

df_percent <- cell_abun %>%
  group_by(Var2) %>%
  mutate(percent = Freq / sum(Freq) * 100) %>%
  ungroup()

df_wide_percent <- df_percent %>%
  dplyr::select(Var1, Var2, percent) %>%
  pivot_wider(names_from = Var2, values_from = percent)


ggscatter(df_wide_percent, x = "spleen", y = "brain",
          add = "loess", conf.int = TRUE)




# Heavy gene usage in brain vs. spleen - statistics including replicates

cell_abun <- as.data.frame(table(heavy_gene$heavy_gene, heavy_gene$organ, heavy_gene$sample ))
cell_abun <- cell_abun[!(cell_abun$Var3 == "Hash12" & cell_abun$Var2 == "brain"),]

df_percent <- cell_abun %>%
  group_by(Var3, Var2) %>%   # Group by sample and tissue
  mutate(Percent = Freq / sum(Freq) * 100) %>%  # Normalize per tissue per sample
  ungroup()

a <- df_percent[df_percent$Var1 == "IGHV4-39",]

df_mean_percent <- df_percent %>%
  group_by(Var1, Var2) %>%
  summarise(MeanPercent = mean(Percent), .groups = "drop") %>%
  pivot_wider(names_from = Var2, values_from = MeanPercent)



rank_sum_results <- df_percent %>%
  dplyr::group_by(Var1) %>%
  dplyr::filter(Var2 %in% c("brain", "spleen")) %>%
  dplyr::summarise(
    rank_sum_test = list(wilcox.test(Percent ~ Var2, exact = FALSE)),  # Mann-Whitney U test (Wilcoxon rank-sum)
    .groups = "drop"
  ) %>%
  mutate(
    P_value = sapply(rank_sum_test, function(x) x$p.value),  # Extract p-value
    Adj_P = p.adjust(P_value, method = "holm"),  # Apply Holm-Sidak correction for multiple comparisons
    Enriched = case_when(
      Adj_P < 0.05 & rank_sum_test[[1]]$statistic[1] > rank_sum_test[[1]]$statistic[2] ~ "brain",  # Enriched in brain
      Adj_P < 0.05 & rank_sum_test[[1]]$statistic[1] < rank_sum_test[[1]]$statistic[2] ~ "spleen",  # Enriched in spleen
      TRUE ~ NA_character_  # Not enriched
    )
  )

df_plot <- df_mean_percent %>%
  left_join(rank_sum_results %>% dplyr::select(Var1, Enriched, Adj_P), by = "Var1")

specific_genes <- c("IGHV4-34","IGHV4-4","IGHV3-7","IGHV1-69D")

df_plot <- df_plot %>%
  mutate(
    Var1 = as.character(Var1),
    Label_Gene = ifelse(Var1 %in% specific_genes, Var1, NA),  # Label only specific genes
    Color_Gene = ifelse(Var1 %in% specific_genes, "highlight", "default")  # Color specific genes differently
  )

df_plot$brain_spleen_ratio <- df_plot$brain / df_plot$spleen


ggplot(df_plot, aes(x = spleen, y = brain)) +
  geom_point(aes(fill = Color_Gene), size = 3, color = "black", shape = 21, stroke = 1) +  # Points with color based on enrichment
  geom_text(data = dplyr::filter(df_plot, !is.na(Label_Gene)), aes(label = Label_Gene), vjust = -0.5, size = 5) +  # Add labels for enriched genes
  scale_fill_manual(values = c("highlight" = "red", "default" = "black")) +  # Coloring
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +  # Linear regression line with CI
  labs(
    x = "Spleen (%)",
    y = "Brain (%)",
    color = "Enriched in"
  ) + ggdist::theme_ggdist() + theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
    axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
    axis.ticks = element_line(color="black"))












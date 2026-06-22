#Import libraries
library(svglite)
library(tidyverse)
library(Seurat)
library(EnhancedVolcano)
library(pheatmap)
library(tibble)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(grid)

##################################
desired_order <- c(
  "Choroid Plexus", "Vascular cells", "Microglia", "Astrocytes", "Oligodendrocytes",
  "Oligo Progenitor Cells", "DG-PIR Immature Neurons", "Dentate Gyrus Neurons", "CTX-CGE GABA",
  "VIP GABA Neurons", "Lamp5 GABA Neurons", "Pvalb GABA Neurons", "Sst GABA Neurons", "NP CT L6b Glut",
  "IT-ET Glut Neurons", "SUB-ProS Glut Neurons", "HPF CR Glut Neurons", "CA3 Glut Neurons",
  "CA2 Glut Neurons", "CA1 Glut Neurons"
)

subtype_palette <- c(
  "Choroid Plexus"          = "#D08A3C",
  "Vascular cells"          = "#5FB7C3",
  "Microglia"               = "#7FBF7B",
  "Astrocytes"              = "#E6A07A",
  "Oligodendrocytes"        = "#6F9FD8",
  "Oligo Progenitor Cells"  = "#D48AC2",
  "DG-PIR Immature Neurons" = "#F2C94C",
  "Dentate Gyrus Neurons"   = "#9B85C9",
  "CTX-CGE GABA"            = "#E58E73",
  "Immature GABA"           = "#7FC8A9",
  "VIP GABA Neurons"        = "#F0B45B",
  "Lamp5 GABA Neurons"      = "#6EC6B2",
  "Pvalb GABA Neurons"      = "#6FA8DC",
  "Sst GABA Neurons"        = "#E58A9A",
  "NP CT L6b Glut"          = "#C4A7D9",
  "IT-ET Glut Neurons"      = "#E9CFA2",
  "SUB-ProS Glut Neurons"   = "#8BB8E8",
  "HPF CR Glut Neurons"     = "#9AD9C5",
  "CA3 Glut Neurons"        = "#E07B73",
  "CA2 Glut Neurons"        = "#F0A3A3",
  "CA1 Glut Neurons"        = "#9FC97A"
)
########################################################

#Fig_2a

DimPlot(KA_seurat_plots,
        group.by = "condition",
        label = FALSE,
        cols = c("Control" = "royalblue1",  # pastel blue
                 "Epi"  = "red1")  # pastel red/pink
) +
  scale_color_manual(
    values = c("Control" = "royalblue1", "Epi" = "red1"),
    breaks = c("Control", "Epi")  # Controls legend order
  ) +
  theme(legend.text = element_text(size = 14))

#Fig_2b

metadata <- KA_seurat_plots@meta.data

cell_counts <- metadata %>%
  filter(classes %in% desired_order) %>%
  group_by(sample.origin, condition, classes) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  mutate(
    condition = factor(condition, levels = c("Control", "Epi")),
    classes   = factor(classes, levels = desired_order)
  )

ggplot(cell_counts, aes(classes, cell_count, fill = condition)) +
  geom_boxplot(
    position = position_dodge(0.8),
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(fill = condition),
    shape = 21,
    color = "black",
    stroke = 0.4,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.8
    ),
    size = 2,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(Control = "royalblue1", Epi = "red1"),
    breaks = c("Control","Epi")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    legend.text = element_text(size = 16),
    axis.title  = element_text(size = 16)
  ) +
  labs(
    x = "Cell subtype",
    y = "Number of Cells",
    fill = "Condition",
    title = "Cell Counts per Subtype and Condition"
  )

#Fig_2c

padj_cutoff <- 0.05
findMarkers_results_full <- de_result_classes_pseudo_DESeq2_list


findMarkers_results <- lapply(findMarkers_results_full, function(df) {
  df <- as.data.frame(df)
  df[df$p_adj_BH < padj_cutoff, , drop = FALSE]
})


findMarkers_results <- findMarkers_results[sapply(findMarkers_results, nrow) > 0]

deg_count <- data.frame(
  subtype    = names(findMarkers_results),
  DEG_number = sapply(findMarkers_results, nrow),
  stringsAsFactors = FALSE
)

deg_count <- deg_count[order(deg_count$DEG_number, decreasing = TRUE), ]
top_subtypes <- deg_count[seq_len(min(10, nrow(deg_count))), ]

findMarkers_results_top <- findMarkers_results[top_subtypes$subtype]

data <- bind_rows(lapply(names(findMarkers_results_top5), function(sub) {
  df <- as.data.frame(findMarkers_results_top5[[sub]])
  df$subtype <- sub
  df
}))

top_subtypes$subtype <- factor(top_subtypes$subtype, levels = top_subtypes$subtype)
data$subtype <- factor(data$subtype, levels = levels(top5_subtypes$subtype))

color_mapping <- subtype_palette[levels(top_subtypes$subtype)]

bar_plot <- ggplot(top_subtypes, aes(x = subtype, y = DEG_number, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = DEG_number), vjust = -0.4, size = 4) +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal(base_size = 15) +
  labs(x = NULL, y = "DEG Count") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

scatter_plot <- ggplot(data, aes(x = subtype, y = log2FC, color = subtype)) +
  geom_jitter(alpha = 0.7, size = 1, width = 0.2) +
  scale_color_manual(values = color_mapping) +
  theme_minimal(base_size = 15) +
  labs(x = NULL, y = "log2 Fold Change") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  )


fig_2c <- bar_plot / scatter_plot
print(fig_2c)


#Fig_2d

top_n <- 12
padj_cutoff <- 0.05
baseMean_cutoff <- 20
log2fc_cap <- 2

selected_subtypes <- c(
  "CA1 Glut Neurons",
  "CA3 Glut Neurons",
  "Dentate Gyrus Neurons",
  "NP CT L6b Glut",
  "Astrocytes"
)

de_subset <- de_result_classes_pseudo_DESeq2_list[selected_subtypes]

# Select top genes per subtype
top_genes_by_subtype <- lapply(de_subset, function(df) {

  df <- df[
    !is.na(df$p_adj_BH) &
      !is.na(df$log2FC) &
      !is.na(df$baseMean),
    ,
    drop = FALSE
  ]

  df_filtered <- df[
    df$p_adj_BH < padj_cutoff &
      df$baseMean > baseMean_cutoff,
    ,
    drop = FALSE
  ]

  df_filtered <- df_filtered[order(df_filtered$p_adj_BH), , drop = FALSE]

  head(rownames(df_filtered), top_n)
})

selected_genes <- unique(unlist(top_genes_by_subtype))

# Build log2FC matrix
log2fc_mat <- sapply(selected_subtypes, function(subtype) {

  df <- de_subset[[subtype]]

  lfc <- df[selected_genes, "log2FC"]
  base_mean <- df[selected_genes, "baseMean"]

  lfc[is.na(lfc)] <- 0
  base_mean[is.na(base_mean)] <- 0

  # Set low-expression genes to 0 so they do not drive the heatmap
  lfc[base_mean <= baseMean_cutoff] <- 0

  lfc
})

rownames(log2fc_mat) <- selected_genes

# Build significance matrix for asterisk labels
sig_mat <- sapply(selected_subtypes, function(subtype) {

  df <- de_subset[[subtype]]

  padj <- df[selected_genes, "p_adj_BH"]
  base_mean <- df[selected_genes, "baseMean"]

  padj[is.na(padj)] <- 1
  base_mean[is.na(base_mean)] <- 0

  ifelse(
    padj < padj_cutoff & base_mean > baseMean_cutoff,
    "*",
    ""
  )
})

rownames(sig_mat) <- selected_genes

# Cap extreme log2FC values for plotting
log2fc_mat[log2fc_mat > log2fc_cap] <- log2fc_cap
log2fc_mat[log2fc_mat < -log2fc_cap] <- -log2fc_cap

# Order genes by best adjusted p-value across selected subtypes
padj_mat <- sapply(selected_subtypes, function(subtype) {

  df <- de_subset[[subtype]]

  padj <- df[selected_genes, "p_adj_BH"]
  padj[is.na(padj)] <- 1

  padj
})

rownames(padj_mat) <- selected_genes

gene_order <- order(apply(padj_mat, 1, min))

log2fc_mat <- log2fc_mat[gene_order, , drop = FALSE]
sig_mat <- sig_mat[gene_order, , drop = FALSE]

# Column annotations
region_colors <- subtype_palette[selected_subtypes]

class_colors <- c(
  "Neuron" = "darkgoldenrod1",
  "Glia" = "lightblue"
)

column_annotation <- HeatmapAnnotation(
  Class = c(rep("Neuron", 4), "Glia"),
  Region = colnames(log2fc_mat),
  col = list(
    Class = class_colors,
    Region = region_colors
  ),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(4, "mm")
)

# Heatmap colour scale
col_fun <- colorRamp2(
  c(-log2fc_cap, 0, log2fc_cap),
  c("#3B4CC0", "#FFFFFF", "#B40426")
)

# Plot heatmap
ht <- Heatmap(
  log2fc_mat,
  name = "log2FC",
  top_annotation = column_annotation,
  col = col_fun,

  cluster_rows = TRUE,
  cluster_columns = TRUE,

  column_split = c(
    rep("Hippocampus/Cortex Neurons", 4),
    "Astrocytes"
  ),

  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 11),
  row_names_gp = gpar(fontsize = 11),

  rect_gp = gpar(col = NA),

  show_row_dend = TRUE,
  show_column_dend = TRUE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] == "*") {
      grid.text(
        "*",
        x,
        y,
        gp = gpar(
          fontsize = 10,
          fontface = "bold",
          col = "black"
        )
      )
    }
  },

  heatmap_legend_param = list(
    at = c(-log2fc_cap, -1, 0, 1, log2fc_cap),
    labels = c(
      paste0("≤-", log2fc_cap),
      "-1",
      "0",
      "1",
      paste0("≥", log2fc_cap)
    )
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)


#Fig_2e

de_result_classes_pseudo_DESeq2_list[["Dentate Gyrus Neurons"]]$gene<- rownames(de_result_classes_pseudo_DESeq2_list$`Dentate Gyrus Neurons`)

dg_col <- subtype_palette[["Dentate Gyrus Neurons"]]

EnhancedVolcano(
  de_result_classes_pseudo_DESeq2_list[["Dentate Gyrus Neurons"]],
  lab = de_result_classes_pseudo_DESeq2_list[["Dentate Gyrus Neurons"]]$gene,
  x = "log2FC",
  y = "p_adj_BH",
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p-value"),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  col = c("grey80", "grey80", "grey80", dg_col),
  pointSize = 2.0,
  labSize = 5.0,
  title = "DEGs in Dentate Gyrus Neurons",
  subtitle = ""
)




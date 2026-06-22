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


#Fig_1b
# Extract and prepare data
cell_counts <- KA_seurat_plots@meta.data %>%
  dplyr::count(sample.origin, condition)

# Set factor levels to control x-axis and legend order
cell_counts$condition <- factor(cell_counts$condition, levels = c("Control", "Epi"))
cell_counts$sample.origin <- factor(cell_counts$sample.origin,
                                    levels = cell_counts$sample.origin[cell_counts$condition == "Control"] %>%
                                      unique() %>%
                                      c(., unique(cell_counts$sample.origin[cell_counts$condition == "Epi"])))

# Plot

ggplot(cell_counts, aes(x = condition, y = n, color = condition, fill = condition)) +
  geom_boxplot(alpha = 0.2, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "royalblue1", "Epi" = "red1")) +
  scale_fill_manual(values = c("Control" = "royalblue1", "Epi" = "red1")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal() +
  ylab("Number of nuclei") +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 18)
  )


#Fig_1c

DimPlot(
  KA_seurat_plots,
  group.by = "sample.origin",
  pt.size = 0.01,
  cols = c(
    "Epi1" = "#E69F00",
    "Epi2" = "#0072B2",
    "Epi3" = "#F0E442",
    "Epi4" = "#56B4E9",
    "Control1" = "#D55E00",
    "Control2" = "#CC79A7",
    "Control3" = "#009E73",
    "Control4" = "#BCE9C5"
  )
) +
  theme(
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 18),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16)
  )

#Fig_1d

  DimPlot(
    KA_seurat_plots,
    group.by = "classes",
    label = TRUE,
    repel = TRUE,
    label.size = 4.5,
    pt.size = 0.01,
    cols = subtype_palette
  ) +
    NoLegend()


#Fig_1e

  DotPlot(KA_seurat_plots, features = c("Slc17a7","Trp73","Pou3f1","Tspan18","Satb2","Foxp2","Gad1",
    "Sst","Pvalb","Lamp5","Vip","Prox1","Mex3a","Pdgfra","Mbp","Gfap","C1qb","Aqp4","Apod","Ttr"), group.by = "classes")
  +scale_y_discrete(limits=desired_order) + theme(axis.text.x = element_text(angle = 45, hjust = 1))



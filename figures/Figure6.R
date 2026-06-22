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
# Figure 6e

df_long <- latrophilin_validation %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Condition",
    values_to = "Value"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = latrophilin_validation$Gene),
    Condition = factor(
      Condition,
      levels = colnames(latrophilin_validation)[-1]
    )
  )

fig_6e <- ggplot(df_long, aes(x = Gene, y = Condition, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C23B3B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("≤ -2", "-1", "0", "1", "≥ 2"),
    na.value = "grey92",
    name = "LogFC"
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = "Latrophilin pathway validation"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

fig_6e
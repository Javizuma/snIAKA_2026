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

#Fig_5b

ipa_astrocytes_pseudo$`z-score` <- as.numeric(ipa_astrocytes_pseudo$`z-score`)
ipa_astrocytes_pseudo$`-log(p-value)` <- as.numeric(ipa_astrocytes_pseudo$`-log(p-value)`)

top_df <- ipa_astrocytes_pseudo %>%
  arrange(desc(`-log(p-value)`)) %>%
  slice_head(n = 25)

top_df$`Ingenuity Canonical Pathways` <- factor(
  top_df$`Ingenuity Canonical Pathways`,
  levels = top_df$`Ingenuity Canonical Pathways`[order(top_df$`-log(p-value)`)]
)

ggplot(
  top_df,
  aes(
    x = `-log(p-value)`,
    y = `Ingenuity Canonical Pathways`,
    size = Ratio,
    fill = `z-score`
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.3,
    alpha = 0.9
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey85",
    name = "Z-score"
  ) +
  scale_size_continuous(
    name = "Ratio",
    range = c(3, 10)
  ) +
  labs(
    x = "-log(p-value)",
    y = "",
    title = "Astrocytes"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13),
    panel.grid.major.y = element_blank()
  )


#Fig_5c

library(dplyr)
library(forcats)
library(ggplot2)

path_summary <- paths_merged %>%
  mutate(
    is_changed = change_category != "unchanged",
    abs_FC = abs(FC)
  ) %>%
  filter(is_changed) %>%
  group_by(pathway_name) %>%
  summarise(
    altered_pairs = n(),
    mean_lprob = mean(FC, na.rm = TRUE),
    mean_abs_lprob = mean(abs_FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_lprob), desc(altered_pairs)) %>%
  slice_head(n = 15) %>%
  mutate(
    pathway_name = fct_reorder(pathway_name, mean_lprob)
  )

fig_5c <- ggplot(
  path_summary,
  aes(
    x = mean_lprob,
    y = pathway_name,
    size = altered_pairs,
    colour = mean_lprob
  )
) +
  geom_point(alpha = 0.9) +
  scale_colour_gradient2(
    name = "Mean prob",
    low = "blue",
    mid = "#F5E8DF",
    high = "red",
    midpoint = 0
  ) +
  scale_size_continuous(
    name = "Altered L-R pairs",
    range = c(2, 10)
  ) +
  labs(
    x = "Mean |.prob|",
    y = NULL,
    title = "Cell-cell communication changed in epileptic astrocytes"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

fig_5c


#Fig_5d

genes <- c("Gpam","Fdtd1","Cyp51","Idi1","Hmgcs1")
get_counts_for_gene <- function(g) {
  df <- plotCounts(
    dds_classes_pseudo_list[["Astrocytes"]],
    gene       = g,
    intgroup   = "group",   # assumes 'group' has Control / Epi
    returnData = TRUE
  )
  df$gene <- g
  df
}

pc <- bind_rows(lapply(genes, get_counts_for_gene))

pc$gene <- factor(pc$gene, levels = genes)

ggplot(pc, aes(x = gene, y = count, fill = group)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.5,
    color = "black",
    outlier.shape = NA,
    lwd = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(group = group),
    size = 2.2,
    alpha = 0.8,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.7
    )
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Expression of Cholesterol biosynthesis genes",
    subtitle= "In Astrocytes",
    x = "Gene",
    y = "Normalized count",
    fill = "Group"
  ) +
  theme_bw(base_size = 16)



#Fig_5e

genes <- c("Ank3","Chl1","Dnm3","Itga5","Nrp2")
get_counts_for_gene <- function(g) {
  df <- plotCounts(
    dds_classes_pseudo_list[["Astrocytes"]],
    gene       = g,
    intgroup   = "group",
    returnData = TRUE
  )
  df$gene <- g
  df
}

pc <- bind_rows(lapply(genes, get_counts_for_gene))

pc$gene <- factor(pc$gene, levels = genes)

ggplot(pc, aes(x = gene, y = count, fill = group)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.5,
    color = "black",
    outlier.shape = NA,
    lwd = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(group = group),
    size = 2.2,
    alpha = 0.8,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.7
    )
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Expression of L1cam interaction genes",
    subtitle= "In Astrocytes",
    x = "Gene",
    y = "Normalized count",
    fill = "Group"
  ) +
  theme_bw(base_size = 16)

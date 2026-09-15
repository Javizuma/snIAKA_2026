#Import libraries
library(svglite)
library(tidyverse)
library(Seurat)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Note: final figure assembly, panel labels and some annotations
# (e.g. significance asterisks and dataset group brackets) were added manually.

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
# Figure 6a & 6b

adgrl_lr <- pairs_merged %>%
  filter(change_category != "unchanged") %>%
  filter(pathway_name_PBS == "ADGRL" | pathway_name_KA == "ADGRL")

adgrl_lr$interaction_name_merged <- ifelse(
  !is.na(adgrl_lr$interaction_name_2_PBS),
  adgrl_lr$interaction_name_2_PBS,
  adgrl_lr$interaction_name_2_KA
)

df <- adgrl_lr

############################################################
# Figure 6a

counts <- df %>%
  mutate(
    interaction_name_merged = fct_explicit_na(interaction_name_merged, "(NA)"),
    target = fct_explicit_na(target, "(NA)")
  ) %>%
  dplyr::count(interaction_name_merged, target, name = "freq")

row_order <- counts %>%
  group_by(interaction_name_merged) %>%
  summarise(total = sum(freq), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(interaction_name_merged)

col_order <- counts %>%
  group_by(target) %>%
  summarise(total = sum(freq), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(target)

counts_ord <- counts %>%
  mutate(
    interaction_name_merged = factor(interaction_name_merged, levels = row_order),
    target = factor(target, levels = col_order)
  )

top_n_rows <- 40

counts_plot <- counts_ord %>%
  filter(as.integer(interaction_name_merged) <= min(top_n_rows, length(row_order)))

fig_6a <- ggplot(counts_plot, aes(x = target, y = interaction_name_merged, fill = freq)) +
  geom_tile(color = "black", linewidth = 0.05) +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(
    x = "Cell subtype",
    y = "L-R Pair",
    fill = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

fig_6a

########################################################
# Figure 6b

sums <- df %>%
  mutate(
    interaction_name_merged = fct_explicit_na(interaction_name_merged, "(NA)"),
    target = fct_explicit_na(target, "(NA)"),
    FC_num = suppressWarnings(as.numeric(FC))
  ) %>%
  filter(!is.na(FC_num)) %>%
  group_by(interaction_name_merged, target) %>%
  summarise(sum_abs_FC = sum(abs(FC_num)), .groups = "drop")

row_order <- sums %>%
  group_by(interaction_name_merged) %>%
  summarise(total = sum(sum_abs_FC), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(interaction_name_merged)

col_order <- sums %>%
  group_by(target) %>%
  summarise(total = sum(sum_abs_FC), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(target)

sums_ord <- sums %>%
  mutate(
    interaction_name_merged = factor(interaction_name_merged, levels = row_order),
    target = factor(target, levels = col_order)
  )

top_n_rows <- 40
keep_rows <- head(row_order, min(top_n_rows, length(row_order)))

sums_plot <- sums_ord %>%
  filter(interaction_name_merged %in% keep_rows)

fig_6b <- ggplot(sums_plot, aes(x = target, y = interaction_name_merged, fill = sum_abs_FC)) +
  geom_tile(color = "black", linewidth = 0.05) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Cell subtype",
    y = "L-R Pairs",
    fill = "|FC|"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

fig_6b

########################################################
# Figure 6c

genes <- c("Adgrl1", "Adgrl2", "Adgrl3", "Adgrl4")

get_counts_for_gene <- function(g) {
  df <- plotCounts(
    dds_classes_pseudo_list[["Dentate Gyrus Neurons"]],
    gene = g,
    intgroup = "group",
    returnData = TRUE
  )
  df$gene <- g
  df
}

pc <- bind_rows(lapply(genes, get_counts_for_gene))

pc$gene <- factor(pc$gene, levels = genes)

fig_6c <- ggplot(pc, aes(x = gene, y = count, fill = group)) +
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
      dodge.width = 0.7
    )
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Expression of the latrophilin family",
    subtitle = "In Dentate Gyrus neurons",
    x = "Gene",
    y = "Normalized count",
    fill = "Group"
  ) +
  theme_bw(base_size = 16)

fig_6c

########################################################
# Figure 6d

# Lower-expression ligands
genes <- c(
  "Tenm1", "Tenm3", "Tenm4", "Nrxn2",
  "Unc5a", "Unc5b", "Unc5c", "Unc5d", "Bdnf"
)

pc <- bind_rows(lapply(genes, get_counts_for_gene))

pc$gene <- factor(pc$gene, levels = genes)

fig_6d_left <- ggplot(pc, aes(x = gene, y = count, fill = group)) +
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
      dodge.width = 0.7
    )
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Expression of the latrophilin family ligands",
    subtitle = "In Dentate Gyrus neurons",
    x = "Gene",
    y = "Normalized count",
    fill = "Group"
  ) +
  theme_bw(base_size = 16)

fig_6d_left


# Higher-expression ligands
genes <- c("Tenm2", "Nrxn1", "Nrxn3")

pc <- bind_rows(lapply(genes, get_counts_for_gene))

pc$gene <- factor(pc$gene, levels = genes)

fig_6d_right <- ggplot(pc, aes(x = gene, y = count, fill = group)) +
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
      dodge.width = 0.7
    )
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Gene",
    y = "Normalized count",
    fill = "Group"
  ) +
  theme_bw(base_size = 16)

fig_6d_right


########################################################
# Figure 6e
# RT-qPCR validation added separately


########################################################
# Figure 6f

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

fig_6f <- ggplot(df_long, aes(x = Gene, y = Condition, fill = Value)) +
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
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

fig_6f
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
# Figure 3a

ipa_dentate_pseudo <- read_csv("../IPA_results/DG/Canonical_pathways.csv")
ipa_ca1_pseudo <- read_csv("../IPA_results/CA1 neurons/Canonical_pathways.csv")
ipa_ca3_pseudo <- read_csv("../IPA_results/CA3 neurons/Canonical_pathways.csv")
ipa_NP_pseudo <- read_csv("../IPA_results/NP CT L6b Glut/canonical_pathways.csv")
ipa_astrocytes_pseudo <- read_csv("../IPA_results/Astrocytes/Canonical_pathways.csv")

df_astrocytes <- ipa_astrocytes_pseudo[, c("Ingenuity Canonical Pathways", "-log(p-value)")]
names(df_astrocytes)[2] <- "Astrocytes"

df_ca1 <- ipa_ca1_pseudo[, c("Ingenuity Canonical Pathways", "-log(p-value)")]
names(df_ca1)[2] <- "CA1 Glut Neurons"

df_ca3 <- ipa_ca3_pseudo[, c("Ingenuity Canonical Pathways", "-log(p-value)")]
names(df_ca3)[2] <- "CA3 Glut Neurons"

df_dentate <- ipa_dentate_pseudo[, c("Ingenuity Canonical Pathways", "-log(p-value)")]
names(df_dentate)[2] <- "Dentate Gyrus Neurons"

df_NP <- ipa_NP_pseudo[, c("Ingenuity Canonical Pathways", "-log(p-value)")]
names(df_NP)[2] <- "NP CT L6b Glut"

merged_df <- Reduce(
  function(x, y) merge(x, y, by = "Ingenuity Canonical Pathways", all = TRUE),
  list(df_astrocytes, df_ca1, df_ca3, df_dentate, df_NP)
)

threshold <- 1.3 # equivalent to p < 0.05 as -log(p-value)

binary_df <- merged_df %>%
  mutate(across(
    c(
      `Astrocytes`,
      `CA1 Glut Neurons`,
      `CA3 Glut Neurons`,
      `Dentate Gyrus Neurons`,
      `NP CT L6b Glut`
    ),
    ~ as.integer(. > threshold)
  ))

binary_df_upset <- binary_df
rownames(binary_df_upset) <- binary_df$`Ingenuity Canonical Pathways`
binary_df_upset$`Ingenuity Canonical Pathways` <- NULL
binary_df_upset[is.na(binary_df_upset)] <- 0

upset(
  binary_df_upset,
  sets = colnames(binary_df_upset),
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "darkred",
  keep.order = TRUE,
  text.scale = c(2, 2, 2, 2, 2, 2)
)

# Figure 3b

ipa_dentate_pseudo$`z-score` <- as.numeric(as.character(ipa_dentate_pseudo$`z-score`))
ipa_dentate_pseudo$`-log(p-value)` <- as.numeric(as.character(ipa_dentate_pseudo$`-log(p-value)`))
ipa_dentate_pseudo$Ratio <- as.numeric(as.character(ipa_dentate_pseudo$Ratio))

top_df <- ipa_dentate_pseudo %>%
  arrange(desc(`-log(p-value)`)) %>%
  slice_head(n = 25)

top_df$`Ingenuity Canonical Pathways` <- factor(
  top_df$`Ingenuity Canonical Pathways`,
  levels = top_df$`Ingenuity Canonical Pathways`[order(top_df$`-log(p-value)`)]
)

blue_labels <- c(
  "Activation of NMDA receptors and postsynaptic events",
  "Glutaminergic Receptor Signaling Pathway (Enhanced)",
  "CREB Signaling in Neurons",
  "Calcium Signaling",
  "DAG and IP3 signaling",
  "Effects of PIP2 hydrolysis"
)

purple_labels <- c(
  "Synaptogenesis Signaling Pathway",
  "Extracellular matrix organization",
  "Axonal Guidance Signaling",
  "Regulation of CDH19 Expression and Function"
)

green_labels <- c(
  "Th1 and Th2 Activation Pathway",
  "Th2 Pathway",
  "IL-10 Signaling"
)

red_labels <- c(
  "Wnt/Ca+ pathway",
  "GABA Receptor Signaling"
)

label_color_map <- c(
  setNames(rep("blue", length(blue_labels)), blue_labels),
  setNames(rep("purple", length(purple_labels)), purple_labels),
  setNames(rep("#2C792D", length(green_labels)), green_labels),
  setNames(rep("darkred", length(red_labels)), red_labels)
)

label_fun <- function(x) {
  cols <- ifelse(x %in% names(label_color_map), label_color_map[x], "black")
  paste0("<span style='color:", cols, "'>", x, "</span>")
}

fig_3b <- ggplot(
  top_df,
  aes(
    x = `-log(p-value)`,
    y = `Ingenuity Canonical Pathways`,
    size = Ratio,
    fill = `z-score`
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey85",
    name = "Z-score"
  ) +
  scale_size_continuous(name = "Ratio", range = c(3, 10)) +
  scale_y_discrete(labels = label_fun) +
  labs(
    x = "-log(p-value)",
    y = "",
    title = "Top 25 Ingenuity Canonical Pathways"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 13),
    panel.grid.major.y = element_blank()
  )

fig_3b


# Figure 3c

ipa_ca1_pseudo$`z-score` <- as.numeric(as.character(ipa_ca1_pseudo$`z-score`))
ipa_ca1_pseudo$`-log(p-value)` <- as.numeric(as.character(ipa_ca1_pseudo$`-log(p-value)`))
ipa_ca1_pseudo$Ratio <- as.numeric(as.character(ipa_ca1_pseudo$Ratio))

top_df <- ipa_ca1_pseudo %>%
  arrange(desc(`-log(p-value)`)) %>%
  slice_head(n = 25)

top_df$`Ingenuity Canonical Pathways` <- factor(
  top_df$`Ingenuity Canonical Pathways`,
  levels = top_df$`Ingenuity Canonical Pathways`[order(top_df$`-log(p-value)`)]
)

blue_labels <- c(
  "Neurotransmitter uptake and metabolism In glial cells",
  "Neurotransmitter release cycle",
  "Glutamate Receptor Signaling",
  "Glutaminergic Receptor Signaling Pathway (Enhanced)",
  "Neurotrophin/TRK Signaling"
)

purple_labels <- c(
  "Degradation of the extracellular matrix",
  "Collagen biosynthesis and modifying enzymes",
  "Collagen chain trimerization",
  "Assembly of collagen fibrils and other multimeric structures"
)

green_labels <- c(
  "Neuroinflammation Signaling Pathway",
  "Wound Healing Signaling Pathway",
  "IL-33 Signaling Pathway",
  "IL-4 Signaling",
  "IL-12 Signaling and Production in Macrophages",
  "CD40 Signaling"
)

red_labels <- c(
  "Myelination Signaling Pathway",
  "Axonal Guidance Signaling"
)

label_color_map <- c(
  setNames(rep("blue", length(blue_labels)), blue_labels),
  setNames(rep("purple", length(purple_labels)), purple_labels),
  setNames(rep("#2C792D", length(green_labels)), green_labels),
  setNames(rep("darkred", length(red_labels)), red_labels)
)

label_fun <- function(x) {
  cols <- ifelse(x %in% names(label_color_map), label_color_map[x], "black")
  paste0("<span style='color:", cols, "'>", x, "</span>")
}

fig_3c <- ggplot(
  top_df,
  aes(
    x = `-log(p-value)`,
    y = `Ingenuity Canonical Pathways`,
    size = Ratio,
    fill = `z-score`
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey85",
    name = "Z-score"
  ) +
  scale_size_continuous(name = "Ratio", range = c(3, 10)) +
  scale_y_discrete(labels = label_fun) +
  labs(
    x = "-log(p-value)",
    y = "",
    title = "Top 25 Ingenuity Canonical Pathways"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 13),
    panel.grid.major.y = element_blank()
  )

fig_3c

# Figure 3d

ipa_ca3_pseudo$`z-score` <- as.numeric(as.character(ipa_ca3_pseudo$`z-score`))
ipa_ca3_pseudo$`-log(p-value)` <- as.numeric(as.character(ipa_ca3_pseudo$`-log(p-value)`))
ipa_ca3_pseudo$Ratio <- as.numeric(as.character(ipa_ca3_pseudo$Ratio))

top_df <- ipa_ca3_pseudo %>%
  arrange(desc(`-log(p-value)`)) %>%
  slice_head(n = 25)

top_df$`Ingenuity Canonical Pathways` <- factor(
  top_df$`Ingenuity Canonical Pathways`,
  levels = top_df$`Ingenuity Canonical Pathways`[order(top_df$`-log(p-value)`)]
)

blue_labels <- c(
  "Neurotransmitter uptake and metabolism In glial cells",
  "Neurotransmitter release cycle",
  "Glutamate Receptor Signaling",
  "Glutaminergic Receptor Signaling Pathway (Enhanced)",
  "Neurotrophin/TRK Signaling"
)

purple_labels <- c(
  "Degradation of the extracellular matrix",
  "Collagen biosynthesis and modifying enzymes",
  "Collagen chain trimerization",
  "Assembly of collagen fibrils and other multimeric structures"
)

green_labels <- c(
  "Neuroinflammation Signaling Pathway",
  "Wound Healing Signaling Pathway",
  "IL-33 Signaling Pathway",
  "IL-4 Signaling",
  "IL-12 Signaling and Production in Macrophages",
  "CD40 Signaling"
)

red_labels <- c(
  "Myelination Signaling Pathway",
  "Axonal Guidance Signaling"
)

label_color_map <- c(
  setNames(rep("blue", length(blue_labels)), blue_labels),
  setNames(rep("purple", length(purple_labels)), purple_labels),
  setNames(rep("#2C792D", length(green_labels)), green_labels),
  setNames(rep("darkred", length(red_labels)), red_labels)
)

label_fun <- function(x) {
  cols <- ifelse(x %in% names(label_color_map), label_color_map[x], "black")
  paste0("<span style='color:", cols, "'>", x, "</span>")
}

fig_3d <- ggplot(
  top_df,
  aes(
    x = `-log(p-value)`,
    y = `Ingenuity Canonical Pathways`,
    size = Ratio,
    fill = `z-score`
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey85",
    name = "Z-score"
  ) +
  scale_size_continuous(name = "Ratio", range = c(3, 10)) +
  scale_y_discrete(labels = label_fun) +
  labs(
    x = "-log(p-value)",
    y = "",
    title = "Top 25 Ingenuity Canonical Pathways"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 13),
    panel.grid.major.y = element_blank()
  )

fig_3d

# Figure 3e

astro_filt <- ipa_astrocytes_pseudo %>%
  filter(`-log(p-value)` > 1.3)

dentate_filt <- ipa_dentate_pseudo %>%
  filter(`-log(p-value)` > 1.3)

common_pathways <- intersect(
  astro_filt$`Ingenuity Canonical Pathways`,
  dentate_filt$`Ingenuity Canonical Pathways`
)

astro_common <- astro_filt %>%
  filter(`Ingenuity Canonical Pathways` %in% common_pathways) %>%
  select(
    Pathway = `Ingenuity Canonical Pathways`,
    z_astro = `z-score`
  ) %>%
  mutate(z_astro = as.numeric(z_astro))

dentate_common <- dentate_filt %>%
  filter(`Ingenuity Canonical Pathways` %in% common_pathways) %>%
  select(
    Pathway = `Ingenuity Canonical Pathways`,
    z_dentate = `z-score`
  ) %>%
  mutate(z_dentate = as.numeric(z_dentate))

merged <- inner_join(astro_common, dentate_common, by = "Pathway") %>%
  mutate(avg_z = rowMeans(abs(select(., z_astro, z_dentate)), na.rm = TRUE)) %>%
  arrange(desc(avg_z)) %>%
  slice(1:20)

pathway_order <- merged %>%
  arrange(desc(avg_z)) %>%
  pull(Pathway)

heatmap_data <- merged %>%
  select(Pathway, z_astro, z_dentate, avg_z) %>%
  pivot_longer(
    cols = c(z_astro, z_dentate),
    names_to = "Cell_Type",
    values_to = "z_score"
  ) %>%
  mutate(
    Cell_Type = recode(
      Cell_Type,
      "z_astro" = "Astrocytes",
      "z_dentate" = "Dentate gyrus"
    ),
    Pathway = factor(Pathway, levels = pathway_order)
  )

lim <- max(abs(heatmap_data$z_score), na.rm = TRUE)

fig_3e <- ggplot(
  heatmap_data,
  aes(
    x = Cell_Type,
    y = Pathway,
    fill = z_score
  )
) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    na.value = "grey85",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = squish,
    name = "IPA z-score"
  ) +
  labs(
    title = "Top 20 Enriched Pathways shared by Astrocytes-Dentate Gyrus",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.text.x = element_text(face = "bold", size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(10, "mm"),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
    plot.margin = margin(6, 12, 6, 6)
  )

fig_3e

# Figure 3f

ipa_upstream_pseudo_comp <- read_csv(
  "output/IPA_results_/comparison/upstream_regulators.csv"
)

heatmap_data <- ipa_upstream_pseudo_comp[1:25, ]

heatmap_data <- heatmap_data %>%
  column_to_rownames("Upstream Regulators")

heatmap_data$`NP CT L6b Glut_pseudo` <- NULL

heatmap_data <- heatmap_data %>%
  rename(
    `Dentate Gyrus Neurons` = `Dendate Gyrus Neurons_pseudo`,
    `CA1 Glut Neurons` = `CA1-ProS Glut Neurons_pseudo`,
    `CA3 Glut Neurons` = `CA3 Glut Neurons_pseudo`,
    `Astrocytes` = `Astrocytes_pseudo`
  )

heatmap_data[] <- lapply(heatmap_data, as.numeric)

regulator_order <- rownames(heatmap_data)

plot_data <- heatmap_data %>%
  rownames_to_column("Regulator") %>%
  pivot_longer(
    cols = -Regulator,
    names_to = "Cell_Type",
    values_to = "z_score"
  ) %>%
  mutate(
    Regulator = factor(Regulator, levels = rev(regulator_order))
  )

fig_3f <- ggplot(
  plot_data,
  aes(
    x = Cell_Type,
    y = Regulator,
    fill = z_score
  )
) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    na.value = "grey85",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "IPA z-score"
  ) +
  labs(
    title = "Top 25 Upstream Regulators",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.text.x = element_text(angle = 60, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(10, "mm"),
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.4)
  )

fig_3f
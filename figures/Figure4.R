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

# Figure 4a


fig_4a_1 <- netVisual_heatmap(cellchat)

fig_4a_2 <- netVisual_heatmap(
  cellchat,
  measure = "weight"
)

fig_4a <- fig_4a_1 + fig_4a_2

fig_4a


# Figure 4b

subtype_total <- full_join(out_sum, in_sum, by = "subtype") %>%
  mutate(
    across(
      c(n_out, n_in, sum_abs_out, sum_abs_in),
      ~ replace_na(., 0)
    ),
    n_edges_total = n_out + n_in,
    sum_abs_FC_total = sum_abs_out + sum_abs_in
  )

subtype_top <- subtype_total %>%
  arrange(desc(sum_abs_FC_total), desc(n_edges_total)) %>%
  slice_head(n = 20)

neuronal <- c(
  "DG-PIR Immature Neurons", "Lamp5 GABA Neurons", "Dentate Gyrus Neurons",
  "Immature GABA Neurons", "CTX-CGE GABA", "VIP GABA Neurons",
  "CA2 Glut Neurons", "CA3 Glut Neurons", "SUB-ProS Glut Neurons",
  "Pvalb GABA Neurons", "NP CT L6b Glut", "CA1 Glut Neurons",
  "Sst GABA Neurons", "IT-ET Glut Neurons"
)

glial <- c(
  "Astrocytes", "Oligodendrocytes", "Microglia", "Oligo Progenitor Cells"
)

other <- c(
  "Choroid Plexus", "Vascular cells"
)

subtype_top <- subtype_top %>%
  mutate(
    class = case_when(
      subtype %in% neuronal ~ "Neuronal",
      subtype %in% glial ~ "Glial",
      subtype %in% other ~ "Other",
      TRUE ~ "Other"
    ),
    subtype = fct_reorder(subtype, sum_abs_FC_total)
  )

cols_class <- c(
  "Neuronal" = "#9ecae1",
  "Glial" = "#fcbba1",
  "Other" = "#d9d9d9"
)

fig_4b <- ggplot(
  subtype_top,
  aes(x = subtype, y = sum_abs_FC_total, fill = class)
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", sum_abs_FC_total)),
    hjust = -0.2,
    size = 3.2
  ) +
  coord_flip() +
  scale_fill_manual(values = cols_class, name = "Category") +
  expand_limits(y = max(subtype_top$sum_abs_FC_total) * 1.15) +
  labs(
    x = NULL,
    y = "Total |Δprob| (incoming + outgoing)",
    title = "Most altered subtypes (sum |Δprob|)"
  ) +
  theme_classic(base_size = 11)

fig_4b


# Figure 4c

path_summary <- paths_merged %>%
  mutate(
    is_changed = change_category != "unchanged",
    abs_FC = abs(FC)
  ) %>%
  group_by(pathway_name) %>%
  summarise(
    n_edges = n(),
    n_changed = sum(is_changed, na.rm = TRUE),
    changed_frac = n_changed / n_edges,
    sum_FC = sum(FC, na.rm = TRUE),
    sum_abs_FC = sum(abs_FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    direction_ratio = sum_FC / sum_abs_FC,
    direction_bias = case_when(
      direction_ratio > 0.4 ~ "Increased",
      direction_ratio < -0.4 ~ "Decreased",
      TRUE ~ "Mixed"
    ),
    direction_bias = factor(
      direction_bias,
      levels = c("Increased", "Mixed", "Decreased")
    )
  )

path_top20 <- path_summary %>%
  arrange(desc(sum_abs_FC), desc(changed_frac), desc(n_changed)) %>%
  slice_head(n = 20) %>%
  mutate(pathway_name = fct_reorder(pathway_name, sum_abs_FC))

fig_4c <- ggplot(
  path_top20,
  aes(
    x = pathway_name,
    y = sum_abs_FC,
    size = changed_frac,
    colour = direction_bias
  )
) +
  geom_point(alpha = 0.9) +
  coord_flip() +
  scale_y_continuous(
    trans = "log1p",
    name = "Cumulative signaling change (log scale)",
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  scale_size_continuous(name = "Fraction of altered L-R pairs") +
  scale_colour_manual(
    name = "Net pathway bias",
    values = c(
      "Increased" = "red",
      "Mixed" = "#F5E8DF",
      "Decreased" = "blue"
    )
  ) +
  labs(
    x = NULL,
    title = "Top altered signaling pathways"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )

fig_4c


# Figure 4d

glut_edges <- paths_merged %>%
  filter(pathway_name == "Glutamate")

glut_out <- glut_edges %>%
  group_by(subtype = source) %>%
  summarise(
    n_out = n(),
    sum_out = sum(FC, na.rm = TRUE),
    .groups = "drop"
  )

glut_in <- glut_edges %>%
  group_by(subtype = target) %>%
  summarise(
    n_in = n(),
    sum_in = sum(FC, na.rm = TRUE),
    .groups = "drop"
  )

glut_summary <- full_join(glut_out, glut_in, by = "subtype") %>%
  mutate(
    across(c(n_out, n_in, sum_out, sum_in), ~ replace_na(., 0)),
    n_edges = n_out + n_in,
    sum_FC = sum_out + sum_in
  ) %>%
  arrange(sum_FC) %>%
  mutate(subtype = factor(subtype, levels = subtype))

fig_4d <- ggplot(glut_summary, aes(x = sum_FC, y = subtype)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(size = n_edges), colour = "grey35", alpha = 0.9) +
  scale_size_continuous(
    name = "Number of L-R pairs",
    range = c(2, 9)
  ) +
  labs(
    x = "Net glutamatergic change (sum Δprob)",
    y = NULL,
    title = "Glutamate signaling alterations by subtype"
  ) +
  theme_classic()

fig_4d

# Figure 4e

class_df <- plot_df %>%
  group_by(receptor_class, partner, direction) %>%
  summarise(
    mean_delta = mean(delta_prob),
    sum_delta = sum(delta_prob),
    sum_abs = sum(abs_delta),
    n_pairs = n(),
    .groups = "drop"
  )

desired_order <- c(
  "Dentate Gyrus Neurons", "DG-PIR Immature Neurons",
  "CA1 Glut Neurons", "CA3 Glut Neurons",
  "HPF CR Glut Neurons", "SUB-ProS Glut Neurons",
  "IT-ET Glut Neurons", "NP CT L6b Glut",
  "Sst GABA Neurons", "Pvalb GABA Neurons",
  "Lamp5 GABA Neurons", "VIP GABA Neurons",
  "CTX-CGE GABA",
  "Oligo Progenitor Cells", "Oligodendrocytes",
  "Astrocytes", "Microglia",
  "Vascular cells", "Choroid Plexus"
)

desired_order_full <- desired_order[
  desired_order %in% unique(as.character(class_df$partner))
]

class_df <- class_df %>%
  filter(partner %in% desired_order_full) %>%
  mutate(
    partner = factor(partner, levels = desired_order_full),
    receptor_class = factor(
      receptor_class,
      levels = c("AMPA", "NMDA", "Kainate", "mGluR")
    )
  )

max_abs <- max(abs(class_df$sum_delta), na.rm = TRUE)

fig_4e <- ggplot(
  class_df,
  aes(x = partner, y = receptor_class, fill = sum_delta)
) +
  geom_tile(color = "grey90", linewidth = 0.25) +
  facet_grid(. ~ direction, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "Δprob (sum)"
  ) +
  labs(
    x = "Interacting cell type",
    y = "Receptor class",
    title = "Glutamatergic signaling changes involving CA2 neurons"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

fig_4e




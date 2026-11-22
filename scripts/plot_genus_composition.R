#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(scales)
  library(ggh4x)
})

root <- "C:/Users/nojus/OneDrive/Desktop/Bubble project"
fig_dir <- file.path(root, "taxonomy_0.1", "figures")
final_dir <- file.path(root, "final_results")
comm <- read_tsv(
  file.path(root, "taxonomy_0.1", "tables", "community_top_genera_relabund_long.tsv"),
  show_col_types = FALSE
)

season_order <- c("P", "V", "R")
habitat_order <- c(
  "Muddy sediment",
  "Sandy sediment",
  "Muddy surface",
  "Muddy bottom",
  "Sandy surface"
)

habitat_of <- function(replicate_group) {
  dplyr::case_when(
    replicate_group %in% c("PN", "VN", "RN") ~ habitat_order[1],
    replicate_group %in% c("PV", "VV", "RV") ~ habitat_order[2],
    replicate_group %in% c("PNS", "VNS", "RNS") ~ habitat_order[3],
    replicate_group %in% c("PNB", "VNB", "RNB") ~ habitat_order[4],
    replicate_group %in% c("PVF", "VVF", "RVF") ~ habitat_order[5],
    TRUE ~ NA_character_
  )
}

genus_levels <- comm %>%
  group_by(genus_plot) %>%
  summarise(m = mean(rel_abund), .groups = "drop") %>%
  arrange(genus_plot == "Other", desc(m)) %>%
  pull(genus_plot)
genus_levels <- c(setdiff(genus_levels, "Other"), "Other")

comm <- comm %>%
  mutate(
    habitat = factor(habitat_of(replicate_group), levels = habitat_order),
    season_label = factor(
      recode(as.character(season), P = "Spring", V = "Summer", R = "Autumn"),
      levels = c("Spring", "Summer", "Autumn")
    ),
    season = factor(season, levels = season_order),
    genus_plot = factor(genus_plot, levels = genus_levels)
  )

sample_levels <- comm %>%
  distinct(sample, habitat, season, season_label, replicate_group, rep_num) %>%
  arrange(habitat, season, rep_num, sample) %>%
  pull(sample)
comm$sample <- factor(comm$sample, levels = sample_levels)

comm_pal <- setNames(
  c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#E6AB02", "#7570B3", "#66A61E", "#E7298A", "#1B9E77",
    "grey75"
  ),
  genus_levels
)
comm_pal["Other"] <- "grey75"

season_fill <- c(Spring = "#43a047", Summer = "#fb8c00", Autumn = "#1e88e5")
habitat_fill <- c(
  "Muddy sediment" = "#d6e0ea",
  "Sandy sediment" = "#f5d6c4",
  "Muddy surface" = "#c5e4f7",
  "Muddy bottom" = "#b7c9dc",
  "Sandy surface" = "#f7dcc8"
)

inner_map <- comm %>%
  distinct(habitat, season_label) %>%
  arrange(habitat, season_label)
outer_fills <- unname(habitat_fill[habitat_order])
inner_fills <- unname(season_fill[as.character(inner_map$season_label)])
n_outer <- length(outer_fills)
n_inner <- length(inner_fills)

p_comm <- ggplot(comm, aes(x = sample, y = rel_abund, fill = genus_plot)) +
  geom_col(width = 0.9) +
  facet_nested(
    ~ habitat + season_label,
    scales = "free_x",
    space = "free_x",
    nest_line = element_line(colour = "grey35", linewidth = 0.3),
    strip = strip_nested(
      clip = "off",
      bleed = FALSE,
      background_x = elem_list_rect(
        fill = c(outer_fills, inner_fills),
        colour = "grey35"
      ),
      text_x = elem_list_text(
        colour = c(rep("grey15", n_outer), rep("white", n_inner)),
        face = rep("bold", n_outer + n_inner),
        size = c(rep(15, n_outer), rep(11.5, n_inner))
      )
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = comm_pal) +
  labs(
    title = "Relative abundance of the 15 most abundant genera",
    x = "Sample",
    y = "Relative abundance",
    fill = "Genus"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11, face = "italic"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(4, "pt"),
    plot.margin = margin(5.5, 4, 5.5, 5.5),
    strip.text.x = element_text(face = "bold")
  )

panel_w <- comm %>%
  distinct(habitat, season_label, sample) %>%
  count(habitat, season_label, name = "n") %>%
  arrange(habitat, season_label) %>%
  mutate(
    w = if_else(
      as.character(habitat) == "Muddy surface" & as.character(season_label) == "Spring",
      1.45,
      as.numeric(n)
    )
  ) %>%
  pull(w)
p_comm <- p_comm + force_panelsizes(cols = panel_w)

outfile <- file.path(fig_dir, "community_top_genera_relabund.png")
ggsave(outfile, p_comm, width = 22, height = 9.6, dpi = 300)
ggsave(sub("\\.png$", ".pdf", outfile), p_comm, width = 22, height = 9.6)
file.copy(outfile, file.path(final_dir, "community_top_genera_relabund.png"), overwrite = TRUE)
file.copy(sub("\\.png$", ".pdf", outfile), file.path(final_dir, "community_top_genera_relabund.pdf"), overwrite = TRUE)
cat("Wrote", outfile, "\n")

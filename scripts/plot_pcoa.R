#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggrepel)
})

root <- "C:/Users/nojus/OneDrive/Desktop/Bubble project"
out_dir <- file.path(root, "taxonomy_0.1")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
final_dir <- file.path(root, "final_results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

bracken_g_dir <- file.path(
  out_dir, "input", "data", "nextflow", "taxprofiler",
  "detaxed", "results", "bracken", "k2_core_nt_G"
)

season_order <- c("P", "V", "R")
replicate_group_order <- c(
  "PN", "PNB", "PNS", "PV", "PVF",
  "VN", "VNB", "VNS", "VV", "VVF",
  "RN", "RNB", "RNS", "RV", "RVF"
)

dropped_samples <- c("PNS1", "PNS3")

paper_name <- function(samples) {
  vapply(samples, function(sample) {
    season <- c(P = "Sp", V = "Su", R = "Au")[[substr(sample, 1, 1)]]
    site <- c(N = "Mu", V = "Sa")[[substr(sample, 2, 2)]]
    rest <- substring(sample, 3)
    if (nchar(rest) > 0 && substr(rest, 1, 1) %in% c("S", "B", "F")) {
      depth <- if (substr(rest, 1, 1) %in% c("S", "F")) "S" else "B"
      paste0(season, site, depth, substring(rest, 2))
    } else {
      paste0(season, site, rest)
    }
  }, character(1), USE.NAMES = FALSE)
}

copy_to_final <- function(src, dest_name) {
  file.copy(src, file.path(final_dir, dest_name), overwrite = TRUE)
  paper_dir <- file.path(final_dir, "PAPER")
  dir.create(paper_dir, showWarnings = FALSE)
  file.copy(src, file.path(paper_dir, dest_name), overwrite = TRUE)
}

sample_meta <- function(sample_ids) {
  tibble(sample = sample_ids) %>%
    mutate(
      replicate_group = str_replace(sample, "[0-9]+$", ""),
      rep_num = suppressWarnings(as.integer(str_extract(sample, "[0-9]+$"))),
      season = factor(substr(sample, 1, 1), levels = season_order),
      replicate_group = factor(replicate_group, levels = replicate_group_order)
    )
}

order_samples <- function(sample_ids) {
  sample_meta(sample_ids) %>%
    arrange(season, replicate_group, rep_num, sample) %>%
    pull(sample)
}

files <- list.files(bracken_g_dir, pattern = "\\.bracken\\.tsv$", full.names = TRUE)
if (length(files) != 45) {
  stop("Expected 45 Bracken G TSVs, found ", length(files))
}

read_one_g <- function(path) {
  sid <- str_remove(basename(path), "_.*")
  df <- read_tsv(path, show_col_types = FALSE)
  df %>%
    filter(taxonomy_lvl == "G") %>%
    transmute(
      sample = sid,
      taxonomy_id = as.character(taxonomy_id),
      taxon_name = str_trim(name),
      reads = as.numeric(new_est_reads)
    )
}

long_counts <- bind_rows(lapply(files, read_one_g)) %>%
  mutate(reads = replace_na(reads, 0))

wide_counts <- long_counts %>%
  group_by(taxonomy_id, taxon_name, sample) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0)

sample_cols <- setdiff(colnames(wide_counts), c("taxonomy_id", "taxon_name"))
sample_cols <- order_samples(sample_cols)
wide_counts <- wide_counts[, c("taxonomy_id", "taxon_name", sample_cols)]

write_tsv(wide_counts, file.path(tab_dir, "genus_abundance_with_names_all_samples.tsv"))

keep_cols <- setdiff(sample_cols, dropped_samples)
wide_counts <- wide_counts[, c("taxonomy_id", "taxon_name", keep_cols)]
sample_cols <- keep_cols
long_counts <- long_counts %>% filter(!sample %in% dropped_samples)
mat <- t(as.matrix(wide_counts[, sample_cols, drop = FALSE]))
storage.mode(mat) <- "integer"
abund_mat <- mat[, colSums(mat) > 0, drop = FALSE]
sample_totals <- rowSums(abund_mat)
if (any(sample_totals == 0)) {
  stop("At least one sample has zero classified genus abundance.")
}
rel_abund_mat <- sweep(abund_mat, 1, sample_totals, "/")
bray <- vegdist(rel_abund_mat, method = "bray")
bray_matrix <- as.matrix(bray)

write.table(bray_matrix, file.path(tab_dir, "bray_curtis_distance_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(rel_abund_mat, file.path(tab_dir, "genus_relative_abundance_classified_only.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

pcoa <- cmdscale(bray, eig = TRUE, k = 2)
coords <- as.data.frame(pcoa$points)
colnames(coords) <- c("PCoA1", "PCoA2")
coords$sample <- rownames(coords)
coords <- coords[, c("sample", "PCoA1", "PCoA2")]
coords$replicate_group <- sub("[0-9]+$", "", coords$sample)
coords$season <- substr(coords$sample, 1, 1)

coords$coassembly_group <- ifelse(
  coords$replicate_group %in% c("PN", "RN", "VN"), "G1_sediment_Nida",
  ifelse(
    coords$replicate_group %in% c("PV", "RV", "VV"), "G2_sediment_Vidmares",
    ifelse(
      coords$replicate_group %in% c("PNB", "RNB", "RNS", "RVF") | coords$sample == "PNS2",
      "G3_water_Nida_autumn_Vidmares",
      ifelse(
        coords$replicate_group %in% c("PVF"), "G4_spring_surface",
        ifelse(
          coords$replicate_group %in% c("VNS", "VNB", "VVF"), "G5_summer_Vidmares_water",
          "Unassigned"
        )
      )
    )
  )
)

positive_eig <- pcoa$eig[pcoa$eig > 0]
var_exp <- 100 * positive_eig / sum(positive_eig)
coords$PCoA2 <- -coords$PCoA2
coords$habitat <- case_when(
  coords$replicate_group %in% c("PN", "VN", "RN") ~ "Muddy sediment",
  coords$replicate_group %in% c("PV", "VV", "RV") ~ "Sandy sediment",
  coords$replicate_group %in% c("PNS", "VNS", "RNS") ~ "Muddy surface",
  coords$replicate_group %in% c("PNB", "VNB", "RNB") ~ "Muddy bottom",
  coords$replicate_group %in% c("PVF", "VVF", "RVF") ~ "Sandy surface",
  TRUE ~ NA_character_
)
coords$habitat <- factor(coords$habitat, levels = c(
  "Sandy surface",
  "Muddy surface",
  "Muddy bottom",
  "Sandy sediment",
  "Muddy sediment"
))
coords$season_label <- factor(
  recode(coords$season, P = "Spring", V = "Summer", R = "Autumn"),
  levels = c("Spring", "Summer", "Autumn")
)
coords$ellipse_group <- case_when(
  coords$habitat == "Muddy sediment" ~ "Muddy sediment",
  coords$habitat == "Sandy sediment" ~ "Sandy sediment",
  coords$replicate_group %in% c("VNS", "VNB", "VVF") ~ "Summer water",
  coords$replicate_group == "PVF" ~ "Sandy spring surface",
  TRUE ~ "Muddy and autumn water"
)
coords$ellipse_group <- factor(coords$ellipse_group, levels = c(
  "Sandy spring surface",
  "Summer water",
  "Muddy and autumn water",
  "Sandy sediment",
  "Muddy sediment"
))
coords$paper_id <- paper_name(coords$sample)

write_tsv(coords, file.path(tab_dir, "pcoa_coordinates_with_groups.tsv"))
write_tsv(
  data.frame(axis = c("PCoA1", "PCoA2"), variance_explained_percent = c(var_exp[1], var_exp[2])),
  file.path(tab_dir, "pcoa_variance_explained.tsv")
)

group_sizes <- table(coords$coassembly_group)
valid_groups <- names(group_sizes[group_sizes >= 2])
coords_perm <- coords[coords$coassembly_group %in% valid_groups, , drop = FALSE]
dist_perm <- as.dist(as.matrix(bray)[coords_perm$sample, coords_perm$sample])
ad_global <- vegan::adonis2(dist_perm ~ coassembly_group, data = coords_perm, permutations = 9999)
ad_global_df <- as.data.frame(ad_global)
ad_global_df$term <- rownames(ad_global_df)
write_tsv(ad_global_df, file.path(tab_dir, "permanova_global.tsv"))

grp_levels <- sort(unique(coords_perm$coassembly_group))
pairs <- combn(grp_levels, 2, simplify = FALSE)
pair_results <- lapply(pairs, function(pr) {
  sub <- coords_perm[coords_perm$coassembly_group %in% pr, , drop = FALSE]
  dsub <- as.dist(as.matrix(bray)[sub$sample, sub$sample])
  a <- vegan::adonis2(dsub ~ coassembly_group, data = sub, permutations = 9999)
  a_df <- as.data.frame(a)
  data.frame(
    group_1 = pr[1], group_2 = pr[2],
    F_model = a_df$F[1], R2 = a_df$R2[1], p_value = a_df$`Pr(>F)`[1]
  )
})
pair_df <- bind_rows(pair_results)
pair_df$p_adj_fdr <- p.adjust(pair_df$p_value, method = "BH")
write_tsv(pair_df, file.path(tab_dir, "permanova_pairwise.tsv"))

bd <- vegan::betadisper(dist_perm, group = as.factor(coords_perm$coassembly_group))
write_tsv(as.data.frame(anova(bd)), file.path(tab_dir, "betadisper_anova.tsv"))

season_colors <- c(
  Spring = "#43a047",
  Summer = "#d73027",
  Autumn = "#1e88e5"
)
habitat_shapes <- c(
  "Sandy surface" = 21,
  "Muddy surface" = 24,
  "Muddy bottom" = 22,
  "Sandy sediment" = 23,
  "Muddy sediment" = 25
)

x_range <- range(coords$PCoA1, na.rm = TRUE)
y_range <- range(coords$PCoA2, na.rm = TRUE)
x_pad <- 0.18 * diff(x_range)
y_pad <- 0.18 * diff(y_range)
ell_ok <- coords %>%
  count(ellipse_group) %>%
  filter(n >= 4) %>%
  pull(ellipse_group)
coords_ell <- filter(coords, ellipse_group %in% ell_ok)

p_pcoa <- ggplot(coords, aes(x = PCoA1, y = PCoA2, fill = season_label, shape = habitat, label = paper_id)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  stat_ellipse(
    data = coords_ell,
    aes(x = PCoA1, y = PCoA2, group = ellipse_group),
    type = "t",
    level = 0.95,
    linetype = 2,
    linewidth = 0.55,
    colour = "grey35",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(size = 3.4, stroke = 0.45, colour = "grey20") +
  geom_text_repel(
    size = 3.2,
    colour = "grey20",
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.45,
    point.padding = 0.28,
    force = 2.2,
    force_pull = 0.6,
    seed = 1,
    bg.color = "white",
    bg.r = 0.12,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = season_colors, name = NULL) +
  scale_shape_manual(values = habitat_shapes, name = NULL) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, colour = "grey20", size = 3.4)),
    shape = guide_legend(override.aes = list(fill = "grey70", colour = "grey20", size = 3.4))
  ) +
  coord_cartesian(
    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    clip = "off"
  ) +
  labs(
    x = sprintf("PCoA1 (%.2f%%)", var_exp[1]),
    y = sprintf("PCoA2 (%.2f%%)", var_exp[2])
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.015, 0.02),
    legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(2, "pt"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.text = element_text(size = 10),
    plot.margin = margin(12, 12, 12, 12)
  )

pcoa_png <- file.path(fig_dir, "pcoa_bray_curtis.png")
pcoa_pdf <- file.path(fig_dir, "pcoa_bray_curtis.pdf")
ggsave(pcoa_png, p_pcoa, width = 11.4, height = 9.0, dpi = 300)
ggsave(pcoa_pdf, p_pcoa, width = 11.4, height = 9.0)
copy_to_final(pcoa_png, "pcoa_bray_curtis.png")
copy_to_final(pcoa_pdf, "pcoa_bray_curtis.pdf")

classified <- long_counts %>%
  group_by(sample) %>%
  summarise(classified_reads = sum(reads), .groups = "drop")

top_n <- 15
genus_rel <- long_counts %>%
  left_join(classified, by = "sample") %>%
  mutate(rel_abund = if_else(classified_reads > 0, reads / classified_reads, 0))
top_genera <- genus_rel %>%
  group_by(taxon_name) %>%
  summarise(mean_rel = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n = top_n) %>%
  pull(taxon_name)
comm <- genus_rel %>%
  mutate(genus_plot = if_else(taxon_name %in% top_genera, taxon_name, "Other")) %>%
  group_by(sample, genus_plot) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>%
  left_join(sample_meta(unique(genus_rel$sample)), by = "sample") %>%
  mutate(
    sample = factor(sample, levels = order_samples(unique(sample))),
    genus_plot = factor(genus_plot, levels = c(top_genera, "Other")),
    season_label = factor(
      recode(as.character(season), P = "Spring", V = "Summer", R = "Autumn"),
      levels = c("Spring", "Summer", "Autumn")
    )
  )
write_tsv(comm, file.path(tab_dir, "community_top_genera_relabund_long.tsv"))

cat("Excluded from PCoA/taxonomy: ", paste(dropped_samples, collapse = ", "), "\n", sep = "")
cat("n samples in PCoA: ", nrow(coords), "\n", sep = "")
cat("PCoA variance: PCoA1 ", round(var_exp[1], 2), "%, PCoA2 ", round(var_exp[2], 2), "%\n", sep = "")
cat("Wrote figures to ", fig_dir, " and ", final_dir, "\n", sep = "")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggh4x)
  library(cowplot)
  library(grid)
  library(gtable)
  library(scales)
})

root <- "C:/Users/nojus/OneDrive/Desktop/Bubble project"
repo <- "C:/Users/nojus/OneDrive/Desktop/nida-vidmares-metagenomics"
fig_dir <- file.path(root, "taxonomy_0.1", "figures")
tab_dir <- file.path(root, "taxonomy_0.1", "tables")
final_dir <- file.path(root, "final_results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

bracken_g_dir <- file.path(
  root, "taxonomy_0.1", "input", "data", "nextflow", "taxprofiler",
  "detaxed", "results", "bracken", "k2_core_nt_G"
)
dropped <- c("PNS1", "PNS3")

curated <- read_tsv(file.path(repo, "data", "methane_cycling_genera.tsv"), show_col_types = FALSE)
canon <- function(x) {
  str_squish(str_remove(x, regex("^Candidatus\\s+", ignore_case = TRUE)))
}
curated <- curated %>%
  mutate(
    genus_canon = recode(canon(genus), Methanosaeta = "Methanothrix"),
    guild_label = recode(
      guild,
      methanogen = "Methanogens",
      methanotroph = "Methanotrophs"
    )
  )
guild_of <- curated %>%
  distinct(genus_canon, guild_label)

files <- list.files(bracken_g_dir, pattern = "\\.bracken\\.tsv$", full.names = TRUE)
read_one_g <- function(path) {
  sid <- str_remove(basename(path), "_.*")
  read_tsv(path, show_col_types = FALSE) %>%
    filter(taxonomy_lvl == "G") %>%
    transmute(
      sample = sid,
      taxon_name = str_trim(name),
      genus_canon = recode(canon(str_trim(name)), Methanosaeta = "Methanothrix"),
      reads = as.numeric(new_est_reads)
    )
}

long_counts <- bind_rows(lapply(files, read_one_g)) %>%
  mutate(reads = replace_na(reads, 0)) %>%
  filter(!sample %in% dropped)

classified <- long_counts %>%
  group_by(sample) %>%
  summarise(classified_reads = sum(reads), .groups = "drop")

df <- long_counts %>%
  inner_join(guild_of, by = "genus_canon") %>%
  group_by(sample, genus_canon, guild_label) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  left_join(classified, by = "sample") %>%
  mutate(rel_abund = if_else(classified_reads > 0, reads / classified_reads, 0))

keep_genera <- df %>%
  group_by(genus_canon) %>%
  summarise(tot = sum(reads), .groups = "drop") %>%
  filter(tot > 0) %>%
  pull(genus_canon)

habitat_order <- c(
  "Muddy sediment",
  "Sandy sediment",
  "Muddy surface",
  "Muddy bottom",
  "Sandy surface"
)
habitat_of <- function(replicate_group) {
  case_when(
    replicate_group %in% c("PN", "VN", "RN") ~ habitat_order[1],
    replicate_group %in% c("PV", "VV", "RV") ~ habitat_order[2],
    replicate_group %in% c("PNS", "VNS", "RNS") ~ habitat_order[3],
    replicate_group %in% c("PNB", "VNB", "RNB") ~ habitat_order[4],
    replicate_group %in% c("PVF", "VVF", "RVF") ~ habitat_order[5],
    TRUE ~ NA_character_
  )
}

df <- df %>%
  mutate(
    taxon_plot = if_else(genus_canon %in% keep_genera, genus_canon, "Other"),
    replicate_group = sub("[0-9]+$", "", sample),
    rep_num = as.integer(sub(".*?([0-9]+)$", "\\1", sample)),
    season = substr(sample, 1, 1),
    habitat = factor(habitat_of(replicate_group), levels = habitat_order),
    season_label = factor(
      recode(season, P = "Spring", V = "Summer", R = "Autumn"),
      levels = c("Spring", "Summer", "Autumn")
    ),
    guild_label = factor(guild_label, levels = c("Methanogens", "Methanotrophs"))
  ) %>%
  filter(taxon_plot != "Other") %>%
  group_by(sample, taxon_plot, guild_label, habitat, season_label, season, replicate_group, rep_num) %>%
  summarise(rel_abund = sum(rel_abund), reads = sum(reads), .groups = "drop")

named_order <- function(guild) {
  df %>%
    filter(guild_label == guild, taxon_plot != "Other") %>%
    group_by(taxon_plot) %>%
    summarise(tot = sum(reads), .groups = "drop") %>%
    arrange(desc(tot)) %>%
    pull(taxon_plot) %>%
    as.character()
}
producer_gen <- named_order("Methanogens")
oxidiser_gen <- named_order("Methanotrophs")
genus_levels <- c(producer_gen, oxidiser_gen)

df <- df %>%
  group_by(sample, taxon_plot, habitat, season_label, season, replicate_group, rep_num) %>%
  summarise(rel_abund = sum(rel_abund), reads = sum(reads), .groups = "drop")

sample_meta <- df %>%
  distinct(sample, habitat, season_label, season, replicate_group, rep_num)
plot_df <- sample_meta %>%
  crossing(taxon_plot = unique(df$taxon_plot)) %>%
  left_join(
    df %>% select(sample, taxon_plot, rel_abund),
    by = c("sample", "taxon_plot")
  ) %>%
  mutate(
    rel_abund = replace_na(rel_abund, 0),
    pct = 100 * rel_abund
  )

sample_levels <- plot_df %>%
  distinct(sample, habitat, season, rep_num) %>%
  arrange(habitat, season, rep_num, sample) %>%
  pull(sample)

plot_df <- plot_df %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    taxon_plot = factor(taxon_plot, levels = genus_levels)
  )

reds <- c(
  "#bc3c29", "#c9b25a", "#d0763a", "#c47a9e",
  "#9c6b4a", "#8c510a", "#d4a8a6", "#d4b88a"
)
blues <- c(
  "#3c5488", "#3a7d5c", "#6b548c", "#7ba8c4", "#8a9a4a",
  "#332288", "#4a8b82", "#4a7c8c", "#2B3D26", "#5a9ba8",
  "#6a9a70", "#9e9ac8", "#01665e", "#a3c48a", "#5a4a8c",
  "#a6cee3", "#7aada8"
)
prod_pal <- setNames(reds[seq_along(producer_gen)], producer_gen)
ox_pal <- setNames(blues[seq_along(oxidiser_gen)], oxidiser_gen)
fill_pal <- c(prod_pal, ox_pal)

season_fill <- c(Spring = "#43a047", Summer = "#fb8c00", Autumn = "#1e88e5")
habitat_fill <- c(
  "Muddy sediment" = "#d6e0ea",
  "Sandy sediment" = "#f5d6c4",
  "Muddy surface" = "#c5e4f7",
  "Muddy bottom" = "#b7c9dc",
  "Sandy surface" = "#f7dcc8"
)
write_tsv(plot_df, file.path(tab_dir, "methane_cycling_genus_relabund_long.tsv"))

null_to_zero <- function(u) {
  u[grepl("null", as.character(u))] <- unit(0, "pt")
  u
}

take_legend <- function(plt) {
  leg <- get_plot_component(plt + theme(legend.position = "right"), "guide-box-right")
  if (is.null(leg) || inherits(leg, "zeroGrob")) {
    return(leg)
  }
  leg[["widths"]] <- null_to_zero(leg[["widths"]])
  leg[["heights"]] <- null_to_zero(leg[["heights"]])
  leg
}

leg_plot <- function(items, pal, title, title_col) {
  d <- data.frame(x = 1, y = 1, g = factor(items, levels = items))
  ggplot(d, aes(x, y, fill = g)) +
    geom_col() +
    scale_fill_manual(name = title, values = pal[items], breaks = items) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 12.5, face = "italic"),
      legend.title = element_text(
        size = 14,
        face = "bold",
        colour = title_col,
        margin = margin(b = 2, t = 0)
      ),
      legend.key.size = unit(0.48, "cm"),
      legend.key.spacing.y = unit(1, "pt"),
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0, "pt"),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

leg_mg <- take_legend(leg_plot(producer_gen, prod_pal, "Methanogens", "#7f0000"))
leg_mt <- take_legend(leg_plot(oxidiser_gen, ox_pal, "Methanotrophs", "#0d47a1"))
gap <- gtable(widths = unit(1, "mm"), heights = unit(4.5, "mm"))
legs <- list(leg_mg, gap, leg_mt)
leg_h <- do.call(unit.c, lapply(legs, function(g) sum(g[["heights"]])))
leg_w <- do.call(unit.pmax, lapply(legs, function(g) sum(g[["widths"]])))
leg_col <- gtable(widths = leg_w, heights = leg_h)
for (i in seq_along(legs)) {
  leg_col <- gtable_add_grob(leg_col, legs[[i]], t = i, l = 1, clip = "off")
}

save_split <- function(dat, title, outfile, width, y_acc) {
  habs <- levels(droplevels(dat$habitat))
  inner_map <- dat %>%
    distinct(habitat, season_label) %>%
    arrange(habitat, season_label)
  outer_fills <- unname(habitat_fill[habs])
  inner_fills <- unname(season_fill[as.character(inner_map$season_label)])
  n_outer <- length(outer_fills)
  n_inner <- length(inner_fills)

  p <- ggplot(dat, aes(x = sample, y = rel_abund, fill = taxon_plot)) +
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
          size = c(rep(15, n_outer), rep(14.5, n_inner))
        )
      )
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = y_acc),
      expand = expansion(mult = c(0, 0.06))
    ) +
    scale_fill_manual(values = fill_pal, breaks = genus_levels, guide = "none") +
    labs(title = title, x = "Sample", y = "Relative abundance") +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(face = "bold")
    )

  panel_w <- dat %>%
    distinct(habitat, season_label, sample) %>%
    count(habitat, season_label, name = "n") %>%
    arrange(habitat, season_label) %>%
    mutate(
      w = if_else(
        as.character(habitat) == "Muddy surface" & as.character(season_label) == "Spring",
        1.7,
        as.numeric(n)
      )
    ) %>%
    pull(w)
  p <- p + force_panelsizes(cols = panel_w)

  p_grob <- ggplotGrob(p)
  out <- gtable(
    widths = unit.c(unit(1, "null"), unit(1.5, "mm"), leg_w, unit(2, "mm")),
    heights = unit(1, "null")
  )
  out <- gtable_add_grob(out, p_grob, t = 1, l = 1)
  out <- gtable_add_grob(
    out,
    grobTree(
      leg_col,
      vp = viewport(
        x = unit(0, "npc"),
        y = unit(0.5, "npc"),
        just = c("left", "center"),
        width = leg_w,
        height = sum(leg_h)
      )
    ),
    t = 1,
    l = 3,
    clip = "off"
  )

  png <- file.path(fig_dir, outfile)
  pdf <- sub("\\.png$", ".pdf", png)
  ggsave(png, out, width = width, height = 7.8, dpi = 300, bg = "white")
  ggsave(pdf, out, width = width, height = 7.8, bg = "white")
  file.copy(png, file.path(final_dir, outfile), overwrite = TRUE)
  file.copy(pdf, file.path(final_dir, basename(pdf)), overwrite = TRUE)
  cat("Wrote", png, "\n")
}

sed_habs <- habitat_order[1:2]
wat_habs <- habitat_order[3:5]
save_split(
  plot_df %>% filter(habitat %in% sed_habs) %>% droplevels(),
  "Methanogen and methanotroph genera (sediment)",
  "methane_cycling_genus_relabund_sediment.png",
  16,
  1
)
save_split(
  plot_df %>% filter(habitat %in% wat_habs) %>% droplevels(),
  "Methanogen and methanotroph genera (water)",
  "methane_cycling_genus_relabund_water.png",
  18,
  0.1
)
cat("Named genera:", paste(setdiff(genus_levels, "Other"), collapse = ", "), "\n")

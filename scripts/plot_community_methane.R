#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(ggh4x)
  library(scales)
  library(cowplot)
  library(grid)
  library(gtable)
})

root <- "C:/Users/nojus/OneDrive/Desktop/Bubble project"
repo <- "C:/Users/nojus/OneDrive/Desktop/nida-vidmares-metagenomics"
fig_dir <- file.path(root, "taxonomy_0.1", "figures")
final_dir <- file.path(root, "final_results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

habitat_order <- c(
  "Sandy surface",
  "Muddy surface",
  "Muddy bottom",
  "Sandy sediment",
  "Muddy sediment"
)
water_habs <- habitat_order[1:3]
sed_habs <- habitat_order[4:5]
season_fill <- c(
  " " = "#43a047",
  Spring = "#43a047",
  Summer = "#fb8c00",
  Autumn = "#1e88e5"
)
habitat_fill <- c(
  "Muddy sediment" = "#d6e0ea",
  "Sandy sediment" = "#f5d6c4",
  "Muddy surface" = "#c5e4f7",
  "Muddy bottom" = "#b7c9dc",
  "Sandy surface" = "#f7dcc8"
)

paper_name <- function(samples) {
  vapply(as.character(samples), function(sample) {
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

habitat_of <- function(replicate_group) {
  case_when(
    replicate_group %in% c("PN", "VN", "RN") ~ "Muddy sediment",
    replicate_group %in% c("PV", "VV", "RV") ~ "Sandy sediment",
    replicate_group %in% c("PNS", "VNS", "RNS") ~ "Muddy surface",
    replicate_group %in% c("PNB", "VNB", "RNB") ~ "Muddy bottom",
    replicate_group %in% c("PVF", "VVF", "RVF") ~ "Sandy surface",
    TRUE ~ NA_character_
  )
}

add_strips <- function(dat) {
  dat %>%
    mutate(
      habitat = factor(habitat, levels = habitat_order),
      season_label = factor(season_label, levels = c("Spring", "Summer", "Autumn")),
      season_strip = if_else(
        habitat == "Muddy surface" & season_label == "Spring",
        " ",
        as.character(season_label)
      ),
      season_strip = factor(season_strip, levels = c(" ", "Spring", "Summer", "Autumn"))
    )
}

hplot <- function(dat, fill_col, pal, breaks, title, y_acc, ncol_leg) {
  dat <- add_strips(dat)
  sample_levels <- dat %>%
    distinct(sample, habitat, season_label, replicate_group, rep_num) %>%
    arrange(habitat, season_label, rep_num, sample) %>%
    pull(sample)
  dat$sample <- factor(paper_name(dat$sample), levels = paper_name(sample_levels))
  dat[[fill_col]] <- factor(dat[[fill_col]], levels = breaks)

  habs <- levels(droplevels(dat$habitat))
  hab_fills <- unname(habitat_fill[habs])
  if (identical(title, "Sediment")) {
    new_habs <- sub(" sediment$", "", habs)
    dat$habitat <- factor(
      sub(" sediment$", "", as.character(dat$habitat)),
      levels = new_habs
    )
    habs <- new_habs
  }
  inner_map <- dat %>%
    distinct(habitat, season_strip) %>%
    arrange(habitat, season_strip)
  season_fills <- unname(season_fill[as.character(inner_map$season_strip)])
  n_hab <- length(hab_fills)
  n_season <- length(season_fills)
  use_comp <- nzchar(title)
  if (use_comp) {
    dat$compartment <- factor(title, levels = title)
    comp_fill <- if (title == "Water") "#b7d4ea" else "#e2d0c2"
    strip_fills <- c(comp_fill, hab_fills, season_fills)
    strip_cols <- c("grey15", rep("grey15", n_hab), rep("white", n_season))
    strip_sizes <- c(12.5, rep(11.5, n_hab), rep(11, n_season))
    facet_form <- ~ compartment + habitat + season_strip
  } else {
    strip_fills <- c(hab_fills, season_fills)
    strip_cols <- c(rep("grey15", n_hab), rep("white", n_season))
    strip_sizes <- c(rep(12, n_hab), rep(11, n_season))
    facet_form <- ~ habitat + season_strip
  }

  p <- ggplot(dat, aes(x = sample, y = rel_abund, fill = .data[[fill_col]])) +
    geom_col(width = 0.9) +
    facet_nested(
      facet_form,
      scales = "free_x",
      space = "free_x",
      nest_line = element_line(colour = "grey35", linewidth = 0.3),
      strip = strip_nested(
        clip = "off",
        bleed = FALSE,
        background_x = elem_list_rect(
          fill = strip_fills,
          colour = "grey35"
        ),
        text_x = elem_list_text(
          colour = strip_cols,
          face = rep("bold", length(strip_fills)),
          size = strip_sizes
        )
      )
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = y_acc),
      expand = expansion(mult = c(0, 0.06))
    ) +
    scale_fill_manual(
      values = setNames(unname(pal[as.character(breaks)]), as.character(breaks)),
      breaks = as.character(breaks),
      drop = FALSE
    ) +
    labs(x = NULL, y = "Relative abundance", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11.5),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.placement = "outside",
      plot.margin = margin(10, 4, 4, 4)
    )

  panel_w <- dat %>%
    distinct(habitat, season_strip, sample) %>%
    count(habitat, season_strip, name = "n") %>%
    arrange(habitat, season_strip) %>%
    mutate(w = as.numeric(n)) %>%
    pull(w)
  p + force_panelsizes(cols = panel_w)
}

vplot <- function(dat, fill_col, pal, breaks, title, y_acc, ncol_leg) {
  dat <- add_strips(dat)
  sample_levels <- dat %>%
    distinct(sample, habitat, season_label, replicate_group, rep_num) %>%
    arrange(habitat, season_label, rep_num, sample) %>%
    pull(sample)
  dat$sample <- factor(paper_name(dat$sample), levels = rev(paper_name(sample_levels)))
  dat[[fill_col]] <- factor(dat[[fill_col]], levels = breaks)

  habs <- levels(droplevels(dat$habitat))
  inner_map <- dat %>%
    distinct(habitat, season_strip) %>%
    arrange(habitat, season_strip)
  outer_fills <- unname(habitat_fill[habs])
  inner_fills <- unname(season_fill[as.character(inner_map$season_strip)])
  n_outer <- length(outer_fills)
  n_inner <- length(inner_fills)

  p <- ggplot(dat, aes(x = rel_abund, y = sample, fill = .data[[fill_col]])) +
    geom_col(width = 0.82) +
    facet_nested(
      habitat + season_strip ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y",
      nest_line = element_line(colour = "grey35", linewidth = 0.3),
      strip = strip_nested(
        clip = "off",
        bleed = FALSE,
        background_y = elem_list_rect(
          fill = c(outer_fills, inner_fills),
          colour = "grey35"
        ),
        text_y = elem_list_text(
          colour = c(rep("grey15", n_outer), rep("white", n_inner)),
          face = rep("bold", n_outer + n_inner),
          size = c(rep(11, n_outer), rep(8.5, n_inner))
        )
      )
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = y_acc),
      expand = expansion(mult = c(0, 0.04))
    ) +
    scale_fill_manual(
      values = setNames(unname(pal[as.character(breaks)]), as.character(breaks)),
      breaks = as.character(breaks),
      drop = FALSE
    ) +
    labs(x = "Relative abundance", y = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 13, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90),
      plot.margin = margin(2, 6, 2, 4)
    )

  panel_h <- dat %>%
    distinct(habitat, season_strip, sample) %>%
    count(habitat, season_strip, name = "n") %>%
    arrange(habitat, season_strip) %>%
    mutate(h = as.numeric(n)) %>%
    pull(h)
  p + force_panelsizes(rows = panel_h)
}

save_one <- function(p, outfile, width, height) {
  png <- file.path(fig_dir, outfile)
  pdf <- sub("\\.png$", ".pdf", png)
  png_out <- tryCatch({
    ggsave(png, p, width = width, height = height, dpi = 300, bg = "white")
    png
  }, error = function(e) {
    alt <- sub("\\.png$", "_export.png", png)
    ggsave(alt, p, width = width, height = height, dpi = 300, bg = "white")
    alt
  })
  pdf_out <- tryCatch({
    ggsave(pdf, p, width = width, height = height, bg = "white")
    pdf
  }, error = function(e) {
    alt <- sub("\\.pdf$", "_export.pdf", pdf)
    ggsave(alt, p, width = width, height = height, bg = "white")
    alt
  })
  try(file.copy(png_out, file.path(final_dir, basename(png_out)), overwrite = TRUE), silent = TRUE)
  try(file.copy(pdf_out, file.path(final_dir, basename(pdf_out)), overwrite = TRUE), silent = TRUE)
  paper_dir <- file.path(final_dir, "PAPER")
  dir.create(paper_dir, showWarnings = FALSE)
  try(file.copy(png_out, file.path(paper_dir, basename(png_out)), overwrite = TRUE), silent = TRUE)
  try(file.copy(pdf_out, file.path(paper_dir, basename(pdf_out)), overwrite = TRUE), silent = TRUE)
  cat("Wrote", png_out, "\nWrote", pdf_out, "\n")
}

null_to_zero <- function(u) {
  u[grepl("null", as.character(u))] <- unit(0, "pt")
  u
}

take_legend <- function(plt, pos = "right") {
  box <- if (pos == "bottom") "guide-box-bottom" else "guide-box-right"
  leg <- get_plot_component(plt + theme(legend.position = pos), box)
  if (is.null(leg) || inherits(leg, "zeroGrob")) {
    return(leg)
  }
  leg[["widths"]] <- null_to_zero(leg[["widths"]])
  leg[["heights"]] <- null_to_zero(leg[["heights"]])
  leg
}

leg_plot <- function(items, pal, title, title_col, pos = "right", ncol_leg = 1) {
  d <- data.frame(x = 1, y = 1, g = factor(items, levels = items))
  ggplot(d, aes(x, y, fill = g)) +
    geom_col() +
    scale_fill_manual(
      name = title,
      values = setNames(unname(pal[as.character(items)]), as.character(items)),
      breaks = as.character(items),
      limits = as.character(items),
      drop = FALSE
    ) +
    guides(fill = guide_legend(ncol = ncol_leg, byrow = FALSE, title.position = "top")) +
    theme_void() +
    theme(
      legend.position = pos,
      legend.text = element_text(size = 10, face = "italic"),
      legend.title = if (is.null(title)) {
        element_blank()
      } else {
        element_text(
          size = 11.5,
          face = "bold",
          colour = title_col,
          margin = margin(b = 2, t = 0)
        )
      },
      legend.key.size = unit(0.4, "cm"),
      legend.key.spacing.y = unit(1, "pt"),
      legend.key.spacing.x = unit(6, "pt"),
      legend.justification = "left",
      legend.margin = margin(0, 2, 0, 2),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

pad_leg <- function(items, pal, ncol = 4, face = "plain") {
  n <- length(items)
  nr <- ceiling(n / ncol)
  fills <- unname(pal[as.character(items)])
  gt <- gtable(
    widths = unit.c(unit(3, "mm"), rep(unit.c(unit(3.6, "mm"), unit(1, "null")), ncol)),
    heights = unit(c(1.2, rep(4.6, nr), 0.8), "mm")
  )
  for (i in seq_len(n)) {
    r <- ((i - 1) %% nr) + 2
    c <- floor((i - 1) / nr)
    l_sq <- 2 + c * 2
    gt <- gtable_add_grob(
      gt,
      rectGrob(
        width = unit(3.2, "mm"),
        height = unit(3.2, "mm"),
        gp = gpar(fill = fills[[i]], col = NA)
      ),
      t = r,
      l = l_sq
    )
    gt <- gtable_add_grob(
      gt,
      textGrob(
        items[[i]],
        x = unit(1.2, "mm"),
        just = c("left", "center"),
        gp = gpar(fontsize = 10, fontface = face)
      ),
      t = r,
      l = l_sq + 1,
      clip = "off"
    )
  }
  gt
}

frame_plot <- function(g, top = 8, right = 5, bottom = 8, left = 5) {
  grob <- if (inherits(g, "gtable")) g else ggplotGrob(g)
  out <- gtable(
    widths = unit.c(unit(left, "mm"), unit(1, "null"), unit(right, "mm")),
    heights = unit.c(unit(top, "mm"), unit(1, "null"), unit(bottom, "mm"))
  )
  gtable_add_grob(out, grob, t = 2, l = 2)
}

stack_left <- function(p, leg) {
  left <- gtable(
    widths = unit.c(unit(9, "mm"), unit(1, "null")),
    heights = unit.c(unit(1, "null"), sum(leg[["heights"]]))
  )
  left <- gtable_add_grob(left, ggplotGrob(p), t = 1, l = 2)
  left <- gtable_add_grob(left, leg, t = 2, l = 2)
  gtable_add_grob(
    left,
    textGrob(
      "A",
      x = unit(0.08, "npc"),
      y = unit(1, "npc") - unit(1.2, "mm"),
      just = c("left", "top"),
      gp = gpar(fontsize = 22, fontface = "bold")
    ),
    t = 1,
    l = 1,
    clip = "off"
  )
}

named_pal <- function(shared, rest, colors) {
  nms <- c(shared, setdiff(rest, shared))
  setNames(colors[seq_along(nms)], nms)
}

dropped <- c("PNS1", "PNS3")
rep_dir <- file.path(
  root, "taxonomy_0.1", "input", "data", "nextflow",
  "taxprofiler", "detaxed", "results", "bracken", "k2_core_nt_G"
)
rep_files <- list.files(
  rep_dir,
  pattern = "\\.bracken\\.kraken2\\.report_bracken\\.txt$",
  full.names = TRUE
)

parse_phyla <- function(path) {
  sample <- sub("_.*", "", basename(path))
  x <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c("pct", "clade", "assigned", "rank", "taxid", "name"),
    stringsAsFactors = FALSE
  )
  x$name <- trimws(x$name)
  domain <- NA_character_
  out <- list()
  for (i in seq_len(nrow(x))) {
    if (x$rank[i] == "R2") {
      domain <- x$name[i]
    }
    if (x$rank[i] == "P" && domain %in% c("Bacteria", "Archaea")) {
      out[[length(out) + 1]] <- data.frame(
        sample = sample,
        phylum = x$name[i],
        reads = as.numeric(x$clade[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

phy_raw <- bind_rows(lapply(rep_files, parse_phyla)) %>%
  filter(!sample %in% dropped)
long <- phy_raw %>%
  mutate(
    replicate_group = gsub("[0-9]+$", "", sample),
    rep_num = as.integer(sub("^.*[^0-9]", "", sample)),
    season = substr(sample, 1, 1),
    habitat = habitat_of(replicate_group),
    compartment = if_else(
      replicate_group %in% c("PN", "VN", "RN", "PV", "VV", "RV"),
      "sediment",
      "water"
    ),
    season_label = recode(season, P = "Spring", V = "Summer", R = "Autumn")
  ) %>%
  group_by(sample) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup()

mean_comp <- long %>%
  group_by(compartment, phylum) %>%
  summarise(m = mean(rel_abund), .groups = "drop")
keep_phy <- mean_comp %>%
  filter(m >= 0.01) %>%
  distinct(phylum) %>%
  pull(phylum)
mean_all <- long %>%
  filter(phylum %in% keep_phy) %>%
  group_by(phylum) %>%
  summarise(m = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(phylum)

comm <- long %>%
  mutate(phylum_plot = if_else(phylum %in% keep_phy, phylum, "Other")) %>%
  group_by(
    sample, phylum_plot, replicate_group, rep_num, season,
    habitat, compartment, season_label
  ) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

write_tsv(
  long,
  file.path(root, "taxonomy_0.1", "tables", "phylum_abundance_classified_only.tsv")
)

ch <- read_tsv(
  file.path(root, "taxonomy_0.1", "tables", "methane_cycling_genus_relabund_long.tsv"),
  show_col_types = FALSE
)
guilds <- read_tsv(file.path(repo, "data", "methane_cycling_genera.tsv"), show_col_types = FALSE) %>%
  mutate(
    taxon_plot = recode(
      gsub("^Candidatus\\s+", "", genus),
      Methanosaeta = "Methanothrix"
    ),
    guild_label = recode(guild, methanogen = "Methanogens", methanotroph = "Methanotrophs")
  ) %>%
  distinct(taxon_plot, guild_label)

ch <- ch %>%
  left_join(guilds, by = "taxon_plot") %>%
  mutate(
    habitat = habitat_of(replicate_group),
    guild_label = factor(guild_label, levels = c("Methanogens", "Methanotrophs"))
  )

named_order <- function(dat, guild) {
  dat %>%
    filter(guild_label == guild) %>%
    group_by(taxon_plot) %>%
    summarise(tot = sum(rel_abund), .groups = "drop") %>%
    filter(tot > 0) %>%
    arrange(desc(tot)) %>%
    pull(taxon_plot) %>%
    as.character()
}

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
mg_wat <- named_order(filter(ch, habitat %in% water_habs), "Methanogens")
mg_sed <- named_order(filter(ch, habitat %in% sed_habs), "Methanogens")
mt_wat <- named_order(filter(ch, habitat %in% water_habs), "Methanotrophs")
mt_sed <- named_order(filter(ch, habitat %in% sed_habs), "Methanotrophs")
mg_all <- named_order(ch, "Methanogens")
mt_all <- named_order(ch, "Methanotrophs")
mg_pal <- named_pal(mg_all[mg_all %in% intersect(mg_wat, mg_sed)], mg_all, reds)
mt_pal <- named_pal(mt_all[mt_all %in% intersect(mt_wat, mt_sed)], mt_all, blues)
ch_pal <- c(mg_pal, mt_pal)

phy_cols <- c(
  Pseudomonadota = "#3C5488",
  Actinomycetota = "#E64B35",
  Cyanobacteriota = "#00A087",
  "Candidatus Melainabacteria" = "#91D1C2",
  Methanobacteriota = "#bc3c29",
  Nitrospirota = "#7E6148",
  Bacteroidota = "#E6AB02",
  Thermodesulfobacteriota = "#6A3D9A",
  Acidobacteriota = "#4DBBD5",
  Myxococcota = "#F39B7F",
  Planctomycetota = "#8491B4",
  Verrucomicrobiota = "#7ba8c4",
  Chloroflexota = "#B09C85",
  Gemmatimonadota = "#7570B3"
)
comm_pal <- phy_cols[intersect(names(phy_cols), as.character(mean_all))]
extra_phy <- setdiff(as.character(mean_all), names(comm_pal))
if (length(extra_phy) > 0) {
  spare <- c(
    "#DC0000", "#66A61E", "#E7298A", "#1B9E77", "#A6761D",
    "#D95F02", "#1F78B4", "#CAB2D6", "#B15928", "#8DD3C7"
  )
  comm_pal <- c(comm_pal, setNames(spare[seq_along(extra_phy)], extra_phy))
}
comm_pal["Other"] <- "grey75"
phy_items <- c(as.character(mean_all), "Other")

p_comm <- stack_left(
  vplot(comm, "phylum_plot", comm_pal, phy_items, "", 1, 4),
  pad_leg(phy_items, comm_pal, 3, "plain")
)

add_guild_legs <- function(p, prod, ox, pal) {
  leg_mg <- take_legend(leg_plot(prod, pal, "Methanogens", "#7f0000"))
  leg_mt <- take_legend(leg_plot(ox, pal, "Methanotrophs", "#0d47a1"))
  gap <- gtable(widths = unit(1, "mm"), heights = unit(5, "mm"))
  legs <- list(leg_mg, gap, leg_mt)
  leg_h <- do.call(unit.c, lapply(legs, function(g) sum(g[["heights"]])))
  leg_w <- do.call(unit.pmax, lapply(legs, function(g) sum(g[["widths"]])))
  leg_col <- gtable(widths = leg_w, heights = leg_h)
  for (i in seq_along(legs)) {
    leg_col <- gtable_add_grob(leg_col, legs[[i]], t = i, l = 1, clip = "off")
  }
  p_grob <- ggplotGrob(p)
  out <- gtable(
    widths = unit.c(unit(1, "null"), unit(2.5, "mm"), leg_w, unit(2, "mm")),
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
  out
}

methane_plot <- function(habs, title, y_acc) {
  dat <- ch %>% filter(habitat %in% habs)
  prod <- named_order(dat, "Methanogens")
  ox <- named_order(dat, "Methanotrophs")
  breaks <- c(prod, ox)
  dat <- dat %>% filter(taxon_plot %in% breaks)
  p <- hplot(dat, "taxon_plot", ch_pal, breaks, title, y_acc, 4)
  add_guild_legs(p, prod, ox, ch_pal)
}

p_sed <- methane_plot(sed_habs, "Sediment", 1)
p_wat <- methane_plot(water_habs, "Water", 0.1)

save_one(p_comm, "community.png", 8.6, 16.2)
save_one(p_sed, "methane_sediment.png", 12, 5.6)
save_one(p_wat, "methane_water.png", 12, 6)

right <- plot_grid(
  p_wat,
  NULL,
  p_sed,
  ncol = 1,
  labels = c("B", "", "C"),
  label_size = 22,
  label_fontface = "bold",
  rel_heights = c(1.1, 0.03, 1)
)
combo <- frame_plot(
  plot_grid(
    p_comm,
    NULL,
    right,
    ncol = 3,
    rel_widths = c(1.04, 0.02, 1.18)
  ),
  top = 9,
  right = 5,
  bottom = 8,
  left = 5
)
save_one(combo, "combined.png", 18.6, 15.4)

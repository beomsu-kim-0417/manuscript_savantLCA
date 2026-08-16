#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

library(showtext)
library(sysfonts)
font_add("Arial",
         regular    = "/System/Library/Fonts/Supplemental/Arial.ttf",
         bold       = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         italic     = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
showtext_auto()
showtext_opts(dpi = 72)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing argument: ", flag)
  args[[i + 1L]]
}
source_dir <- normalizePath(arg_value("--analysis-source"), mustWork = TRUE)
fig_dir <- normalizePath(arg_value("--fig-dir"), mustWork = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  ink = "#262626", grey = "#737373", light = "#D9D9D9",
  blue = "#3C78A8", teal = "#4C9B8F", orange = "#D38B42",
  red = "#C65A52", purple = "#8073AC"
)
class_cols <- c(Memory = palette[["blue"]], Visuospatial = palette[["teal"]],
                Academic = palette[["orange"]], Arts = palette[["purple"]])

theme_release <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = palette[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = palette[["ink"]]),
      axis.text = element_text(colour = palette[["ink"]]),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.2, face = "bold"),
      plot.title = element_text(size = 7.2, face = "bold"),
      plot.subtitle = element_text(size = 6.0, colour = palette[["grey"]]),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}
theme_set(theme_release())

save_pub <- function(plot, stem, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  # PDF output uses outlined glyphs.
  showtext_opts(dpi = 72)
  showtext_auto(TRUE)
  pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = w, height = h, useDingbats = FALSE)
  print(plot); dev.off()

  # SVG retains editable text; PNG uses ragg.
  showtext_auto(FALSE)
  svglite(file.path(fig_dir, paste0(stem, ".svg")), width = w, height = h)
  print(plot); dev.off()

  agg_png(file.path(fig_dir, paste0(stem, ".png")), width = w, height = h,
          units = "in", res = 300)
  print(plot); dev.off()
  showtext_auto(TRUE)
}

read_candidate <- function(name) fread(file.path(source_dir, name), sep = "\t")
domain_levels <- c("Memory", "Reading", "Visuospatial", "Computational", "Music", "Drawing")
pgs_levels <- c("Autism", "Educational attainment", "Intelligence", "Schizophrenia", "ADHD",
                "Bipolar disorder", "Major depression", "Obsessive-compulsive disorder",
                "Epilepsy", "Seizure")

class_counts <- read_candidate("release_lca_class_counts_v2.9.tsv")
lca_fit <- read_candidate("release_lca_fit_reference_v2.9.tsv")
lca_profile <- read_candidate("release_lca_profiles_reference_v2.9.tsv")
item_primary <- read_candidate("release_item_pgs_ados_primary_pc_v3.0.tsv")
item_srs <- read_candidate("release_item_pgs_srs_pc_sensitivity_v3.0.tsv")
lca_primary <- read_candidate("release_lca_pgs_ados_primary_pc_v3.0.tsv")
dnv <- read_candidate("release_dnv_burden_ados_primary_no_pc_v3.0.tsv")

# Figure 1
item_primary[, domain := factor(domain, levels = domain_levels)]
item_primary[, PGS := factor(PGS, levels = pgs_levels)]
item_primary[, log2OR := log2(OR)]
item_primary[, sig := fifelse(p_fdr < 0.05, "*", "")]
p1a <- ggplot(item_primary, aes(domain, PGS, fill = log2OR)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig), size = 2.8, fontface = "bold") +
  scale_fill_gradient2(low = "#477DAF", mid = "white", high = "#C65A52", midpoint = 0,
                       name = expression(log[2](OR))) +
  labs(title = "Item-level polygenic-score associations",
       subtitle = "Primary ADOS CSS + VABS + PC1-PC5 model; * within-analysis FDR < 0.05",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

focus_domains <- data.table(
  domain = c("Reading", "Reading"),
  PGS = c("Educational attainment", "Intelligence")
)
axis_rows <- function(x, label) merge(x, focus_domains, by = c("domain", "PGS"))[, axis := label]
robust <- rbindlist(list(
  axis_rows(item_primary, "Primary"),
  axis_rows(item_srs, "SRS + VABS")
), use.names = TRUE, fill = TRUE)
robust[, signal := paste(domain, PGS, sep = " — ")]
robust[, signal := factor(signal, levels = rev(c(
  "Reading — Educational attainment", "Reading — Intelligence"
)))]
robust[, axis := factor(axis, levels = c("Primary", "SRS + VABS"))]
p1b <- ggplot(robust, aes(OR, signal, colour = axis)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3, colour = palette[["grey"]]) +
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), orientation = "y", width = 0.15,
                position = position_dodge(width = 0.52), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.52), size = 1.6) +
  scale_colour_manual(values = c("Primary" = palette[["ink"]],
                                 "SRS + VABS" = palette[["orange"]]), name = NULL) +
  coord_cartesian(xlim = c(1.0, 2.02)) +
  labs(title = "Planned sensitivity analyses",
       x = "Odds ratio per 1 s.d. PS (95% CI)", y = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.8))

fig1 <- (p1a | p1b) + plot_layout(widths = c(1.18, 1)) + plot_annotation(tag_levels = "a")
save_pub(fig1, "Figure1_v3.0", 183, 90)

# Figure 2
lca_profile[, domain := factor(domain, levels = domain_levels)]
lca_profile[, class_label := factor(class_label, levels = names(class_cols))]
p2a <- ggplot(lca_profile, aes(domain, endorsement_probability, colour = class_label, group = class_label)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.7) +
  scale_colour_manual(values = class_cols, name = "Latent class") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Four-class savant profiles", x = NULL, y = "Endorsement probability") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

class_counts[, class_label := factor(class_label, levels = names(class_cols))]
class_counts[, pct := 100 * N / sum(N)]
p2b <- ggplot(class_counts, aes(class_label, N, fill = class_label)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, pct)), vjust = -0.2, size = 2.0, lineheight = 0.9) +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(title = "Class sizes", x = NULL, y = "Individuals") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fit_long <- melt(lca_fit[, .(classes, BIC, aBIC)], id.vars = "classes",
                 variable.name = "metric", value.name = "value")
p2c <- ggplot(fit_long, aes(classes, value, colour = metric)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.5) +
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.3, colour = palette[["grey"]]) +
  annotate("text", x = 4.12, y = max(fit_long$value) - 45,
           label = "Four-class\ninterpretive model", hjust = 0, size = 2.0, colour = palette[["grey"]]) +
  scale_colour_manual(values = c(BIC = palette[["red"]], aBIC = palette[["teal"]]), name = NULL) +
  scale_x_continuous(breaks = 1:6) +
  labs(title = "Information criteria", x = "Number of latent classes", y = "Criterion") +
  theme(legend.position = "bottom")

fig2 <- (p2a | p2b | p2c) + plot_layout(widths = c(1.2, 0.85, 0.9)) + plot_annotation(tag_levels = "a")
save_pub(fig2, "Figure2_v3.0", 183, 70)

# Figure 3
lca_focus <- lca_primary[PGS %in% c("Autism", "Educational attainment", "Intelligence")]
lca_focus[, comparison := factor(comparison, levels = c("Visuospatial vs Memory", "Academic vs Memory", "Arts vs Memory"))]
lca_focus[, PGS := factor(PGS, levels = c("Autism", "Educational attainment", "Intelligence"))]
lca_focus[, is_fdr := p_fdr < 0.05]
p3a <- ggplot(lca_focus, aes(OR, comparison, colour = PGS, alpha = is_fdr)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3, colour = palette[["grey"]]) +
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), orientation = "y", width = 0.14,
                linewidth = 0.38, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.55, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Autism" = palette[["blue"]],
                                 "Educational attainment" = palette[["orange"]],
                                 "Intelligence" = palette[["red"]]), name = "PS") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.42), guide = "none") +
  labs(title = "Latent-class polygenic-score associations",
       subtitle = "Primary ADOS CSS + VABS + PC1-PC5 model", x = "Odds ratio per 1 s.d. PS (95% CI)", y = NULL)

dnv_min <- dnv[, .SD[which.min(p_value)], by = domain]
dnv_min[, domain := factor(domain, levels = rev(domain_levels))]
p3b <- ggplot(dnv_min, aes(OR, domain)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3, colour = palette[["grey"]]) +
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), orientation = "y", width = 0.14,
                linewidth = 0.4, colour = palette[["grey"]]) +
  geom_point(size = 1.6, colour = palette[["ink"]]) +
  coord_cartesian(xlim = c(0.6, 1.85)) +
  labs(title = "DNV burden screen", subtitle = "Minimum P by domain; no FDR-significant test",
       x = "Odds ratio per burden unit (95% CI)", y = NULL)

fig3 <- (p3a | p3b) + plot_layout(widths = c(1.15, 0.85)) + plot_annotation(tag_levels = "a")
save_pub(fig3, "Figure3_v3.0", 183, 72)


writeLines(capture.output(sessionInfo()), file.path(fig_dir, "figure_session_info_v3.0.txt"))

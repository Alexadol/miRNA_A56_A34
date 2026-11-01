#!/usr/bin/env Rscript
# fig_arm_asymmetry.R
#
# Two-panel figure of arm-usage asymmetry.
#   a  distribution of the dominant-arm share across loci
#   b  the same, resolved by family, with discordant loci highlighted
#
# Input files (produced earlier in the pipeline):
#   arm_group_counts.tsv       per-arm read counts, one column per library
#   arm_group_members.tsv      membership of each arm group (to drop shared arms)
#   annotation_functional.tsv  functional-arm assignment, for family and
#                              discordance; lines beginning with # are a legend
#
# Output: Fig_arm_asymmetry.svg and .png
#
# Requires: tidyverse, patchwork
#   install.packages(c("tidyverse", "patchwork"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

MIN_READS <- 50          # minimum total reads across both arms

## ---- load -----------------------------------------------------------------
counts <- read_tsv("arm_group_counts.tsv", show_col_types = FALSE)
members <- read_tsv("arm_group_members.tsv", show_col_types = FALSE)
func <- read_tsv("annotation_functional.tsv", comment = "#", show_col_types = FALSE)

lib_cols <- setdiff(names(counts), "group")
tot <- counts %>%
  mutate(total = rowSums(across(all_of(lib_cols)))) %>%
  select(group, total)

## arms that share a sequence with a paralogue cannot be quantified per locus
shared <- members %>% filter(n_members > 1) %>% pull(group)

## one row per locus with its two arm totals
arm_tbl <- tot %>%
  separate(group, into = c("cluster", "side"),
           sep = "\\.(?=mature$|star$)", remove = FALSE) %>%
  filter(!group %in% shared) %>%
  select(cluster, side, total) %>%
  pivot_wider(names_from = side, values_from = total) %>%
  filter(!is.na(mature), !is.na(star), mature + star >= MIN_READS) %>%
  mutate(dominant = pmax(mature, star) / (mature + star),
         coverage = mature + star)

## family and discordance from the functional annotation
meta <- func %>%
  transmute(cluster,
            family = ifelse(is.na(family) | family == "", "unassigned", family),
            discordant = functional_strand == "star")

loci <- arm_tbl %>% left_join(meta, by = "cluster")

med_all <- median(loci$dominant) * 100
message(sprintf("loci: %d; median dominant-arm share: %.1f%%",
                nrow(loci), med_all))

## ---- panel a: histogram ---------------------------------------------------
pa <- ggplot(loci, aes(x = dominant * 100)) +
  geom_histogram(breaks = seq(50, 100, 5),
                 fill = "#8E8C84", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = med_all, colour = "#D85A30",
             linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = med_all, y = Inf, vjust = 1.5,
           label = sprintf("median %.0f%%", med_all),
           colour = "#D85A30", size = 3.4) +
  scale_x_continuous(breaks = seq(50, 100, 10)) +
  labs(x = "reads on the dominant arm (%)", y = "number of loci", tag = "a") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

## ---- panel b: by family ---------------------------------------------------
## families with two or more loci shown individually; the rest pooled
multi <- loci %>% count(family) %>% filter(n >= 2, family != "unassigned") %>% pull(family)
loci_b <- loci %>%
  mutate(group_lab = case_when(
    family %in% multi            ~ family,
    family == "unassigned"       ~ "unassigned",
    TRUE                         ~ "single-locus families"))

fam_order <- loci_b %>%
  group_by(group_lab) %>%
  summarise(med = median(dominant), .groups = "drop") %>%
  arrange(med) %>% pull(group_lab)
loci_b$group_lab <- factor(loci_b$group_lab, levels = fam_order)

fam_med <- loci_b %>%
  group_by(group_lab) %>%
  summarise(med = median(dominant) * 100, .groups = "drop")

pb <- ggplot(loci_b, aes(x = dominant * 100, y = group_lab)) +
  geom_segment(data = fam_med,
               aes(x = med, xend = med,
                   y = as.numeric(factor(group_lab, levels = fam_order)) - 0.3,
                   yend = as.numeric(factor(group_lab, levels = fam_order)) + 0.3),
               colour = "#B4B2A9", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(size = coverage, colour = discordant), alpha = 0.6) +
  scale_colour_manual(values = c(`TRUE` = "#D85A30", `FALSE` = "#8E8C84"),
                      labels = c(`TRUE` = "conserved sequence on the minor arm",
                                 `FALSE` = "on the dominant arm"),
                      name = NULL) +
  scale_size_continuous(range = c(2, 7), trans = "log10", guide = "none") +
  scale_x_continuous(breaks = seq(50, 100, 10)) +
  labs(x = "reads on the dominant arm (%)", y = NULL, tag = "b") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

## ---- combine --------------------------------------------------------------
fig <- pa | pb + plot_layout(widths = c(1, 1.1))

ggsave("Fig_arm_asymmetry.svg", fig, width = 10.5, height = 5.2)
ggsave("Fig_arm_asymmetry.png", fig, width = 10.5, height = 5.2, dpi = 300)
message("written: Fig_arm_asymmetry.svg and .png")

#!/usr/bin/env Rscript
# fig2_size_classes.R
#
# Figure 2: read-length composition of the functional arm of each locus.
#   (a) loci whose functional arm is predominantly 21 nt
#   (b) loci dominated by another size class
#   (c) summary count of the dominant class across all loci (drawn as a small
#       table beside panel b in the manuscript; here printed to the console)
#
# Input
#   size_by_arm_long.tsv       cluster, arm_label, arm_5p3p, read_length, count
#   annotation_functional.tsv  functional strand per locus (comment lines #)
#   mirna_names.tsv            names (comment lines #)
#
# Requires: tidyverse

suppressPackageStartupMessages(library(tidyverse))

size <- read_tsv("size_by_arm_long.tsv", show_col_types = FALSE)
func <- read_tsv("annotation_functional.tsv", comment = "#", show_col_types = FALSE)
nm   <- read_tsv("mirna_names.tsv", comment = "#", show_col_types = FALSE)

## use the functional arm; where none is assigned, fall back to the mature arm
fstrand <- func %>%
  transmute(cluster,
            use = if_else(functional_strand %in% c("mature", "star"),
                          functional_strand, "mature"))

## six size classes
bin_len <- function(L) case_when(L <= 19 ~ "<=19", L >= 24 ~ ">=24",
                                 TRUE ~ as.character(L))
lvls <- c("<=19", "20", "21", "22", "23", ">=24")

arm <- size %>%
  inner_join(fstrand, by = "cluster") %>%
  filter(arm_label == use) %>%
  mutate(class = factor(bin_len(read_length), levels = lvls)) %>%
  group_by(cluster, class) %>%
  summarise(count = sum(count), .groups = "drop")

totals <- arm %>% group_by(cluster) %>% summarise(total = sum(count), .groups = "drop")
arm <- arm %>% left_join(totals, by = "cluster") %>% mutate(pct = 100 * count / total)

## dominant class per locus
dom <- arm %>% group_by(cluster) %>% slice_max(count, n = 1, with_ties = FALSE) %>%
  transmute(cluster, dominant = class, total)

## names and family, group A vs B
meta <- nm %>% transmute(cluster, name = coalesce(name, precursor_name), family)
dom  <- dom %>% left_join(meta, by = "cluster") %>%
  mutate(group = if_else(dominant == "21", "A", "B"))

message("Panel c summary (dominant class, all loci):")
print(count(dom, dominant))

cols <- c("<=19" = "#D9D7CE", "20" = "#A7C5A0", "21" = "#2E6E4E",
          "22" = "#C0472B", "23" = "#E0A06A", ">=24" = "#8E8C84")

make_panel <- function(which_group, title) {
  cl <- dom %>% filter(group == which_group) %>% arrange(family, name) %>% pull(cluster)
  d  <- arm %>% filter(cluster %in% cl) %>%
    left_join(meta, by = "cluster") %>%
    mutate(cluster = factor(cluster, levels = cl))
  ggplot(d, aes(cluster, pct, fill = class)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = cols, name = "read length (nt)") +
    scale_x_discrete(labels = function(x) meta$name[match(x, meta$cluster)]) +
    labs(title = title, x = NULL, y = "reads (%)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid.major.x = element_blank(),
          legend.position = "top")
}

pa <- make_panel("A", "a  Functional arm predominantly 21 nt")
pb <- make_panel("B", "b  Functional arm not predominantly 21 nt")

ggsave("Fig2a_size_21nt.pdf",  pa, width = 11, height = 4.5)
ggsave("Fig2b_size_other.pdf", pb, width = 7,  height = 4.5)
message("written: Fig2a_size_21nt.pdf, Fig2b_size_other.pdf")

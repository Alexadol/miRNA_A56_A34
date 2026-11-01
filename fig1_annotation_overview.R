#!/usr/bin/env Rscript
# fig1_annotation_overview.R
#
# Figure 1: pea root microRNA annotation in the context of published datasets.
#   (a) sampled tissues per study, drawn by hand or in a vector editor
#   (b) UpSet plot of sequence-level intersections across the four datasets
#   (c) composition of the root annotation by locus count and read share
#
# Panel (a) is a schematic and is not generated here. Panels (b) and (c) are
# produced from the tables below.
#
# Input
#   byloc.tsv        per-locus recovery in each published dataset (yes/no),
#                    from the set-comparison script
#   annotation_functional.tsv, mirna_names.tsv   for family and read share
#
# Requires: tidyverse, ComplexUpset (for panel b), patchwork

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexUpset)
  library(patchwork)
})

## ---- panel b: UpSet -------------------------------------------------------
## byloc.tsv has one row per root locus and a yes/no column per published set.
byloc <- read_tsv("byloc.tsv", show_col_types = FALSE) %>%
  mutate(across(c(Kreplak_2019, Mazurier_2022, Khan_2023),
                ~ .x == "yes"))

sets <- c("Kreplak_2019", "Mazurier_2022", "Khan_2023")
pb <- ComplexUpset::upset(
  byloc, sets,
  name = "dataset",
  min_size = 1,
  width_ratio = 0.2,
  base_annotations = list('Intersection size' =
    intersection_size(counts = TRUE)))

## ---- panel c: composition -------------------------------------------------
nm   <- read_tsv("mirna_names.tsv", comment = "#", show_col_types = FALSE)
func <- read_tsv("annotation_functional.tsv", comment = "#", show_col_types = FALSE)
## reads per locus must be supplied; here taken from a counts column if present
## (adjust to your file). Families contributing >= 1 % of reads are shown.

## This panel is a paired horizontal bar (loci vs read share). The exact
## geometry is given in the manuscript; build it from a per-family summary:
##   family, n_loci, pct_reads
## then two geom_col() layers mirrored around a central axis. Left as a
## template because the read-share column depends on your counts file.

message("Panel b written; panel c is a template, fill the per-family read share.")
ggsave("Fig1b_upset.pdf", pb, width = 7, height = 4)

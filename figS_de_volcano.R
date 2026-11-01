#!/usr/bin/env Rscript
# figS_de_volcano.R
#
# Supplementary volcano figure: differential arm-group abundance for the three
# contrasts, points coloured by functional / partner / unassigned arm.
# miR396b is labelled in the two mock contrasts.
#
# Input
#   DE_arm_groups_summary.tsv   DESeq2 output (group, log2FoldChange, pvalue,
#                               padj, contrast)
#   annotation_functional.tsv   functional strand per locus (comment lines #)
#   mirna_names.tsv             names (comment lines #)
#
# Requires: tidyverse, patchwork

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

de   <- read_tsv("DE_arm_groups_summary.tsv", show_col_types = FALSE)
func <- read_tsv("annotation_functional.tsv", comment = "#", show_col_types = FALSE)
nm   <- read_tsv("mirna_names.tsv", comment = "#", show_col_types = FALSE)

fstrand <- func %>% select(cluster, functional_strand)
names_t <- nm %>% transmute(cluster, pname = coalesce(name, precursor_name))

d <- de %>%
  mutate(cluster = str_replace(str_remove(group, "\\.(mature|star).*$"), "\\+.*$", ""),
         side    = str_extract(group, "(mature|star)")) %>%
  left_join(fstrand, by = "cluster") %>%
  left_join(names_t, by = "cluster") %>%
  mutate(role = case_when(
           is.na(functional_strand) | !functional_strand %in% c("mature","star") ~ "unassigned",
           functional_strand == side ~ "guide",
           TRUE ~ "passenger"),
         neglogp = -log10(pvalue),
         arm = str_remove(pname, "-[35]p$"),
         label = paste0(str_remove(arm, "^psa-"), " ", if_else(side == "mature", "", "")))

contrasts <- c(A34_vs_Con = "effective strain vs mock",
               A56_vs_Con = "nodC mutant vs mock",
               A34_vs_A56 = "effective strain vs nodC mutant")
d$contrast <- factor(d$contrast, levels = names(contrasts))
cols <- c(guide = "#1D7A54", passenger = "#E8622A", unassigned = "#B9B7AE")

## loci to label: both arms of miR396b in the two mock contrasts
lab <- d %>% filter(str_detect(pname, "MIR396b"),
                    contrast %in% c("A34_vs_Con", "A56_vs_Con"), pvalue < 0.05)

panel <- function(ct) {
  dd <- filter(d, contrast == ct) %>% arrange(role != "unassigned")
  ll <- filter(lab, contrast == ct)
  ggplot(dd, aes(log2FoldChange, neglogp, colour = role)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "#8E8C84", linewidth = 0.4) +
    geom_point(aes(size = role, alpha = role)) +
    ggrepel::geom_text_repel(data = ll, aes(label = str_remove(pname, "^psa-")),
                             size = 4, colour = "black", min.segment.length = 0) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_size_manual(values = c(guide = 1.5, passenger = 1.5, unassigned = 1),
                      guide = "none") +
    scale_alpha_manual(values = c(guide = 0.8, passenger = 0.8, unassigned = 0.4),
                       guide = "none") +
    coord_cartesian(xlim = c(-4.2, 4.2), ylim = c(0, 5.2)) +
    labs(title = contrasts[[ct]], x = "log2 fold change", y = "-log10 P") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 13, hjust = 0.5),
          panel.grid.minor = element_blank())
}

fig <- (panel("A34_vs_Con") | panel("A56_vs_Con") | panel("A34_vs_A56")) +
  plot_annotation(tag_levels = "a")

ggsave("FigS_de_volcano.svg", fig, width = 11.5, height = 4.4)
ggsave("FigS_de_volcano.png", fig, width = 11.5, height = 4.4, dpi = 300)
message("written: FigS_de_volcano.svg and .png (needs ggrepel for labels)")

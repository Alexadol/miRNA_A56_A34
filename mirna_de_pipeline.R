# =============================================================================
# Differential expression analysis of plant miRNAs from ShortStack4 output
#
# Input:  Counts_404.txt and Results.txt from ShortStack4
#         cluster_annotations.tsv with columns: Cluster, family
# Output: per-contrast DE tables (MLE + ashr-shrunken LFCs), QC plots,
#         family-level dotplot/barplot, paralog profile plots, FASTA
#
# ashr is preferred over apeglm for small-n datasets where apeglm's
# empirical Bayes prior over-shrinks real signals (e.g. miR164d MLE = -2.0
# shrunk to -0.04 with apeglm). ashr also accepts contrast=, so all three
# contrasts come from a single fit.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ashr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

set.seed(42)

# ---- Configuration ---------------------------------------------------------

INPUT_DIR  <- "/home/alexadol/Work/miRNA"
OUTPUT_DIR <- file.path(INPUT_DIR, "mirna_DE_results")
ANNOT_PATH <- file.path(INPUT_DIR, "cluster_annotations.tsv")

# Technical outliers identified from PCA and sample-distance heatmap
OUTLIERS <- c("") #Add outliers

MIN_TOTAL_COUNT   <- 10
ALPHA             <- 0.05
SHRINK_TYPE       <- "ashr"
INCLUDE_REPLICATE <- FALSE   # set TRUE only if _N indices reflect prep batches

# Visualization thresholds (looser than ALPHA: exploration, not significance)
PVAL_CUT_VIS    <- 0.05
LFC_CUT_VIS     <- 0.5
INCLUDE_UNKNOWN <- TRUE
MIN_PARALOGS    <- 2          # min paralogs to show in profile plot

DIR_COLORS    <- c(UP = "#D7191C", DOWN = "#2C7BB6", NS = "grey85")
COND_COLORS   <- c(Cont = "#E41A1C", A34 = "#4DAF4A", A56 = "#377EB8")
STRIPE_COLORS <- c("white", "grey94")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ---------------------------------------------------------------

# Coerce list-columns (returned by some merge paths) to atomic character.
to_character <- function(x) {
  if (is.list(x)) {
    vapply(x, function(v) {
      if (is.null(v) || length(v) == 0L) NA_character_
      else as.character(v)[[1L]]
    }, character(1L))
  } else as.character(x)
}

# Build a results table with both MLE and shrunken LFCs for one contrast.
build_results <- function(dds_obj, contrast_vec, results_meta) {
  res_mle <- results(dds_obj, contrast = contrast_vec, alpha = ALPHA)
  res_shr <- lfcShrink(dds_obj, contrast = contrast_vec,
                       type = SHRINK_TYPE, quiet = TRUE)

  meta_keep <- intersect(c("Name", "Coords", "MajorRNA",
                           "DicerCall", "Length", "Reads"),
                         colnames(results_meta))

  out <- data.frame(
    Name          = rownames(res_mle),
    baseMean      = res_mle$baseMean,
    log2FC_MLE    = res_mle$log2FoldChange,
    lfcSE_MLE     = res_mle$lfcSE,
    log2FC_shrunk = res_shr$log2FoldChange,
    pvalue        = res_mle$pvalue,
    padj          = res_mle$padj,
    stringsAsFactors = FALSE
  )
  out <- merge(out, results_meta[, meta_keep], by = "Name", sort = FALSE)
  out <- out[order(out$padj), ]

  # Flag loci where shrinkage drastically diverges from MLE (low confidence).
  out$shrinkage_flag <- ifelse(
    !is.na(out$log2FC_MLE) & !is.na(out$log2FC_shrunk) &
      abs(out$log2FC_MLE) > 0.5 &
      abs(out$log2FC_MLE / pmax(abs(out$log2FC_shrunk), 0.01)) > 4,
    "HIGH_SE_unstable", ""
  )

  list(mle = res_mle, shrunken = res_shr, table = out)
}

# ---- Load and subset to confident miRNAs (MIRNA = Y) -----------------------

counts_raw   <- read.table(file.path(INPUT_DIR, "Counts_404.txt"),
                           header = TRUE, sep = "\t", check.names = FALSE,
                           stringsAsFactors = FALSE)
results_meta <- read.table(file.path(INPUT_DIR, "Results.txt"),
                           header = TRUE, sep = "\t", check.names = FALSE,
                           stringsAsFactors = FALSE)

mirna <- if ("MIRNA" %in% colnames(counts_raw))
  counts_raw[counts_raw$MIRNA == "Y", , drop = FALSE] else counts_raw

meta_cols   <- intersect(c("Coords", "Name", "MIRNA"), colnames(mirna))
sample_cols <- setdiff(colnames(mirna), meta_cols)

count_mat <- as.matrix(mirna[, sample_cols])
rownames(count_mat) <- mirna$Name
mode(count_mat) <- "integer"

# ---- Sample metadata and outlier removal -----------------------------------

coldata <- data.frame(
  sample    = sample_cols,
  condition = factor(sub("_.*", "", sample_cols),
                     levels = c("Cont", "A34", "A56")),
  replicate = factor(sub("^[^_]+_([0-9]+)_.*", "\\1", sample_cols)),
  row.names = sample_cols
)

keep <- !Reduce("|", lapply(OUTLIERS, function(p) grepl(p, sample_cols)))
count_mat <- count_mat[, keep]
coldata   <- coldata[keep, ]
coldata$condition <- droplevels(coldata$condition)
coldata$replicate <- droplevels(coldata$replicate)
coldata$condition <- relevel(coldata$condition, ref = "Cont")

stopifnot(all(rownames(coldata) == colnames(count_mat)))

cat("Samples per condition after outlier removal:\n")
print(table(coldata$condition))

# ---- Fit DESeq2 ------------------------------------------------------------

design_formula <- if (INCLUDE_REPLICATE) ~ replicate + condition else ~ condition

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData   = coldata,
                              design    = design_formula)
dds <- dds[rowSums(counts(dds)) >= MIN_TOTAL_COUNT, ]
dds <- DESeq(dds)

# ---- QC plots --------------------------------------------------------------

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("condition", "replicate"),
                    returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition,
                     shape = replicate, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  xlab(sprintf("PC1: %d%% variance", pct_var[1])) +
  ylab(sprintf("PC2: %d%% variance", pct_var[2])) +
  scale_color_manual(values = COND_COLORS) +
  theme_bw()
ggsave(file.path(OUTPUT_DIR, "QC_PCA.pdf"), width = 7, height = 5)

sample_dist <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dist),
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = coldata[, "condition", drop = FALSE],
         annotation_row = coldata[, "condition", drop = FALSE],
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         filename = file.path(OUTPUT_DIR, "QC_sample_distances.pdf"),
         width = 7, height = 6)

pdf(file.path(OUTPUT_DIR, "QC_dispersion.pdf"), width = 6, height = 5)
plotDispEsts(dds)
dev.off()

pdf(file.path(OUTPUT_DIR, "QC_cooks_distance.pdf"), width = 7, height = 5)
par(mar = c(8, 5, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), las = 2,
        ylab = "log10(Cook's distance)")
dev.off()

# ---- Three contrasts -------------------------------------------------------

contrasts_list <- list(
  A34_vs_Cont = build_results(dds, c("condition", "A34", "Cont"), results_meta),
  A56_vs_Cont = build_results(dds, c("condition", "A56", "Cont"), results_meta),
  A34_vs_A56  = build_results(dds, c("condition", "A34", "A56"),  results_meta)
)

for (nm in names(contrasts_list)) {
  tab <- contrasts_list[[nm]]$table
  write.table(tab, file.path(OUTPUT_DIR, paste0("DE_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  sig <- tab[!is.na(tab$padj) & tab$padj < ALPHA, , drop = FALSE]
  write.table(sig, file.path(OUTPUT_DIR, paste0("DE_", nm, "_significant.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  pdf(file.path(OUTPUT_DIR, paste0("MA_", nm, ".pdf")), width = 6, height = 5)
  DESeq2::plotMA(contrasts_list[[nm]]$shrunken, ylim = c(-5, 5),
                 main = nm, alpha = ALPHA)
  dev.off()
}

# ---- Heatmap of top DE loci across all samples -----------------------------

top_set <- unique(unlist(lapply(contrasts_list,
                                function(x) head(x$table$Name, 15))))
top_set <- top_set[!is.na(top_set) & top_set %in% rownames(vsd)]

if (length(top_set) >= 2) {
  mat <- t(scale(t(assay(vsd)[top_set, , drop = FALSE])))
  ord <- order(coldata$condition, coldata$replicate)
  pheatmap(mat[, ord],
           annotation_col = coldata[ord, "condition", drop = FALSE],
           cluster_cols   = FALSE,
           show_rownames  = TRUE,
           fontsize_row   = 7,
           filename = file.path(OUTPUT_DIR, "top_DE_heatmap.pdf"),
           width = 7, height = max(5, length(top_set) * 0.18))
}

# ---- FASTA of significant miRNAs (for miRBase / PmiREN lookup) -------------

all_de <- bind_rows(lapply(names(contrasts_list), function(nm) {
  contrasts_list[[nm]]$table %>% mutate(contrast = nm)
}))

sig_loci <- all_de %>%
  filter(!is.na(padj), padj < ALPHA, !is.na(MajorRNA)) %>%
  distinct(Name, MajorRNA)

if (nrow(sig_loci) > 0) {
  fa_lines <- c(rbind(paste0(">", sig_loci$Name), sig_loci$MajorRNA))
  writeLines(fa_lines, file.path(OUTPUT_DIR, "significant_miRNAs.fasta"))
}

# ---- Family-level annotation -----------------------------------------------

ann <- read.table(ANNOT_PATH, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
ann$family  <- to_character(ann$family)
ann$Cluster <- to_character(ann$Cluster)
ann <- ann %>%
  distinct(Cluster, family) %>%
  group_by(Cluster) %>% slice(1) %>% ungroup() %>%
  as.data.frame()

de_annot <- merge(all_de %>% select(Name, log2FC_MLE, pvalue, contrast),
                  ann, by.x = "Name", by.y = "Cluster",
                  all.x = TRUE, sort = FALSE)
de_annot$family <- to_character(de_annot$family)
de_annot$family[is.na(de_annot$family) | de_annot$family == ""] <- "unknown"

de_annot$direction <- with(de_annot, ifelse(
  is.na(pvalue) | is.na(log2FC_MLE), "NS",
  ifelse(pvalue < PVAL_CUT_VIS &  log2FC_MLE >  LFC_CUT_VIS, "UP",
         ifelse(pvalue < PVAL_CUT_VIS & log2FC_MLE < -LFC_CUT_VIS, "DOWN",
                "NS"))
))
de_annot$contrast <- factor(de_annot$contrast,
                            levels = c("A34_vs_Cont", "A56_vs_Cont",
                                       "A34_vs_A56"))

if (!INCLUDE_UNKNOWN) de_annot <- de_annot[de_annot$family != "unknown", ]

# Order families by total DE activity (most-active on top)
fam_active <- de_annot[de_annot$direction %in% c("UP", "DOWN"), , drop = FALSE]
fam_order  <- names(sort(table(fam_active$family), decreasing = TRUE))
fam_order  <- c(fam_order, setdiff(unique(de_annot$family), fam_order))
de_annot$family <- factor(de_annot$family, levels = rev(fam_order))
fam_levels <- levels(de_annot$family)

stripe_df <- data.frame(
  family   = fam_levels,
  fill_col = STRIPE_COLORS[seq_along(fam_levels) %% 2 + 1]
)

write.table(de_annot, file.path(OUTPUT_DIR, "family_DE_per_locus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Family-level dotplot (per-locus LFCs with row stripes) ----------------

lfc_lim <- ceiling(max(abs(de_annot$log2FC_MLE), na.rm = TRUE) * 1.05)

p_dots <- ggplot(de_annot, aes(x = log2FC_MLE, y = family)) +
  geom_tile(data = stripe_df, inherit.aes = FALSE,
            aes(x = 0, y = family, fill = fill_col),
            width = Inf, height = 1) +
  scale_fill_identity() +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_CUT_VIS, LFC_CUT_VIS),
             linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_jitter(aes(color = direction, size = direction, alpha = direction),
              height = 0.22, width = 0) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_color_manual(values = DIR_COLORS,
                     breaks = c("UP", "DOWN", "NS"), name = NULL) +
  scale_size_manual( values = c(UP = 2.8, DOWN = 2.8, NS = 1.6),
                     guide = "none") +
  scale_alpha_manual(values = c(UP = 0.95, DOWN = 0.95, NS = 0.55),
                     guide = "none") +
  scale_x_continuous(limits = c(-lfc_lim, lfc_lim),
                     breaks = pretty_breaks(5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "log2 fold change (MLE)", y = NULL,
       title = "Per-locus log2 fold changes by miRNA family",
       subtitle = sprintf("UP/DOWN: raw p < %g and |LFC| > %g",
                          PVAL_CUT_VIS, LFC_CUT_VIS)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(face = "bold", size = 10),
        legend.position    = "bottom",
        strip.background   = element_rect(fill = "grey92"),
        strip.text         = element_text(face = "bold"),
        plot.title         = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "family_DE_dotplot.pdf"), p_dots,
       width = 11, height = max(4.5, length(fam_levels) * 0.36 + 1.5))

# ---- Family-level stacked barplot ------------------------------------------

fam_summary <- de_annot %>%
  dplyr::count(family, contrast, direction, name = "n_loci") %>%
  complete(family, contrast, direction, fill = list(n_loci = 0))
fam_summary$direction <- factor(fam_summary$direction,
                                levels = c("DOWN", "NS", "UP"))
fam_summary$family <- factor(fam_summary$family, levels = fam_levels)

p_stack <- ggplot(fam_summary,
                  aes(x = n_loci, y = family, fill = direction)) +
  geom_tile(data = stripe_df, inherit.aes = FALSE,
            aes(x = 0, y = family, fill = fill_col),
            width = Inf, height = 1) +
  scale_fill_identity() +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = DIR_COLORS, breaks = c("UP", "DOWN", "NS"),
    name = sprintf("Direction (raw p < %g, |LFC| > %g)",
                   PVAL_CUT_VIS, LFC_CUT_VIS)
  ) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_x_continuous(breaks = pretty_breaks(4),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of loci", y = NULL,
       title = "miRNA family DE summary across contrasts") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(face = "bold", size = 10),
        legend.position    = "bottom",
        strip.background   = element_rect(fill = "grey92"),
        plot.title         = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "family_DE_stacked_barplot.pdf"), p_stack,
       width = 9, height = max(4.5, length(fam_levels) * 0.36 + 1.5))

fam_wide <- pivot_wider(fam_summary, names_from = direction,
                        values_from = n_loci, values_fill = 0)
fam_wide <- fam_wide[order(fam_wide$contrast, fam_wide$family), ]
fam_wide <- fam_wide[, c("family", "contrast", "UP", "DOWN", "NS")]
write.table(fam_wide, file.path(OUTPUT_DIR, "family_DE_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Paralog profile plots (z-scored VST per locus, family-faceted) --------

sample_order <- rownames(coldata)[order(coldata$condition, coldata$replicate)]
vst_z <- t(scale(t(assay(vsd)[, sample_order])))

vst_long <- as.data.frame(vst_z) %>%
  mutate(Cluster = rownames(vst_z)) %>%
  pivot_longer(-Cluster, names_to = "sample", values_to = "value") %>%
  left_join(ann, by = "Cluster") %>%
  mutate(family    = ifelse(is.na(family) | family == "", "unknown", family),
         x         = match(sample, sample_order),
         condition = sub("_.*", "", sample))

# Restrict to families with multiple paralogs (single-locus families plot
# trivially as one line, providing little visual value in this comparison).
fam_counts <- vst_long %>% distinct(Cluster, family) %>% dplyr::count(family)
multi <- fam_counts$family[fam_counts$n >= MIN_PARALOGS]
vst_long <- vst_long[vst_long$family %in% multi, ]

family_mean <- vst_long %>%
  group_by(family, x) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Background condition-block shading
cond_seq <- coldata$condition[match(sample_order, rownames(coldata))]
cond_run <- rle(as.character(cond_seq))
cs <- cumsum(c(0, cond_run$lengths))
cond_rect <- data.frame(
  condition = cond_run$values,
  xmin      = cs[-length(cs)] + 0.5,
  xmax      = cs[-1] + 0.5
)

# Color paralogs by their direction in the A34_vs_A56 contrast
locus_dir <- de_annot %>%
  filter(contrast == "A34_vs_A56") %>%
  select(Cluster = Name, direction)
vst_long <- left_join(vst_long, locus_dir, by = "Cluster")
vst_long$direction[is.na(vst_long$direction)] <- "NS"

p_profile <- ggplot(vst_long, aes(x = x, y = value, group = Cluster)) +
  geom_rect(data = cond_rect, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf, fill = condition),
            alpha = 0.13) +
  geom_line(aes(color = direction), alpha = 0.65, linewidth = 0.55) +
  geom_point(aes(color = direction), alpha = 0.75, size = 1.3) +
  geom_line(data = family_mean, inherit.aes = FALSE,
            aes(x = x, y = mean_value, group = family),
            color = "black", linewidth = 1.1) +
  scale_fill_manual(values = COND_COLORS, name = "Condition") +
  scale_color_manual(values = DIR_COLORS,
                     breaks = c("UP", "DOWN", "NS"),
                     name = "DE in A34_vs_A56") +
  scale_x_continuous(breaks = seq_along(sample_order),
                     labels = sample_order, expand = c(0, 0)) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = NULL, y = "VST (z-score per locus)",
       title = "Expression profiles of miRNA paralogs by family",
       subtitle = sprintf("Bold black = family mean. Loci: %d in %d families.",
                          length(unique(vst_long$Cluster)), length(multi))) +
  theme_bw(base_size = 10) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold", size = 9),
        plot.title       = element_text(face = "bold"))

n_fam  <- length(multi)
ncol_p <- min(4, n_fam)
nrow_p <- ceiling(n_fam / ncol_p)
ggsave(file.path(OUTPUT_DIR, "family_profile_plots.pdf"), p_profile,
       width  = 3.5 * ncol_p + 1,
       height = 2.6 * nrow_p + 1.5,
       limitsize = FALSE)



# ---- Reproducibility -------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(OUTPUT_DIR, "sessionInfo.txt"))

cat("Done. Results in:", OUTPUT_DIR, "\n")

#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Title: run_mmuphin.R
# Description: Batch effect correction of MetaPhlAn abundance
# tables with MMUPHin::adjust_batch
#
# Author: Mike Martinez
# Lab: GDSC
# Project: GDSC-MGX-Pipeline
# Date created: 08/27/26
#
# Takes ONE rank, as written by stratify_abundance.py, not the stacked
# merged_abundance_table.txt. That file carries every rank at once, so it sums
# to ~100 per rank and a phylum row is its child species counted again. Both
# break adjust_batch. The column sum check below catches it.
#
# Note adjust_batch does NOT drop features, which is the obvious assumption and
# it is wrong. It models only the features passing its own prevalence filter and
# carries the rest through, renormalising each sample afterwards. The corrected
# table therefore has exactly the rows the input had.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
suppressMessages(library(MMUPHin))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(vegan))

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Covariates come last and are variadic, so an empty {params.covariates} in
# the Snakemake shell block simply disappears and the four required arguments
# still land in the right positions. Nothing needs quoting. Putting them
# anywhere else would mean an empty expansion silently shifting out_prefix.
if (length(args) < 4) {
  stop(paste(
    "Usage: Rscript run_mmuphin.R <abundance_table> <sample_sheet> <batch_col> <out_prefix> [covariate_cols ...]",
    "",
    "  abundance_table  one rank, from stratify_abundance.py",
    "  sample_sheet     needs sample_id plus the columns below",
    "  batch_col        technical variable to correct for, exactly one",
    "  out_prefix       path prefix for the corrected table and figure",
    "  covariate_cols   optional biological variable(s) to preserve. Takes",
    "                   separate arguments, a comma separated list, or neither.",
    "                   Each must be its own sample sheet column, MMUPHin builds",
    "                   a design matrix column per level so they cannot share",
    "                   one field.",
    sep = "\n"
  ))
}

abundanceFile <- args[1]
sampleSheet   <- args[2]
batchCol      <- args[3]
outPrefix     <- args[4]

#----- Accept 'group sex' and 'group,sex' alike, so it does not matter whether
# the config holds a list or a single string
covariateCols <- trimws(unlist(strsplit(args[-(1:4)], ",")))
covariateCols <- covariateCols[nzchar(covariateCols)]

#----- No covariates is a supported mode, adjust_batch takes covariates = NULL.
# It is also the risky one, so say what it costs rather than letting it pass
# silently: with nothing to hold onto, anything biological that happens to track
# the batch is removed along with the batch effect.
if (length(covariateCols) == 0) {
  message("No covariates given, correcting on '", batchCol, "' alone. Any ",
          "biological signal that tracks the batch will be removed with it.")
}

correctedFile <- paste0(outPrefix, "_mmuphin_corrected_pct.txt")
figureFile    <- paste0(outPrefix, "_mmuphin_pcoa.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA AND CHECK IT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- clade_name in column 1, one column per sample. Column 1 moves to rownames
# so what is left is the features x samples matrix adjust_batch wants.
abdRaw <- as.data.frame(data.table::fread(abundanceFile))
indexName <- colnames(abdRaw)[1]
abd <- as.matrix(abdRaw[, -1, drop = FALSE])
rownames(abd) <- abdRaw[[1]]

#----- fread types a column as character the moment one cell is not a number
if (!is.numeric(abd)) {
  stop("Non numeric or missing value(s) in ", abundanceFile)
}

#----- Prepare metadata and sanity check
meta <- as.data.frame(data.table::fread(sampleSheet))
if (!"sample_id" %in% colnames(meta)) {
  stop(sampleSheet, " has no 'sample_id' column")
}
rownames(meta) <- meta$sample_id

if (!setequal(rownames(meta), colnames(abd))) {
  stop("Sample names in ", abundanceFile, " do not match ", sampleSheet, ".\n",
       "  In sample sheet but not in table: ",
       paste(setdiff(rownames(meta), colnames(abd)), collapse = ", "), "\n",
       "  In table but not in sample sheet: ",
       paste(setdiff(colnames(abd), rownames(meta)), collapse = ", "))
}
meta <- meta[colnames(abd), , drop = FALSE]

#----- Build model dataframe and sanity check
modelCols <- c(batchCol, covariateCols)
absentCols <- setdiff(modelCols, colnames(meta))
if (length(absentCols) > 0) {
  stop("Column(s) not in ", sampleSheet, ": ", paste(absentCols, collapse = ", "),
       "\n  Available: ", paste(setdiff(colnames(meta), "sample_id"), collapse = ", "))
}

for (col in modelCols) {
  values <- trimws(as.character(meta[[col]]))
  unfilled <- is.na(values) | values == ""
  if (any(unfilled)) {
    stop("'", col, "' is blank for sample(s): ",
         paste(rownames(meta)[unfilled], collapse = ", "))
  }
  meta[[col]] <- factor(values)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GUARD THE CASES MMUPHIN CANNOT HANDLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

batchCounts <- table(meta[[batchCol]])
message("Batches: ", paste0(names(batchCounts), " n=", batchCounts, collapse = ", "))

#----- If there is only one batch, do no correction. 
if (length(batchCounts) < 2) {
  message("Only one level in '", batchCol, "', nothing to correct. ",
          "Writing the input through unchanged.")
  fwrite(abdRaw, correctedFile, sep = "\t")
  quit(save = "no", status = 0)
}

#----- adjust_batch cannot estimate parameters from a batch of one or confounded design.
if (any(batchCounts < 2)) {
  stop("Batch level(s) with a single sample: ",
       paste(names(batchCounts)[batchCounts < 2], collapse = ", "))
}

for (col in covariateCols) {
  crossTab <- table(meta[[batchCol]], meta[[col]])
  if (all(rowSums(crossTab > 0) == 1)) {
    stop("'", col, "' is perfectly confounded with '", batchCol,
         "': every batch holds only one level of it.\n",
         paste(capture.output(print(crossTab)), collapse = "\n"))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RUN MMUPHIN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Ensure a single taxonomic rank is being supplied, not full metaphlan profile
colTotals <- colSums(abd)
if (any(abs(colTotals - 100) > 5)) {
  stop("Sample columns in ", abundanceFile, " sum to ",
       paste(sprintf("%.1f", range(colTotals)), collapse = " - "), ", not ~100.\n",
       "  If this is merged_abundance_table.txt it holds every taxonomic rank at ",
       "once. Run stratify_abundance.py first and pass one rank.")
}
abdFrac <- abd / 100

#----- Run MMUPHin and extract corrected counts
fit <- MMUPHin::adjust_batch(
  feature_abd = abdFrac,           
  batch       = batchCol,
  covariates  = if (length(covariateCols) > 0) covariateCols else NULL,     
  data        = meta,
  control     = list(diagnostic_plot = paste0(outPrefix, "_mmuphin_diagnostic.pdf"),
                     verbose = TRUE)
)
corrected <- fit$feature_abd_adj

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ASSESS THE CORRECTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Before/after PCoA on Bray-Curtis is the quickest read on whether the
# correction did anything, and PERMANOVA puts a number on it: the batch R2
# should fall, the covariate R2 should not.
set.seed(1234)
distBefore <- vegdist(t(abdFrac), method = "bray")
distAfter  <- vegdist(t(corrected), method = "bray")

batchR2 <- function(d) {
  res <- try(adonis2(as.formula(paste("d ~", paste(modelCols, collapse = " + "))),
                     data = meta, by = "margin"), silent = TRUE)
  if (inherits(res, "try-error")) NA_real_ else as.data.frame(res)[batchCol, "R2"]
}
r2Before <- batchR2(distBefore)
r2After  <- batchR2(distAfter)
message(sprintf("%s R2: %.3f -> %.3f", batchCol, r2Before, r2After))

#----- Negative eigenvalues are normal for Bray-Curtis, so scale the percent
# variance against the positive ones only
ordinate <- function(d, label) {
  pcoa <- cmdscale(d, k = 2, eig = TRUE)
  pct <- 100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0])
  data.frame(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2],
             Panel = sprintf("%s (PCo1 %.1f%%, PCo2 %.1f%%)", label, pct[1], pct[2]))
}

ord <- rbind(ordinate(distBefore, "Before correction"),
             ordinate(distAfter, "After correction"))

#----- Keep the panel order, otherwise the factor sorts alphabetically and
# "After" lands on the left. Both blocks are in meta's row order.
ord$Panel <- factor(ord$Panel, levels = unique(ord$Panel))
ord$Batch <- rep(meta[[batchCol]], 2)

#----- Shape carries the first covariate, so there is nothing to shape by when
# none were given
hasCovariates <- length(covariateCols) > 0
if (hasCovariates) ord$Group <- rep(meta[[covariateCols[1]]], 2)

pcoaPlot <- ggplot(ord, aes(PCo1, PCo2, colour = Batch)) +
  (if (hasCovariates) geom_point(aes(shape = Group), size = 3.5, alpha = 0.9)
   else geom_point(size = 3.5, alpha = 0.9)) +
  facet_wrap(~ Panel, nrow = 1, scales = "free") +
  theme_classic(base_size = 15) +
  labs(
    title = "Bray-Curtis PCoA before and after MMUPHin correction",
    subtitle = sprintf("%d features, %d samples, %s R2 %.3f -> %.3f",
                       nrow(corrected), ncol(corrected), batchCol, r2Before, r2After),
    colour = batchCol, shape = if (hasCovariates) covariateCols[1] else NULL
  ) +
  scale_colour_brewer(palette = "Set2") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(figureFile, pcoaPlot, width = 12, height = 5.5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# WRITE OUTPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Back to percentages on the way out, so the corrected table reads with the
# same code as the uncorrected one. adjust_batch renormalises, so the columns
# sum to 100 again. The scale is in the filename.
out <- data.frame(rownames(corrected), corrected * 100, check.names = FALSE)
colnames(out)[1] <- indexName
fwrite(out, correctedFile, sep = "\t")

message("Wrote ", correctedFile)
message("Wrote ", figureFile)

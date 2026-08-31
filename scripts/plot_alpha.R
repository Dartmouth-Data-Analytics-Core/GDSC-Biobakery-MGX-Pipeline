#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: plot_alpha.R
# Description: Generate alpha diversity plots
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/16/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(vegan)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggrepel)

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
# The directories are passed in rather than assumed, because the pipeline
# writes under out_dir from config.yaml, which is not the directory snakemake
# runs from. Hardcoding "results/" here looked for them next to the Snakefile.
if (length(args) != 3) {
  stop(paste("Usage: Rscript plot_alpha.R <sample_sheet.csv> <abundance_dir> <plots_dir>",
             "  abundance_dir  metaphlan merged_abundances folder, read and written to",
             "  plots_dir      figures, plus the per-sample diversity tables",
             sep = "\n"))
}

#----- Set variables based on command line args
sampleSheet <- args[1]
abundDir  <- args[2]
plotDir      <- args[3]

#----- Make sure the paths end in a separator, since they are pasted onto
# filenames below
addSlash <- function(p) if (grepl("/$", p)) p else paste0(p, "/")
abundDir <- addSlash(abundDir)
plotDir     <- addSlash(plotDir)

#----- Create the figure directory if it is not already there
if (!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA, TIDY, AND ARRANGE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in the sample sheet
sampleList <- fread(sampleSheet) |>
    as.data.frame()

#----- Get vector of expected samples from the sample sheet
samples <- sampleList$sample_id

#----- Read a stratified abundance table written by stratify_abundance.py.
# These tables now carry real sample names in the header, so the names are
# checked against the sample sheet instead of being assigned from it. The old
# code overwrote the column names with the sample sheet order, which silently
# mislabelled every sample if the two orders ever disagreed.
readAbundance <- function(path, expected) {
  df <- fread(path) |>
    as.data.frame()
  rownames(df) <- df[,1]
  df[,1] <- NULL

  missing <- setdiff(expected, colnames(df))
  extra <- setdiff(colnames(df), expected)
  if (length(missing) > 0 || length(extra) > 0) {
    stop(paste0(
      "Sample names in ", path, " do not match the sample sheet.\n",
      "  In sample sheet but not in table: ", paste(missing, collapse = ", "), "\n",
      "  In table but not in sample sheet: ", paste(extra, collapse = ", ")
    ))
  }

  #----- Put columns in sample sheet order
  df[, expected, drop = FALSE]
}

#----- Read in the species only data
spec <- readAbundance(paste0(abundDir, "merged_species_abundance_table.txt"), samples)
fwrite(spec, file = paste0(abundDir, "named_merged_species_abundance_table.txt"),
       row.names = TRUE)

#----- Transpose the data (taxa as columns)
specT <- as.data.frame(t(spec))

#----- Read in the genus only data
genus <- readAbundance(paste0(abundDir, "merged_genus_abundance_table.txt"), samples)
fwrite(genus, file = paste0(abundDir, "named_merged_genus_abundance_table.txt"),
       row.names = TRUE)

#----- Transpose the data (taxa as columns)
genusT <- as.data.frame(t(genus))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALCULATE ALPHA DIVERSITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- On species
specShan <- vegan::diversity(specT, index = "shannon")
specSimp <- vegan::diversity(specT, index = "simpson")
specObs <- rowSums(specT > 0)

#----- Combine into dataframe
specAlpha <- data.frame(
    Sample = rownames(specT),
    Shannon = specShan,
    Simpson = specSimp,
    Observed = specObs
)

#----- On Genus
genShan <- vegan::diversity(genusT, index = "shannon")
genSimp <- vegan::diversity(genusT, index = "simpson")
genObs <- rowSums(genusT > 0)

#----- Combine into dataframe
genusAlpha <- data.frame(
    Sample = rownames(genusT),
    Shannon = genShan,
    Simpson = genSimp,
    Observed = genObs
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Pivot longer for plotting (species)
specAlphaLong <- specAlpha |>
  pivot_longer(cols = c(Shannon, Simpson, Observed),
               names_to = "Metric",
               values_to = "Value") |>
    as.data.frame()
fwrite(specAlphaLong, file = paste0(plotDir, "Species_alpha_diversity_for_plotting.csv"))

#----- Pivot longer for plotting (genus)
genusAlphaLong <- genusAlpha |>
  pivot_longer(cols = c(Shannon, Simpson, Observed),
               names_to = "Metric",
               values_to = "Value") |>
    as.data.frame()
fwrite(genusAlphaLong, file = paste0(plotDir, "Genus_alpha_diversity_for_plotting.csv"))

#----- Plot species shannon
x <- ggplot(specAlpha, aes(x = Sample, y = Shannon)) +
  geom_col(color = "black", width = 0.7, alpha = 0.8) +
  geom_text(aes(label = round(Shannon, 2)), vjust = -0.4, size = 3.5) +
  theme_classic(base_size = 15) +
  labs(
    title = "Shannon Alpha Diversity per Sample",
    x = "",
    y = "Shannon Index"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )
ggsave(paste0(plotDir, "Species_Shannon_Diversity_Per_Sample.png"), x)

#----- Plot Genus Shannon
y <- ggplot(genusAlpha, aes(x = Sample, y = Shannon)) +
    geom_col(color = "black", width = 0.7, alpha = 0.8) +
    geom_text(aes(label = round(Shannon, 2)), vjust = -0.4, size = 3.5) +
    theme_classic(base_size = 15) +
    labs(
        title = "Shannon Alpha Diversity per Sample",
        x = "",
        y = "Shannon Index"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank()
    )
ggsave(paste0(plotDir, "Genus_Shannon_Diversity_Per_Sample.png"), y)
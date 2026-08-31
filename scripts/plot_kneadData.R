#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: plot_kneadData.R
# Description: Generate diagnostic figure that shows how
# many reads were trimmed with KneadData
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

#----- Libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(purrr))

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
# The directories are passed in rather than assumed, because the pipeline
# writes under out_dir from config.yaml, which is not the directory snakemake
# runs from. Hardcoding "filtered_fq" here silently found no logs at all.
if (length(args) != 2) {
  stop("Usage: Rscript plot_kneadData.R <filtered_dir> <figures_dir>")
}

#----- Set variables based on command line args
filteredDir <- args[1]
figDir      <- args[2]

#----- Make sure the paths end in a separator, since they are pasted onto
# filenames below
addSlash <- function(p) if (grepl("/$", p)) p else paste0(p, "/")
filteredDir <- addSlash(filteredDir)
figDir      <- addSlash(figDir)

#----- Create figure directory
if (!dir.exists(figDir)) dir.create(figDir, recursive = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA AND PARSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Find all KneadData log files recursively
# kneaddata leaves {sample}.log inside filtered_fq/{sample}/, so this needs
# to recurse
logFileVec <- list.files(filteredDir, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)

#----- Stop if there are none.
# Without this the script carries on with an empty table and writes a blank
# figure, which looks like a successful run.
if (length(logFileVec) == 0) {
  stop(paste0("No kneaddata .log files found under ", filteredDir,
              ". Nothing to plot."))
}

message(paste0("Found ", length(logFileVec), " kneaddata log file(s)"))

#----- Function to parse one KneadData log
parseKD <- function(x) {
  text <- readLines(x)

  #== Original reads
  raw1 <- as.numeric(str_extract(str_subset(text, "READ COUNT: raw pair1"), "\\d+\\.\\d+"))
  raw2 <- as.numeric(str_extract(str_subset(text, "READ COUNT: raw pair2"), "\\d+\\.\\d+"))

  #== After trimming
  trim1 <- as.numeric(str_extract(str_subset(text, "READ COUNT: trimmed pair1"), "\\d+\\.\\d+"))
  trim2 <- as.numeric(str_extract(str_subset(text, "READ COUNT: trimmed pair2"), "\\d+\\.\\d+"))

  #== After decontamination
  decon1 <- as.numeric(str_extract(str_subset(text, "READ COUNT: decontaminated .* orphan1"), "\\d+\\.\\d+"))
  decon2 <- as.numeric(str_extract(str_subset(text, "READ COUNT: decontaminated .* orphan2"), "\\d+\\.\\d+"))

  #== Collate as data.frame
  data.frame(
    sample = sub("\\.log$", "", basename(x)),
    raw_reads = sum(raw1, raw2, na.rm = TRUE),
    trimmed_reads = sum(trim1, trim2, na.rm = TRUE),
    decon_reads = sum(decon1, decon2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#----- Apply the function to all logs
kdSummary <- map_dfr(logFileVec, parseKD)

#----- Optional: add summary percentages
kdSummary <- kdSummary %>%
  mutate(
    trim_percent = 100 * trimmed_reads / raw_reads,
    decon_percent = 100 * decon_reads / raw_reads
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PIVOT DATA FOR PLOTTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Pivot to long format
kd_long <- kdSummary %>%
  select(sample, raw_reads, trimmed_reads, decon_reads) %>%
  pivot_longer(
    cols = c(raw_reads, trimmed_reads, decon_reads),
    names_to = "stage",
    values_to = "reads"
  ) %>%
  mutate(
    stage = factor(stage,
                   levels = c("raw_reads", "trimmed_reads", "decon_reads"),
                   labels = c("Raw", "Trimmed", "Decontaminated"))
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PLOT DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

kdPlot <- ggplot(kd_long, aes(x = sample, y = reads, fill = stage)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Sample",
    y = "Number of Reads",
    fill = "Read Type",
    title = "") +
  scale_fill_manual(values = c("Raw" = "#333333", "Trimmed" = "#56B4E9", "Decontaminated" = "#E69F00")) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#----- Save plot
ggsave(paste0(figDir, "KneadData_Diagnostics.png"), kdPlot, width = 12, height = 10)
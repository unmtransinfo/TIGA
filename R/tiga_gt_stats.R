#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT STATS
### tiga_gt_stats.R - Produce gt_stats.tsv, for TIGA Shiny app.
### Input from tiga_gt_variables.R.
#############################################################################
# Now using mean-rank or median-rank instead of mu_scores.
# Based on DISEASES benchmark results, using only variables:
#   * n_study
#   * pvalue_mlog_median (MAYBE CHANGING TO pvalue_mlog_max)
#   * rcras
#############################################################################
### Previously we computed meanRank for (1) genes for a given trait, and
### (2) traits for a given gene. Another way providing additional 
### functionality is to compute meanRank for gene-trait pairs (GTs) relative
### to all GTs. Since the input variables all correspond to GTs.
#############################################################################
### By expressing ranks as [1,100] percentiles rank_ptl, 100/rank_ptl is
### normalized to [1,100] weighted more for higher rank intervals.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
#
if (interactive()) {
  rel_y <- as.integer(readline(prompt="Enter RELEASE_YEAR: "))
  rel_m <- as.integer(readline(prompt="Enter RELEASE_MONTH: "))
  rel_d <- as.integer(readline(prompt="Enter RELEASE_DAY: "))
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (length(args)==3) {
  rel_y <- as.integer(args[1])
  rel_m <- as.integer(args[2])
  rel_d <- as.integer(args[3])
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (file.exists("LATEST_RELEASE_GWC.txt")) {
  GC_REL <- trimws(read_file("LATEST_RELEASE_GWC.txt"))
  rel_y <- as.integer(sub("\\-.*$", "", GC_REL))
  rel_m <- as.integer(sub("\\d+\\-(\\d+)\\-.*$", "\\1", GC_REL))
  rel_d <- as.integer(sub("\\d+\\-\\d+\\-(\\d+).*$", "\\1", GC_REL))
  message(sprintf("LATEST_RELEASE_GWC: %s", GC_REL))
  ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))
} else {
  message("ERROR: Syntax:  tiga_gt_stats.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
#
ifile	<- paste0(ODIR, "/gt_variables.tsv.gz")
ofile	<- paste0(ODIR, "/gt_stats.tsv.gz")
#
writeLines(sprintf("Input file: %s", ifile))
writeLines(sprintf("Output file: %s", ofile))
#
###
# Read computed variables from file.
gt_stats <- read_delim(ifile, "\t", col_types = cols(.default = col_character(),
n_study = col_integer(),
n_snp = col_integer(),
n_snpw = col_double(),
geneNtrait = col_integer(),
geneNstudy = col_integer(),
traitNgene = col_integer(),
traitNstudy = col_integer(),
pvalue_mlog_median = col_double(),
pvalue_mlog_max = col_double(),
or_median = col_double(),
n_beta = col_integer(),
study_N_mean = col_double(),
rcras = col_double(),
geneIdgList = col_logical()))
setDT(gt_stats)
#
gt_stats <- gt_stats[!is.na(efoId)] #Should be in tiga_gt_prepfilter.R
gt_stats <- gt_stats[!is.na(pvalue_mlog_max)] #Should be in tiga_gt_prepfilter.R
#
###
# Mean-rank computation. For each variable in defined set, 
# compute rank with ties having same rank, then compute mean of 
# variable-ranks.
# base::rank cannot equate NAs so need custom my_rank().
# Best rank is 1. All variables should be bigger is better.
###
my_rank <- function(v) {
  ranks <- base::rank(v, na.last="keep", ties.method="min")
  rank_max <- max(ranks, na.rm=T)
  ranks <- ifelse(is.na(ranks), rank_max+1, ranks)
  return(ranks)
}
#
###
TAGS_FOR_RANKING <- c("pvalue_mlog_max", "rcras", "n_snpw")
###
ranks_this <- list()
for (tag in TAGS_FOR_RANKING) {
  ranks_this[[tag]] <- my_rank(-gt_stats[[tag]])
}
setDT(ranks_this)
ranks_this[, meanRank := rowMeans(.SD)]
gt_stats[, meanRank := ranks_this$meanRank]
#
###
# Normalize to meanRankPtl (percentile).
# ecdf() is empirical cumulative distribution function
ranks2pctiles <- function(ranks) {
  pctiles <- 100 * (1.0 - ecdf(ranks)(ranks))
  return(pctiles)
}
#
gt_stats[, meanRankScore := ranks2pctiles(meanRank)]
#
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
message("Testing and checking behavior of meanRankScore vs meanRank:")
message(sprintf("meanRank values: %d; unique values: %d", nrow(gt_stats[!is.na(meanRank)]), gt_stats[, uniqueN(meanRank)]))
message(sprintf("meanRankScore values: %d; unique values: %d", nrow(gt_stats[!is.na(meanRankScore)]), gt_stats[, uniqueN(meanRankScore)]))
gt_stats <- gt_stats[order(meanRank)]
message(sprintf("meanRank monotonically increasing: %s", all(gt_stats$meanRank == cummax(gt_stats$meanRank))))
message(sprintf("meanRankScore monotonically decreasing: %s", all(gt_stats$meanRankScore == cummin(gt_stats$meanRankScore))))
frequent_meanRank <- gt_stats[, .(N = as.integer(.N)), by=meanRank][order(-N)]
writeLines(sprintf("Frequent meanRank (tied values) #%2d: %12.2f (N=%4d)", 1:10, frequent_meanRank[1:10, meanRank], frequent_meanRank[1:10, N]))
message(sprintf("Total tied values: %d / %d (%.1f%%)", frequent_meanRank[N>1, sum(N)], frequent_meanRank[, sum(N)], 100*frequent_meanRank[N>1, sum(N)]/frequent_meanRank[, sum(N)]))
ecdf_meanRank <- ecdf(gt_stats$meanRank)
print(summary(ecdf_meanRank))
ecdf_meanRankScore <- ecdf(gt_stats$meanRankScore)
print(summary(ecdf_meanRankScore))
if (interactive()) {
  plot(ecdf_meanRank)
  plot(ecdf_meanRankScore)
}
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT STATS
### tiga_gt_stats.R - Produce gt_stats.tsv, for TIGA Shiny app.
### Input from tiga_gt_variables.R.
#############################################################################
# Now using mean-rank or median-rank instead of mu_scores.
# Based on DISEASES benchmark results, using only variables:
#   * n_study
#   * pvalue_mlog_median
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
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==2) {
  (ifile	<- args[1])
  (ofile	<- args[2])
} else if (length(args)==0) {
  ifile <- "data/gt_variables.tsv.gz"
  ofile <- "data/gt_stats.tsv.gz"
} else {
  message("ERROR: Syntax: tiga_gt_stats.R VARIABLESFILE OFILE\n...or... no args for defaults")
  quit()
}
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
or_median = col_double(),
n_beta = col_integer(),
study_N_mean = col_double(),
rcras = col_double(),
geneIdgList = col_logical()))
setDT(gt_stats)
#
gt_stats <- gt_stats[!is.na(efoId)] #Should be in tiga_gt_prepfilter.R
gt_stats <- gt_stats[!is.na(pvalue_mlog_median)] #Should be in tiga_gt_prepfilter.R
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
TAGS_FOR_RANKING <- c("pvalue_mlog_median", "rcras", "n_snpw")
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
# Normalize to meanRankPtl (percentile). Inverted so high percentile is high (numerically low) rank.
# PROBLEM: This is only approximate percentile.
#ranks2pctiles <- function(rs) {
#  q <- quantile(rs, seq(0, .999, .001)) #1000 bins.
#  for (i in 1:length(rs)) {
#    rs[i] <- 100 - max(which(q <= rs[i])) / 10
#  }
#  return(rs)
#}
#
# ecdf() is empirical cumulative distribution function
ranks2pctiles <- function(rs) {
  pctiles <- 100 * (1.0 - ecdf(rs)(rs))
  return(pctiles)
}
#
gt_stats[, meanRankScore := ranks2pctiles(meanRank)]
#
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
message(sprintf("cor(meanRank, meanRankScore): %f", cor(gt_stats$meanRank, gt_stats$meanRankScore)))
gt_stats <- gt_stats[order(meanRank)]
message(sprintf("meanRank monotonically increasing: %s", all(gt_stats$meanRank == cummax(gt_stats$meanRank))))
message(sprintf("meanRankScore monotonically decreasing: %s", all(gt_stats$meanRankScore == cummin(gt_stats$meanRankScore))))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT stats
### tiga_gt_stats.R - Produce gt_stats.csv, for TIGA Shiny app.
### Input from tiga_gt_variables.R.
### ~30min
#############################################################################
# Now using mean-rank or median-rank instead of mu_scores.
# Based on DISEASES benchmark results, using only variables:
#   * n_study
#   * pvalue_mlog_median
#   * rcras
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==5) {
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
###
TAGS_FOR_RANKING <- c("n_study", "pvalue_mlog_median", "rcras")
###
# Gene meanrank:
gt_stats[, geneMeanRank := as.numeric(NA)]
ii <- 0
for (efoId_this in unique(gt_stats$efoId)) {
  ii <- ii + 1
  trait_this <- gt_stats[efoId==efoId_this, trait][1]
  n_gene_this <- gt_stats[efoId==efoId_this, uniqueN(ensemblId)]
  message(sprintf("[%d / %d] (N_gene: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$efoId), n_gene_this, efoId_this, trait_this))
  if (n_gene_this<2) { #singletons
    gt_stats[efoId==efoId_this]$geneMeanRank <- 1
    next;
  }
  ranks_this <- list()
  for (tag in TAGS_FOR_RANKING) {
    ranks_this[[tag]] <- my_rank(-gt_stats[efoId==efoId_this][[tag]])
  }
  setDT(ranks_this)
  ranks_this[, meanRank := rowMeans(.SD)]
  gt_stats[efoId==efoId_this]$geneMeanRank <- ranks_this$meanRank
}
#
###
# Trait meanrank:
gt_stats[, traitMeanRank := as.numeric(NA)]
ii <- 0
for (ensemblId_this in unique(gt_stats$ensemblId)) {
  n_trait_this <- gt_stats[ensemblId==ensemblId_this, uniqueN(efoId)]
  ii <- ii + 1
  geneSymbol_this <- gt_stats[ensemblId==ensemblId_this, geneSymbol][1]
  message(sprintf("[%d / %d] (N_trait: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$ensemblId), n_trait_this, ensemblId_this, geneSymbol_this))
  if (n_trait_this<2) { #singletons
    gt_stats[ensemblId==ensemblId_this]$traitMeanRank <- 1
    next;
  }
  ranks_this <- list()
  for (tag in TAGS_FOR_RANKING) {
    ranks_this[[tag]] <- my_rank(-gt_stats[ensemblId==ensemblId_this][[tag]])
  }
  setDT(ranks_this)
  ranks_this[, meanRank := rowMeans(.SD)]
  gt_stats[ensemblId==ensemblId_this]$traitMeanRank <- ranks_this$meanRank
}
###
gt_stats[, `:=`(geneMeanRankScore = 100/geneMeanRank, traitMeanRankScore = 100/traitMeanRank)]
#
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2f %s", t_elapsed, attr(t_elapsed, "units")))

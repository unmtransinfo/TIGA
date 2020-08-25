#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT stats
### tiga_gt_stats_mu.R - Produce gt_stats.csv, for TIGA Shiny app.
#############################################################################
# Multivariable non-parametric ranking via &mu; scores.
###
# Non-dominated solutions are not inferior to any other case at any variable.
# A &mu; score is defined as the number of lower cases minus the number of higher.
# The resulting ranking is the useful result, not so much the score itself.
# Using muStat package.
###
# mu.GE: greater than or equal to
# mu.AND: Logical AND GEs for all variables
###
# Issue: cases globally superior in one variable and inferior in one variable
# have nAbove=0 and nBelow=0 and muScore=0. What should be the rank?
# If a case has nAbove=0 and nBelow=0 then weight=0 and mu.Sums() returns score = NA.
# This is arbitrary so we assign score = nBelow = nAbove = 0 - 0 = 0.
# ?mu.Sums:
# score  = (nB-nA) * ifelse(weight==0, NA, 1)
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
library(muStat)
#
t_start <- Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==2) {
  (ifile        <- args[1])
  (ofile        <- args[2])
} else if (length(args)==0) {
  ifile <- "data/gt_variables.tsv.gz"
  ofile <- "data/gt_stats_mu.tsv.gz"
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
###
# &mu; scores/rankings
#   * __for a given trait__
#   * __for a given gene__
###
#
# Use inverse of n\_traits\_g and n\_genes\_t so bigger is better as needed for muStat.
gt_stats[, geneNtrait_inv := 1 / geneNtrait]
gt_stats[, traitNgene_inv := 1 / traitNgene]

# For each (trait|gene), convert to matrix for muStat::mu.GE().
# The (i,j) entry of GE matrix is 1 if \code{x_i >= x_j}, 0 otherwise.
# The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.
###
# Gene &mu; scores:
# Some or_median will be NA, since beta included.
gt_stats[, `:=`(geneMuScore=as.integer(NA), geneMuRank=as.integer(NA))]
ii <- 0
for (efoId_this in unique(gt_stats$efoId)) {
  gtmat <- as.matrix(gt_stats[efoId==efoId_this, .(n_study, n_snp, n_snpw, geneNtrait_inv, traitNgene_inv, pvalue_mlog_median, or_median, n_beta, rcras)])
  ii <- ii + 1
  trait_this <- gt_stats[efoId==efoId_this, trait][1]
  message(sprintf("[%d / %d] (N_gene: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$efoId), dim(gtmat)[1], efoId_this, trait_this))
  if (dim(gtmat)[1]<2) { #Singletons
    gt_stats[efoId==efoId_this]$geneMuScore <- 0
    gt_stats[efoId==efoId_this]$geneMuRank <- 1
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt_stats[efoId==efoId_this, .(geneName)]
  sums[is.na(score) & weight==0, score := 0] # Override NA scores
  sums[order(-score, -weight), rank := 1:nrow(sums)]
  gt_stats[efoId==efoId_this]$geneMuScore <- sums$score
  gt_stats[efoId==efoId_this]$geneMuRank <- sums$rank
}
#
###
# Trait &mu; scores:
gt_stats[, `:=`(traitMuScore=as.integer(NA), traitMuRank=as.integer(NA))]
ii <- 0
for (ensemblId_this in unique(gt_stats$ensemblId)) {
  gtmat <- as.matrix(gt_stats[ensemblId==ensemblId_this, .(n_study, n_snp, n_snpw, geneNtrait_inv, traitNgene_inv, pvalue_mlog_median, or_median, n_beta, rcras)])
  ii <- ii + 1
  geneSymbol_this <- gt_stats[ensemblId==ensemblId_this, geneSymbol][1]
  message(sprintf("[%d / %d] (N_trait: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$ensemblId), dim(gtmat)[1], ensemblId_this, geneSymbol_this))
  if (dim(gtmat)[1]<2) { #Singletons
    gt_stats[ensemblId==ensemblId_this]$traitMuScore <- 0
    gt_stats[ensemblId==ensemblId_this]$traitMuRank <- 1
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt_stats[ensemblId==ensemblId_this, .(geneName)]
  sums[is.na(score) & weight==0, score := 0] # Override NA scores
  sums[order(-score, -weight), rank := 1:nrow(sums)]
  gt_stats[ensemblId==ensemblId_this]$traitMuScore <- sums$score
  gt_stats[ensemblId==ensemblId_this]$traitMuRank <- sums$rank
}
###
gt_stats[, `:=`(geneNtrait_inv = NULL, traitNgene_inv = NULL)]
#
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

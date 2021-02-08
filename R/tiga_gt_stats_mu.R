#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT STATS (MU VERSION)
### tiga_gt_stats_mu.R - Produce mu scores.
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
## Problem: muStat requires huge memory (>100Gb?)
## Error: vector memory exhausted (limit reached?)
## Error: negative length vectors are not allowed (join > 2^31-1 rows)
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
library(muStat)
#
t_start <- Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
#ifile   <- ifelse(length(args)>0, args[1], "data/gt_variables.tsv.gz")
#ofile   <- ifelse(length(args)>1, args[2], "data/gt_stats_mu.tsv.gz")
ifile   <- ifelse(length(args)>0, args[1], "data/20201216/gt_variables.tsv.gz")
ofile   <- ifelse(length(args)>1, args[2], "data/20201216/gt_stats_mu.tsv.gz")
#
if (length(args)>2) {
  message("ERROR: Syntax: tiga_gt_stats.R VARIABLESFILE OFILE\n...or... no args for defaults")
  quit()
}
message(sprintf("Input file: %s", ifile))
message(sprintf("Output file: %s", ofile))
#
###
# Read computed variables from file.
gt_stats <- read_delim(ifile, "\t", col_types=cols(.default=col_character(),
	n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
	geneNtrait=col_integer(), geneNstudy=col_integer(),
	traitNgene=col_integer(), traitNstudy=col_integer(),
	pvalue_mlog_median=col_double(),
	or_median=col_double(), n_beta=col_integer(),
	rcras=col_double(),
	study_N_mean=col_double(), geneIdgList=col_logical()))
setDT(gt_stats)
#
###
# &mu; scores/rankings __for a given gene-trait__
# (Now computing globally, for all gene-trait pairs.)
###
#
# Use inverse of n\_traits\_g and n\_genes\_t so bigger is better as needed for muStat.
#gt_stats[, geneNtrait_inv := 1 / geneNtrait]
#gt_stats[, traitNgene_inv := 1 / traitNgene]
#
TAGS_FOR_RANKING <- c("pvalue_mlog_max", "rcras", "n_snpw")
message(paste("TAGS_FOR_RANKING: ", paste0(TAGS_FOR_RANKING, collapse=", ")))
#
# Convert to matrix for muStat::mu.GE().
# The (i,j) entry of GE matrix is 1 if \code{x_i >= x_j}, 0 otherwise.
# The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.
###
# Some or_median will be NA, since beta included.
gt_stats[, `:=`(muScore=as.integer(NA), muRank=as.integer(NA))]
message(sprintf("nrow(gt_stats): %d", nrow(gt_stats)))
gtmat <- as.matrix(gt_stats[, ..TAGS_FOR_RANKING])
#gt_stats[, `:=`(geneNtrait_inv = NULL, traitNgene_inv = NULL)]
###
message("Running mu.GE (memory hog, e.g. 50GB for 50k rows x 3vars)...")
ge <- mu.GE(gtmat)
###
sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
setDT(sums)
sums$name <- gt_stats[, .(geneName)]
sums[is.na(score) & weight==0, score := 0] # Override NA scores
sums[order(-score, -weight), rank := 1:nrow(sums)]
#
gt_stats[, muScore := sums$score]
gt_stats[, muRank := sums$rank]
#
write_delim(gt_stats, ofile, delim="\t")
message(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

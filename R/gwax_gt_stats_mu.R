#!/usr/bin/env Rscript
###
# Run after gwax_gt_stats.R (could be combined), with input from that code,
# gt_stats.tsv.gz.
###
# Multivariable non-parametric ranking via &mu; scores.
# &mu; scores
###
# Non-dominated solutions are not inferior to any other case at any variable.
# A &mu; score is defined as the number of lower cases minus the number of higher.
# The resulting ranking is the useful result, not so much the score itself.
# Using muStat package.
###
# Issue: cases globally superior in one variable and inferior in one variable
# have nAbove=0 and nBelow=0 and muScore=0. What should be the rank?
###
library(readr, quietly=T)
library(data.table, quietly=T)
library(muStat)

t0 <- proc.time()

# Read gene-trait file.
gt <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_wsnp=col_double(), n_traits_g=col_integer(), n_genes_t=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double()))
setDT(gt)
#
gt <- gt[!is.na(or_median) & !is.na(gsymb) &!is.na(tdl)]
###
#traits <- gt[, .(N_gene = .N),  by=c("trait", "trait_uri")]
#MIN_ASSN <- 20
#traits <- traits[N_gene>=MIN_ASSN]
#gt <- gt[trait_uri %in% traits$trait_uri]
###
#gt$name <- paste(gt$gsymb, sub("^.*/", "", gt$trait_uri), sep=":")
#
message(sprintf("Genes in dataset: %d", uniqueN(gt$gsymb)))
message(sprintf("Traits in dataset: %d", uniqueN(gt$trait_uri)))
message(sprintf("G-T associations in dataset: %d", nrow(gt)))

# Use inverse of n\_traits\_g and n\_genes\_t so bigger is better as needed for muStat.
gt[, n_traits_g_inv := 1 / n_traits_g]
gt[, n_genes_t_inv := 1 / n_genes_t]

# We are interested in rankings __for a given trait__. So for each trait,
# convert to matrix for muStat::mu.GE().
# The (i,j) entry of GE matrix is 1 if \code{x_i >= x_j}, 0 otherwise. The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.

gt[, `:=`(mu_score=as.integer(NA), nAbove=as.integer(NA), nBelow=as.integer(NA), mu_rank=as.integer(NA))]

ii <- 0
for (trait_this in unique(gt$trait)) {
  gtmat <- as.matrix(gt[trait==trait_this, .(n_study, n_snp, n_wsnp, n_traits_g_inv, n_genes_t_inv, pvalue_mlog_median, or_median, rcras)])
  ii <- ii + 1
  message(sprintf("[%d / %d] (N_gene: %3d) \"%s\"", ii, uniqueN(gt$trait), dim(gtmat)[1], trait_this))
  if (dim(gtmat)[1]<2) { #No ranking for singleton.
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt[trait==trait_this, .(name)]
  # Bug in muStat?
  if (sum(is.na(sums$score))>0) {
    message(sprintf("DEBUG: missing score count: %d", sum(is.na(sums$score))))
    badrows <- is.na(sums$score)
    print(sums[badrows]) #DEBUG
    sums[is.na(score), score := nAbove - nBelow]
    message(sprintf("DEBUG: missing score count: %d (FIXED?)", sum(is.na(sums$score))))
    print(sums[badrows]) #DEBUG
  }
  sums <- setorder(sums, -score, na.last=T)
  sums[, rank := 1:nrow(sums)]
  gt[trait==trait_this]$mu_score <- sums$score
  gt[trait==trait_this]$mu_rank <- sums$rank
  gt[trait==trait_this]$nAbove <- sums$nAbove
  gt[trait==trait_this]$nBelow <- sums$nBelow
}
#
gt[, `:=`(n_traits_g_inv = NULL, n_genes_t_inv = NULL)]
write_delim(gt, "data/gt_stats_mu.tsv", delim="\t")
#
message(sprintf("NOTE: elapsed time: %.2fs",(proc.time()-t0)[3]))

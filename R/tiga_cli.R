#!/usr/bin/env Rscript
#############################################################################
### CLI for GWAX
### gt_stats.tsv.gz produced by gwax_gt_stats.R.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)

###
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  (efoId   <- args[1])
  (ofile   <- "")
} else if (length(args)==2) {
  (efoId   <- args[1])
  (ofile   <- args[2])
} else {
  message("ERROR: Syntax: gwax_cli.R EFOID [OFILE]")
  quit()
}
#
t0 <- proc.time()

#
gt <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), 
n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
geneNtrait=col_integer(), geneNstudy=col_integer(),
traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(),
geneMuScore=col_double(), geneMuRank=col_integer(),
traitMuScore=col_double(), traitMuRank=col_integer()))
setDT(gt)
setnames(gt, old=c("geneIdgTdl"), new=c("TDL"))
#
message(sprintf("Gene count, IDs: %d; symbols: %d", uniqueN(gt$ensemblId), uniqueN(gt$geneSymbol)))
message(sprintf("Trait count (total): %d", uniqueN(gt$efoId)))

###
Hits <- function(qry) {
  gt_this <- gt[efoId==qry]
  if (nrow(gt_this)==0) { return(NULL) }
  gt_this$TDL <- factor(gt_this$TDL, levels=c("NA", "Tdark", "Tbio", "Tchem", "Tclin"), ordered=T)
  gt_this <- gt_this[, .(ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, geneNtrait, pvalue_mlog_median, or_median, rcras, geneMuScore, geneMuRank)]
  setnames(gt_this, old=c("geneMuScore", "geneMuRank"), new=c("muScore", "muRank"))
  setorder(gt_this, muRank)
  return(gt_this)
}

gt_this <- Hits(efoId)
message(sprintf("Hits: %d", nrow(gt_this)))
fwrite(gt_this, ofile, sep="\t")
t_elapsed <- (proc.time()-t0)[3]
message(sprintf("Elapsed: %.1fs", t_elapsed))

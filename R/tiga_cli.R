#!/usr/bin/env Rscript
#############################################################################
### CLI for TIGA
### gt_stats.tsv.gz produced by tiga_gt_stats.R.
#############################################################################
# EFO_0001360 (T2DM)
# EFO_0000249 (Alzheimer)
# EFO_0000289 (bipolar disorder)
# EFO_0000249 (Alzheimers disease)
# EFO_0000305 (breast carcinoma)
# EFO_0000270 (asthma)
# EFO_0000692 (schizophrenia)
# EFO_0005842 (colorectal cancer)
# EFO_0001663 (prostate carcinoma)
# EFO_0000685 (rheumatoid arthritis)
# EFO_0003761 (unipolar depression)
# EFO_0002508 (Parkinsons disease)
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)

###
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  (efoId   <- args[1])
  (ofile   <- "")
} else if (length(args)==2) {
  (efoId   <- args[1])
  (ofile   <- args[2])
} else {
  message("ERROR: Syntax: tiga_cli.R EFOID [OFILE]\n(e.g. EFO_0000270)")
  quit()
}
#
gt <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(), geneNtrait=col_integer(), geneNstudy=col_integer(), traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), n_beta=col_double(), study_N_mean=col_double(), rcras=col_double(), meanRank=col_double(), meanRankScore=col_double()))
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
  gt_this <- gt_this[, .(ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, pvalue_mlog_median, rcras, meanRankScore, study_N_mean, n_snp, n_snpw, geneNtrait, or_median, n_beta)]
  setorder(gt_this, -meanRankScore)
  return(gt_this)
}

gt_this <- Hits(efoId)
message(sprintf("Hits: %d", nrow(gt_this)))
fwrite(gt_this, ofile, sep="\t")

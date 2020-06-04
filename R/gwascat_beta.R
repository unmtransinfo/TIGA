#!/usr/bin/env Rscript
#############################################################################
### Input file from gwascat_assn.R,
### which separated OR from BETA, via heuristic: If all OR_or_BETA values
### for a study are >=1, assume OR.
#############################################################################
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- "data/gwascat_assn.tsv"
  ofile <- "data/gwascat_beta.tsv"
} else {
  message("ERROR: Syntax: gwascat_assn.R ASSNFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), oddsratio=col_double(), beta=col_double()))
setDT(assn)

writeLines(sprintf("Studies with OR: %d", assn[!is.na(oddsratio), uniqueN(STUDY_ACCESSION)]))
writeLines(sprintf("Studies with BETA: %d", assn[!is.na(beta), uniqueN(STUDY_ACCESSION)]))

###
qs <- quantile(assn[["beta"]][!is.na(assn[["beta"]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("beta %4s-ile: %9.1f", names(qs), qs))

qs <- quantile(assn[["oddsratio"]][!is.na(assn[["oddsratio"]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("oddsratio %4s-ile: %9.1f", names(qs), qs))

###
# Heuristic: all units include (in|de)crease
assn[, beta_units := sub("^.*\\] *", "", `95%_CI_(TEXT)`)]
assn$beta_units[!grepl("(in|de)crease", assn$beta_units)] <- NA
tbl <- sort(table(assn$beta_units), decreasing = T)
writeLines(sprintf("%d. (N=%d) %s", 1:100, tbl[1:100], names(tbl)[1:100]))

#Top BETA units
tbl <- as.data.frame(table(assn$beta_units))
colnames(tbl) <- c("beta_units", "Freq")
tbl <- tbl[order(-tbl$Freq),]
writeLines(sprintf("%18s: %3d", tbl$beta_units[1:20], tbl$Freq[1:20]))


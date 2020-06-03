#!/usr/bin/env Rscript
#############################################################################
### Clean studies file.
#############################################################################
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- paste0(Sys.getenv("HOME"), "/../data/gwascatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv")
  ofile <- "data/gwascat_gwas.tsv"
} else {
  message("ERROR: Syntax: gwascat_gwas.R GWASFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

gwas <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), DATE=col_date(), `ASSOCIATION COUNT`=col_integer(), `DATE ADDED TO CATALOG`=col_date()))
setDT(gwas)

setnames(gwas, gsub("[ \\./]", "_", colnames(gwas)))
setnames(gwas, gsub("__", "_", colnames(gwas)))
setnames(gwas, gsub("_$", "", colnames(gwas)))

#Convert special chars.
for (tag in colnames(gwas)) {
  if (typeof(gwas[[tag]])=="character") {
    writeLines(sprintf("NOTE: cleaning: %s", tag))
    gwas[[tag]] <- iconv(gwas[[tag]], from="latin1", to="UTF-8")
  }
}

# Parse study_N from INITIAL_SAMPLE_SIZE. Sum sub-samples.
gwas[, study_N := INITIAL_SAMPLE_SIZE] #data.table warning spurious. 
gwas[, study_N := gsub("PMID:[0-9]+", "", study_N)]
gwas[, study_N := gsub("[^0-9,]+", "\t", study_N)]
gwas[, study_N := gsub(",", "", study_N)]
gwas[, study_N := sub("^\t", "", study_N)]
gwas[, study_N := gsub("\t\t*", "\t", study_N)]
gwas[, study_N := sub("\t$", "", study_N)]
gwas[, study_N := gsub("\t", "+", study_N)]
gwas[, study_N := sapply(sapply(study_N, str2lang), eval)]
gwas[, study_N := as.integer(study_N)]
writeLines(sprintf("UNPARSED INITIAL_SAMPLE_SIZE: %s", gwas[is.na(study_N), INITIAL_SAMPLE_SIZE]))

message(sprintf("Total gwas count: %6d", nrow(gwas)))

write_delim(gwas, ofile, delim="\t")
###

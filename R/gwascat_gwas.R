#!/usr/bin/env Rscript
#############################################################################
### Clean studies file.
#############################################################################
library(readr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- "/home/data/gwascatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv"
  ofile <- "data/gwascat_gwas.tsv"
} else {
  message("ERROR: Syntax: gwascat_gwas.R GWASFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

gwas <- read_delim(ifile, "\t")

colnames(gwas) <- gsub("[ \\./]","_",colnames(gwas))
colnames(gwas) <- gsub("__","_",colnames(gwas))
colnames(gwas) <- gsub("_$","",colnames(gwas))

gwas <- gwas[complete.cases(gwas[,c("STUDY_ACCESSION","PUBMEDID","DISEASE_TRAIT")]),]

#Convert special chars.
for (tag in colnames(gwas)) {
  if (typeof(gwas[[tag]])=="character") {
    writeLines(sprintf("NOTE: cleaning: %s", tag))
    gwas[[tag]] <- iconv(gwas[[tag]], from="latin1", to="UTF-8")
  }
}

writeLines(sprintf("Total gwas count: %6d", nrow(gwas)))

write_delim(gwas, ofile, delim="\t")
###

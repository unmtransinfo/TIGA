#!/usr/bin/env Rscript
#############################################################################
### Clean studies file.
#############################################################################
library(readr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0)
{
  (ifile <- args[1])
} else {
  ifile <- "/home/data/gwascatalog/data/gwas_catalog_v1.0.1-studies_r2017-10-10.tsv"
}
if (length(args)>1)
{
  (ofile <- args[2])
} else {
  ofile <- "data/gwascat_gwas.tsv"
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

gwas <- read_delim(ifile, "\t")

colnames(gwas) <- gsub("[ \\./]","_",colnames(gwas))
colnames(gwas) <- gsub("__","_",colnames(gwas))
colnames(gwas) <- gsub("_$","",colnames(gwas))

gwas <- gwas[complete.cases(gwas[,c("STUDY_ACCESSION","PUBMEDID","DISEASE_TRAIT")]),]

#Convert special chars.
gwas$DISEASE_TRAIT <- iconv(gwas$DISEASE_TRAIT, from="latin1", to="UTF-8")

writeLines(sprintf("Total gwas count: %6d", nrow(gwas)))

write_delim(gwas, ofile, delim="\t")
###

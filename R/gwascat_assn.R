#!/usr/bin/env Rscript
#############################################################################
### Separate OR from beta, create new columns and write TSV.
### Heuristic: If all OR_or_beta values for a study are >=1, assume OR. (Not necessarily true as beta may be >1.) 
###
### Aha.  Column `95%_CI_(TEXT)` has units for betas.  But must be parsed.
### From http://www.ebi.ac.uk/gwas/docs/fileheaders:
### "95% CI (TEXT)*: Reported 95% confidence interval associated with strongest SNP risk allele, 
### along with unit in the case of beta-coefficients. If 95% CIs are not published, we estimate 
### these using the standard error, where available."
#############################################################################
library(readr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0)
{
  (ifile <- args[1])
} else {
  ifile <- "/home/data/gwascatalog/data/gwas_catalog_v1.0.1-associations_e90_r2017-10-10.tsv"
}
if (length(args)>1)
{
  (ofile <- args[2])
} else {
  ofile <- "data/gwascat_assn.tsv"
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

assn <- read_delim(ifile, "\t")

colnames(assn) <- gsub("[ \\./]","_",colnames(assn))
colnames(assn) <- gsub("__","_",colnames(assn))
colnames(assn) <- gsub("_$","",colnames(assn))

assn <- assn[complete.cases(assn[,c("STUDY_ACCESSION","SNPS","DISEASE_TRAIT")]),]

assn$DISEASE_TRAIT <- iconv(assn$DISEASE_TRAIT, from="latin1", to="UTF-8")

writeLines(sprintf("Total assn count: %6d", nrow(assn)))
writeLines(sprintf("OR_or_beta MISSING: %6d", nrow(assn[is.na(assn$OR_or_BETA),])))
writeLines(sprintf("OR_or_beta values: %6d", nrow(assn[!is.na(assn$OR_or_BETA),])))

tag="OR_or_BETA"
assn$OR_or_BETA <- as.numeric(assn$OR_or_BETA)
qs <- quantile(assn[[tag]][!is.na(assn[[tag]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("%s %4s-ile: %9.1f",tag,names(qs),qs))

assn$oddsratio <- NA
assn$beta <- NA

accs <- sort(unique(assn$STUDY_ACCESSION))
writeLines(sprintf("Studies, total: %6d", length(accs)))
i <- 0
n_all_or <- 0
n_all_beta <- 0
n_both <- 0
for (acc in accs) {
  i <- i + 1
  or_or_beta_this <- assn$OR_or_BETA[assn$STUDY_ACCESSION==acc]
  all_or <- as.logical(min(or_or_beta_this, na.rm=T)>1)
  all_beta <- as.logical(max(or_or_beta_this, na.rm=T)<=1)
  if (all_or) {
    n_all_or <- n_all_or + 1
    assn$oddsratio[assn$STUDY_ACCESSION==acc] <- or_or_beta_this
  }
  else if (all_beta) {
    n_all_beta <- n_all_beta + 1
    assn$beta[assn$STUDY_ACCESSION==acc] <- or_or_beta_this
  } else {
    n_both <- n_both + 1
  }
}
writeLines(sprintf("Studies with all values >1 (OR?): %d", n_all_or))
writeLines(sprintf("Studies with all values<=1 (beta?): %d", n_all_beta))
writeLines(sprintf("Studies with both values <=1 and >1: %d", n_both))

tag="beta"
qs <- quantile(assn[[tag]][!is.na(assn[[tag]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("%s %4s-ile: %9.1f",tag,names(qs),qs))

tag="oddsratio"
qs <- quantile(assn[[tag]][!is.na(assn[[tag]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("%s %4s-ile: %9.1f",tag,names(qs),qs))

write_delim(assn, ofile, delim="\t")

###
# Heuristic: all units include (in|de)crease
assn[["beta_units"]] <- sub("^.*\\] *", "", assn[["95%_CI_(TEXT)"]])
assn$beta_units[!grepl("(in|de)crease", assn$beta_units)] <- NA
tbl <- sort(table(assn$beta_units), decreasing = T)
writeLines(sprintf("%d. (N=%d) %s", 1:100, tbl[1:100], names(tbl)[1:100]))

#Top beta units
tbl <- as.data.frame(table(assn$beta_units))
colnames(tbl) <- c("beta_units", "Freq")
tbl <- tbl[order(-tbl$Freq),]
writeLines(sprintf("%18s: %3d", tbl$beta_units[1:20], tbl$Freq[1:20]))

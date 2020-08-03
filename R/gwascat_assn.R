#!/usr/bin/env Rscript
#############################################################################
### Separate OR from BETA, create new columns and write TSV.
### Heuristic: If all OR_or_BETA values for a study are >=1, assume OR.
### (Not necessarily true as BETA may be >1.) 
### Aha.  Column `95%_CI_(TEXT)` has units for BETAs.  But must be parsed.
### From http://www.ebi.ac.uk/gwas/docs/fileheaders:
### "95% CI (TEXT)*: Reported 95% confidence interval associated with
### strongest SNP risk allele, along with unit in the case of BETA-coefficients.
### If 95% CIs are not published, we estimate these using the standard error,
### where available."
#############################################################################
### Also get UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE for mapping
### confidence scoring. Also get CONTEXT, functionalClass in API. 
### INTERGENIC? "RISK ALLELE FREQUENCY"?  MERGED?
#############################################################################
### Output file next processed by tiga_gt_stats.R
#############################################################################
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  #ifile <- paste0(Sys.getenv("HOME"), "/../data/gwascatalog/data/gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv")
  ifile <- paste0(Sys.getenv("HOME"), "/../data/gwascatalog/data/gwas_catalog_v1.0.2-associations_e100_r2020-07-14.tsv")
  ofile <- "data/gwascat_assn.tsv"
} else {
  message("ERROR: Syntax: gwascat_assn.R ASSNFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), INTERGENIC=col_logical(), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(), `OR or BETA`=col_double(), UPSTREAM_GENE_DISTANCE=col_integer(), DOWNSTREAM_GENE_DISTANCE=col_integer(), DATE=col_date(format="%Y-%m-%d"), `DATE ADDED TO CATALOG`=col_date(format="%Y-%m-%d")))
setDT(assn)

setnames(assn, gsub("[ \\./]" ,"_", colnames(assn)))
setnames(assn, gsub("__" ,"_", colnames(assn)))
setnames(assn, gsub("_$" ,"", colnames(assn)))

assn <- assn[complete.cases(assn[, .(STUDY_ACCESSION, SNPS, DISEASE_TRAIT)]), ]

#Convert special chars.
for (tag in colnames(assn)) {
  if (typeof(assn[[tag]])=="character") {
    message(sprintf("NOTE: cleaning: %s", tag))
    assn[[tag]] <- iconv(assn[[tag]], from="latin1", to="UTF-8")
  }
}

message(sprintf("Total assn count: %6d", nrow(assn)))
message(sprintf("OR_or_BETA MISSING: %6d", nrow(assn[is.na(assn$OR_or_BETA),])))
message(sprintf("OR_or_BETA values: %6d", nrow(assn[!is.na(assn$OR_or_BETA),])))

tag="OR_or_BETA"
assn[[tag]] <- as.numeric(assn[[tag]])
qs <- quantile(assn[[tag]][!is.na(assn[[tag]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("%s %4s-ile: %9.1f", tag,names(qs), qs))

assn$oddsratio <- as.numeric(NA)
assn$beta <- as.numeric(NA)

staccs <- sort(unique(assn[, STUDY_ACCESSION]))
writeLines(sprintf("Studies, total: %6d", length(staccs)))
i <- 0
n_all_or <- 0
n_all_beta <- 0
n_both <- 0
for (stacc in staccs) {
  i <- i + 1
  or_or_beta_this <- assn[STUDY_ACCESSION==stacc, OR_or_BETA]
  all_or <- as.logical(min(or_or_beta_this, na.rm=T)>1)
  all_beta <- as.logical(max(or_or_beta_this, na.rm=T)<=1)
  if (all_or) {
    n_all_or <- n_all_or + 1
    set(assn, which(assn$STUDY_ACCESSION==stacc), "oddsratio",  or_or_beta_this)
  }
  else if (all_beta) {
    n_all_beta <- n_all_beta + 1
    set(assn, which(assn$STUDY_ACCESSION==stacc), "beta",  or_or_beta_this)
  } else {
    n_both <- n_both + 1
  }
}
writeLines(sprintf("Studies with all values >1 (OR?): %d", n_all_or))
writeLines(sprintf("Studies with all values<=1 (BETA?): %d", n_all_beta))
writeLines(sprintf("Studies with both values <=1 and >1: %d", n_both))
#

###
# Descriptive only, no more changes to assn.
# See also: gwascat_beta.R
###

# CONTEXT aka functionalClass (via API)
tbl <- sort(table(assn$CONTEXT), decreasing = T)
writeLines(sprintf("%d. (N=%d) %s", 1:100, tbl[1:100], names(tbl)[1:100]))

###
# Also see GENOTYPING_TECHNOLOGY
tbl <- as.data.frame(table(assn$GENOTYPING_TECHNOLOGY))
colnames(tbl) <- c("GENOTYPING_TECHNOLOGY", "Freq")
tbl <- tbl[order(-tbl$Freq),]
writeLines(sprintf("%5d: %s", tbl$Freq, tbl$GENOTYPING_TECHNOLOGY))

###
# Many missing (UP|DOWN)STREM_GENE_DISTANCE.
writeLines(sprintf("Associations with MAPPED_GENE: %d (%.1f%%)", nrow(assn[!is.na(MAPPED_GENE)]), 100*nrow(assn[!is.na(MAPPED_GENE)])/nrow(assn)))
writeLines(sprintf("Associations with MAPPED_GENE and (UP|DOWN)STREAM_GENE_DISTANCE: %d (%.1f%%)", 
  nrow(assn[!is.na(MAPPED_GENE) & (!is.na(UPSTREAM_GENE_DISTANCE) | !is.na(DOWNSTREAM_GENE_DISTANCE))]), 100*nrow(assn[!is.na(MAPPED_GENE) & (!is.na(UPSTREAM_GENE_DISTANCE) | !is.na(DOWNSTREAM_GENE_DISTANCE))])/nrow(assn)))
###
# Write file:
write_delim(assn, ofile, delim="\t")

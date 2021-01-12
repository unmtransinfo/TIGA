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
#############################################################################
### No filtering by this code.
#############################################################################
# See also: gwascat_assn_describe.R (descriptive only)
# See also: gwascat_beta.R
#############################################################################
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- paste0(Sys.getenv("HOME"), "/../data/GWASCatalog/releases/2020/12/16/gwas-catalog-associations_ontology-annotated.tsv")
  ofile <- "data/gwascat_assn.tsv"
} else {
  message("ERROR: Syntax: gwascat_assn.R ASSNFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), INTERGENIC=col_logical(), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(), `OR or BETA`=col_double(), UPSTREAM_GENE_DISTANCE=col_integer(), DOWNSTREAM_GENE_DISTANCE=col_integer(), DATE=col_date(format="%Y-%m-%d"), `DATE ADDED TO CATALOG`=col_date(format="%Y-%m-%d")))
setDT(assn)
message(sprintf("File (%s) rows: %d", sub("^.*/", "", ifile), nrow(assn)))

setnames(assn, gsub("[ \\./]" ,"_", colnames(assn)))
setnames(assn, gsub("__" ,"_", colnames(assn)))
setnames(assn, gsub("_$" ,"", colnames(assn)))

#Clean: convert special chars.
for (tag in colnames(assn)) {
  if (typeof(assn[[tag]])=="character") {
    assn[[tag]] <- iconv(assn[[tag]], from="latin1", to="UTF-8")
  }
}

###
message(sprintf("Studies in raw associations file: %d", assn[, uniqueN(STUDY_ACCESSION)]))
message(sprintf("Studies with SNPS: %d", assn[!is.na(SNPS), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Studies with DISEASE_TRAIT: %d", assn[!is.na(DISEASE_TRAIT), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Studies with OR_or_BETA: %d", assn[!is.na(OR_or_BETA), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Studies with SNPS and DISEASE_TRAIT and OR_or_BETA: %d", assn[(!is.na(DISEASE_TRAIT) & !is.na(SNPS) & !is.na(OR_or_BETA)), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Associations with OR_or_BETA values: %d (%.1f%%)", nrow(assn[!is.na(OR_or_BETA)]), 100*nrow(assn[!is.na(OR_or_BETA)])/nrow(assn)))
#
###
pval_threshold <- 5e-8
pval_mlog_threshold <- -log10(5e-8)
message("NOTE: PVALUE_MLOG (not P-VALUE) used for TIGA.")
message(sprintf("Studies with P-VALUE: %d", assn[!is.na(`P-VALUE`), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Associations with P-VALUE: %d; missing: %d", nrow(assn[!is.na(`P-VALUE`)]), nrow(assn[is.na(`P-VALUE`)])))
message(sprintf("Studies with PVALUE_MLOG: %d", assn[!is.na(PVALUE_MLOG), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Associations with PVALUE_MLOG: %d; missing: %d", nrow(assn[!is.na(PVALUE_MLOG)]), nrow(assn[is.na(PVALUE_MLOG)])))
message(sprintf("Studies with PVALUE_MLOG>=%.1f: %d", pval_mlog_threshold, assn[PVALUE_MLOG>pval_mlog_threshold, uniqueN(STUDY_ACCESSION)]))
message(sprintf("Associations with PVALUE_MLOG>=%.1f: %d", pval_mlog_threshold, nrow(assn[PVALUE_MLOG>pval_mlog_threshold])))
#
###
assn$risk_allele_freq <- assn$RISK_ALLELE_FREQUENCY #data.table bug/warning workaround.
assn[, risk_allele_freq := sub(" (.*)$", "", risk_allele_freq)]
assn[risk_allele_freq=="NR", risk_allele_freq := NA]
assn[, risk_allele_freq := as.double(risk_allele_freq)]
#
###
# Parsing oddsratio and beta from OR_or_BETA
assn[, oddsratio := as.double(NA)]
assn[, beta := as.double(NA)]
staccs <- sort(unique(assn[, STUDY_ACCESSION]))
#
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
message(sprintf("Studies with all values >1 (OR?): %d", n_all_or))
message(sprintf("Studies with all values<=1 (BETA?): %d", n_all_beta))
message(sprintf("Studies with both values <=1 and >1: %d", n_both))
message(sprintf("Studies with either OR or beta: %6d", assn[(!is.na(oddsratio) | !is.na(beta)), uniqueN(STUDY_ACCESSION)]))
message(sprintf("Studies with both OR and beta: %6d", assn[(!is.na(oddsratio) & !is.na(beta)), uniqueN(STUDY_ACCESSION)]))

###
# Write file:
write_delim(assn, ofile, delim="\t")


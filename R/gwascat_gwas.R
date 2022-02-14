#!/usr/bin/env Rscript
#############################################################################
### Clean studies file. No filtering by this code.
#############################################################################
library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  message("ERROR: Syntax: gwascat_gwas.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
rel_y <- args[1]
rel_m <- args[2]
rel_d <- args[3]
ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
#
ifile <- paste0(Sys.getenv("HOME"), sprintf("/../data/GWASCatalog/releases/%d/%02d/%02d/gwas-catalog-studies_ontology-annotated.tsv", rel_y, rel_m, rel_d))
ofile <- paste0(ODIR, "/gwascat_gwas.tsv")


message(sprintf("Input: %s", ifile))
message(sprintf("Output: %s", ofile))

# escape_double=F needed to parse all lines!
gwas <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), DATE=col_date(), `ASSOCIATION COUNT`=col_integer(), `DATE ADDED TO CATALOG`=col_date()), escape_double=F)
setDT(gwas)

setnames(gwas, gsub("[ \\./]", "_", colnames(gwas)))
setnames(gwas, gsub("__", "_", colnames(gwas)))
setnames(gwas, gsub("_$", "", colnames(gwas)))

#Clean: convert special chars.
for (tag in colnames(gwas)) {
  if (typeof(gwas[[tag]])=="character") {
    gwas[[tag]] <- iconv(gwas[[tag]], from="latin1", to="UTF-8")
  }
}

message(sprintf("Unique STUDY_ACCESSIONs: %d (%f per row)", gwas[!is.na(STUDY_ACCESSION), uniqueN(STUDY_ACCESSION)],  gwas[!is.na(STUDY_ACCESSION), uniqueN(STUDY_ACCESSION)]/nrow(gwas)))

# PUBMEDIDs
message(sprintf("PUBMEDIDs: %d (%.1f%% of studies)", gwas[!is.na(PUBMEDID), .N], 100*gwas[!is.na(PUBMEDID), .N]/nrow(gwas)))
message(sprintf("Unique PUBMEDIDs: %d (%f per study)", gwas[!is.na(PUBMEDID), uniqueN(PUBMEDID)], gwas[!is.na(PUBMEDID), uniqueN(PUBMEDID)]/gwas[!is.na(STUDY_ACCESSION), uniqueN(STUDY_ACCESSION)]))
# DATE is publication date.
pmid2gwas_counts <- gwas[, .(n_study = uniqueN(STUDY_ACCESSION)), by=c("PUBMEDID", "FIRST_AUTHOR", "JOURNAL", "DATE")]
setorder(pmid2gwas_counts, n_study)
for (n in unique(pmid2gwas_counts$n_study)) {
  if (n<=10) {
  message(sprintf("%4d publications contain %2d studies.", pmid2gwas_counts[n_study==n, .N], n))
  } else {
    message(sprintf("%4d publications contain [%2d-%d] studies.",
	pmid2gwas_counts[n_study>n, sum(.N)], n, max(pmid2gwas_counts$n_study)))
    break
  }
}
gwas2pmid_counts <- gwas[, .(n_pub = uniqueN(PUBMEDID)), by="STUDY_ACCESSION"]
message(sprintf("Studies with multiple publications: %d", gwas2pmid_counts[n_pub>1, uniqueN(STUDY_ACCESSION)]))

###
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
if (nrow(gwas[is.na(study_N)])>0) {
  writeLines(sprintf("UNPARSED INITIAL_SAMPLE_SIZE: \"%s\"", gwas[is.na(study_N), INITIAL_SAMPLE_SIZE]))
} else {
  message(sprintf("UNPARSED INITIAL_SAMPLE_SIZE: (None)"))
}

message(sprintf("Total gwas count: %6d", gwas[, uniqueN(STUDY_ACCESSION)]))

write_delim(gwas, ofile, delim="\t")
###

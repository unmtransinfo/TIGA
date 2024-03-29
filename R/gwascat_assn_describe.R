#!/usr/bin/env Rscript
###
# Descriptive only, no changes to assn.
# See also: gwascat_beta.R, gwascat_SNPs.Rmd
###
#############################################################################
library(readr)
library(data.table)
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)

GC_REL <- trimws(read_file("LATEST_RELEASE_GWC.txt"))
message(sprintf("LATEST_RELEASE_GWC: %s", GC_REL))
ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))

ifile <- ifelse((length(args)>0), args[1], paste0(ODIR, "/gwascat_assn.tsv"))
message(sprintf("Input ASSN file: %s", ifile))

assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(),
	INTERGENIC=col_logical(), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(),
	OR_or_BETA=col_double(), UPSTREAM_GENE_DISTANCE=col_integer(),
	DOWNSTREAM_GENE_DISTANCE=col_integer(), DATE=col_date(format="%Y-%m-%d"),
	DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"),
	risk_allele_freq=col_double(), oddsratio=col_double(), beta=col_double()))
setDT(assn)
message(sprintf("File (%s) rows: %d", sub("^.*/", "", ifile), nrow(assn)))

#
tag="OR_or_BETA"
qs <- quantile(assn[[tag]][!is.na(assn[[tag]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("%s %4s-ile: %9.1f", tag, names(qs), qs))
#

# CONTEXT (aka functionalClass, and cleaner, via API)
# Multiple values correspond with multiple MAPPED_GENE values.
fix_context <- function(context) {
  paste(sort(unique(unlist(strsplit(context, '\\s*;\\s*')))), collapse = " ; ")
}
assn[, CONTEXT := sapply(CONTEXT, fix_context)]
context_counts <- assn[, .(.N), by="CONTEXT"][order(-N)]
for (i in 1:nrow(context_counts)) {
  if (!grepl(" [;x] ", context_counts[i]$CONTEXT))
    message(sprintf("%d. (N=%d) %s", i, context_counts[i]$N, context_counts[i]$CONTEXT))
}
message(sprintf("Multiple: (N=%d)", sum(context_counts[grepl(" [;x] ", CONTEXT), .(N)])))

###
# GENOTYPING_TECHNOLOGY
tech_counts <- assn[, .(N_study = uniqueN(STUDY_ACCESSION)), by="GENOTYPING_TECHNOLOGY"][order(-N_study)]
print(tech_counts)

###
# Missing (UP|DOWN)STREM_GENE_DISTANCE implies within_gene.
message(sprintf("Associations with MAPPED_GENE: %d (%.1f%%)", nrow(assn[!is.na(MAPPED_GENE)]), 100*nrow(assn[!is.na(MAPPED_GENE)])/nrow(assn)))
message(sprintf("Associations within MAPPED_GENE: %d (%.1f%%)", 
  nrow(assn[!is.na(MAPPED_GENE) & is.na(UPSTREAM_GENE_DISTANCE) & is.na(DOWNSTREAM_GENE_DISTANCE)]), 100*nrow(assn[!is.na(MAPPED_GENE) & is.na(UPSTREAM_GENE_DISTANCE) & is.na(DOWNSTREAM_GENE_DISTANCE)])/nrow(assn)))
message(sprintf("Associations with MAPPED_GENE and (UP|DOWN)STREAM_GENE_DISTANCE: %d (%.1f%%)", 
  nrow(assn[!is.na(MAPPED_GENE) & (!is.na(UPSTREAM_GENE_DISTANCE) | !is.na(DOWNSTREAM_GENE_DISTANCE))]), 100*nrow(assn[!is.na(MAPPED_GENE) & (!is.na(UPSTREAM_GENE_DISTANCE) | !is.na(DOWNSTREAM_GENE_DISTANCE))])/nrow(assn)))
###

###
# SNP_GENE_IDS, UPSTREAM_GENE_ID and DOWNSTREAM_GENE_ID are Ensembl gene IDs, previously only available via API.
# Now we could use these for snp2gene mappings and avoid use of API. (To do.)

###
# SNPS formats (rs prefix means RefSNP. What of others?)
prefix_counts <- data.table(prefix = sub("[0-9].*$", "", assn$SNPS))[, .(.N), by=prefix][order(-N)]
message("SNP ID non-numeric-prefix counts:")
writeLines(sprintf("%2d. %6d: %s", 1:12, prefix_counts[1:12]$N, prefix_counts[1:12]$prefix))
for (i in 2:10) {
  pfx <- prefix_counts[i, prefix]
  message(sprintf("SNP examples, prefix \"%s\":\n\t%s", pfx, paste(collapse = "\n\t", sample(assn[sub("[0-9].*$", "", SNPS)==pfx, SNPS], 10))))
}

# SNPS delimiters (semantic difference?)
message(sprintf("SNPS multiple with delimiter \";\": %d", nrow(assn[grepl(";", SNPS)])))
message(sprintf("SNPS multiple with delimiter \"x\": %d", nrow(assn[grepl(" x ", SNPS)])))
message(sprintf("SNPS multiple with delimiter \",\": %d", nrow(assn[grepl(",", SNPS)])))


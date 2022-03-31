#!/usr/bin/env Rscript
#############################################################################
### Extract and map gene IDs from association file, produced by
### gwascat_assn.R. Report all counts. Accounting for missing genes
### important for users.
#############################################################################
library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  message("ERROR: Syntax: gwascat_gene.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
rel_y <- as.integer(args[1])
rel_m <- as.integer(args[2])
rel_d <- as.integer(args[3])
ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
#
ifile <- paste0(ODIR, "/gwascat_assn.tsv")
ofile <- paste0(ODIR, "/gwascat_gene.tsv")

writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))
#
assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), 
                                               INTERGENIC=col_logical(), 
                                               DATE=col_date(), 
                                               DATE_ADDED_TO_CATALOG=col_date(),
                                               UPSTREAM_GENE_DISTANCE=col_double(),
                                               DOWNSTREAM_GENE_DISTANCE=col_double(),
                                               `P-VALUE`=col_double(),
                                               PVALUE_MLOG=col_double(),
                                               OR_or_BETA=col_double(),
                                               oddsratio=col_double(),
                                               beta=col_double()))
setDT(assn)
#
gene_m <- assn[, .(STUDY_ACCESSION, MAPPED_GENE)]
gene_r <- assn[, .(STUDY_ACCESSION, `REPORTED_GENE(S)`)]
#
setnames(gene_m, c("STUDY_ACCESSION", "GENE"))
setnames(gene_r, c("STUDY_ACCESSION", "GENE"))
#
message(sprintf("Unique MAPPED_GENE values: %d", gene_m[, uniqueN(GENE)]))
message(sprintf("Unique MAPPED_GENE single-values: %d", gene_m[grepl("^[A-z0-9]+$", GENE), uniqueN(GENE)]))
message(sprintf("Unique MAPPED_GENE multi-values: %d", gene_m[!is.na(GENE) & !grepl("^[A-z0-9]+$", GENE), uniqueN(GENE)]))
gene_m_chars <- data.table(ch = unique(strsplit(paste(collapse="", gene_m[, unique(GENE)]), split="", fixed=T)[[1]]))
message(sprintf("MAPPED_GENE non-alphanumeric chars: '%s'", paste(collapse="' '", gene_m_chars[!grepl("[A-z0-9]", ch)]$ch)))
#
message(sprintf("Unique REPORTED_GENE values: %d", gene_r[, uniqueN(GENE)]))
message(sprintf("Unique REPORTED_GENE single-values: %d", gene_r[grepl("^[A-z0-9]+$", GENE), uniqueN(GENE)]))
message(sprintf("Unique REPORTED_GENE multi-values: %d", gene_r[!is.na(GENE) & !grepl("^[A-z0-9]+$", GENE), uniqueN(GENE)]))
gene_r_chars <- data.table(ch = unique(strsplit(paste(collapse="", gene_r[, unique(GENE)]), split="", fixed=T)[[1]]))
message(sprintf("REPORTED_GENE non-alphanumeric chars: '%s'", paste(collapse="' '", gene_r_chars[!grepl("[A-z0-9]", ch)]$ch)))
#
# Split separated vals on apparent delimiters. Must escape '-' and other special regex tokens.
# Problem: cases like: LOC105375010-3.8-1.4, H3.Y-LOC105374666, CAND1.11
setkey(gene_m, "STUDY_ACCESSION")
gene_m <- gene_m[, list(GENE=unlist(strsplit(GENE, "[ ,;@\\-]"))), by=STUDY_ACCESSION]
gene_m <- gene_m[!(is.na(gene_m$GENE) | (gene_m$GENE=="")),]
gene_m[['MAPPED_OR_REPORTED']] <- 'M'
#
setkey(gene_r, "STUDY_ACCESSION")
gene_r <- gene_r[, list(GENE=unlist(strsplit(GENE, "[ ,;@\\-]"))), by=STUDY_ACCESSION]
gene_r <- gene_r[!(is.na(gene_r$GENE) | (gene_r$GENE=="")),]
gene_r[['MAPPED_OR_REPORTED']] <- 'R'
#
writeLines(sprintf("Studies reporting pseudogene associations: %d", length(unique(gene_r$STUDY_ACCESSION[gene_r$GENE=="pseudogene"]))))
writeLines(sprintf("Studies reporting intergenic associations: %d", length(unique(gene_r$STUDY_ACCESSION[gene_r$GENE=="intergenic"]))))
#
gene_r <- gene_r[GENE != "pseudogene",]
gene_r <- gene_r[GENE != "intergenic",]
#
gene <- rbindlist(list(gene_r, gene_m))
setorder(gene, "STUDY_ACCESSION")

# Suspicious 1-char gene symbols:
writeLines(sprintf("Suspicious 1-char gene symbols (N=%d): '%s'", uniqueN(gene$GENE[!(nchar(gene$GENE)>1)]),
                   paste(collapse="', '", unique(gene$GENE[!(nchar(gene$GENE)>1)]))))
gene <- gene[nchar(gene$GENE)>1,]
#
writeLines(sprintf("Studies: %d", gene[, uniqueN(STUDY_ACCESSION)]))
writeLines(sprintf("Genes: %d", gene[, uniqueN(GENE)]))
writeLines(sprintf("Genes (MAPPED): %d", gene[MAPPED_OR_REPORTED=="M", uniqueN(GENE)]))
writeLines(sprintf("Genes (REPORTED): %d", gene[MAPPED_OR_REPORTED=="R", uniqueN(GENE)]))
#
write_delim(gene, ofile, delim="\t")


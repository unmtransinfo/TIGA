#!/usr/bin/env Rscript
#############################################################################
### Extract and map gene IDs from association file, produced by
### gwascat_assn.R. Report all counts. Accounting for missing genes
### important for users.
#############################################################################
library(readr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- "data/gwascat_assn.tsv"
  ofile <- "data/gwascat_gene.tsv"
} else {
  message("ERROR: Syntax: gwascat_gene.R ASSNFILE OFILE\n\t...or no args for defaults.")
  quit()
}
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
#
gene_m <- assn[, c("STUDY_ACCESSION", "MAPPED_GENE")]
gene_r <- assn[, c("STUDY_ACCESSION", "REPORTED_GENE(S)")]
#rm(assn)
colnames(gene_m) <- c("STUDY_ACCESSION", "GENE")
colnames(gene_r) <- c("STUDY_ACCESSION", "GENE")
gene_m$GENE <- gsub(" ", "", gene_m$GENE)
gene_r$GENE <- gsub(" ", "", gene_r$GENE)
#
# Split comma_or_hyphen_or_semicolon separated vals.
delim_regex <- "[,-;]"
setDT(gene_r, key="STUDY_ACCESSION")
gene_r <- gene_r[, list(GENE=unlist(strsplit(GENE, delim_regex))), by=STUDY_ACCESSION]
#
# Split hyphen_or_comma_or_semicolon separated vals.
delim_regex <- "[,-;]"
setDT(gene_m, key="STUDY_ACCESSION")
gene_m <- gene_m[,list(GENE=unlist(strsplit(GENE, delim_regex))), by=STUDY_ACCESSION]
#
gene_r <- gene_r[!(is.na(gene_r$GENE) | (gene_r$GENE=="")),]
gene_m <- gene_m[!(is.na(gene_m$GENE) | (gene_m$GENE=="")),]
#
gene_r[['MAPPED_OR_REPORTED']] <- 'R'
gene_m[['MAPPED_OR_REPORTED']] <- 'M'
#
writeLines(sprintf("Studies reporting pseudogene associations: %d", length(unique(gene_r$STUDY_ACCESSION[gene_r$GENE=="pseudogene"]))))
writeLines(sprintf("Studies reporting intergenic associations: %d", length(unique(gene_r$STUDY_ACCESSION[gene_r$GENE=="intergenic"]))))
#
gene_r <- gene_r[gene_r$GENE!="pseudogene",]
gene_r <- gene_r[gene_r$GENE!="intergenic",]
#
gene <- rbind(gene_r, gene_m)
gene <- gene[order(gene$STUDY_ACCESSION),]

# Suspicious 1-char gene symbols:
writeLines(sprintf("Suspicious 1-char gene symbols (N=%d): '%s'", uniqueN(gene$GENE[!(nchar(gene$GENE)>1)]),
                   paste(collapse="', '", unique(gene$GENE[!(nchar(gene$GENE)>1)]))))
gene <- gene[nchar(gene$GENE)>1,]
#
writeLines(sprintf("Studies: %d", length(unique(gene$STUDY_ACCESSION))))
writeLines(sprintf("Genes: %d", length(unique(gene$GENE))))
writeLines(sprintf("Genes (MAPPED): %d", length(unique(gene$GENE[gene$MAPPED_OR_REPORTED=="M"]))))
writeLines(sprintf("Genes (REPORTED): %d", length(unique(gene$GENE[gene$MAPPED_OR_REPORTED=="R"]))))
#
#
write_delim(gene, ofile, delim="\t")

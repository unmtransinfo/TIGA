#!/usr/bin/env Rscript
#############################################################################
### Extract and map gene IDs from association file, produced by
### gwascat_assn.R.
#############################################################################
library(readr)

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

assn <- read_delim(ifile, "\t")

gene_m <- assn[,c("STUDY_ACCESSION", "MAPPED_GENE")]
gene_r <- assn[,c("STUDY_ACCESSION", "REPORTED_GENE(S)")]
#rm(assn)
colnames(gene_m) <- c("STUDY_ACCESSION", "GENE")
colnames(gene_r) <- c("STUDY_ACCESSION", "GENE")
gene_m$GENE <- sub(" ", "", gene_m$GENE)
gene_r$GENE <- sub(" ", "", gene_r$GENE)
#
# Split comma_or_hyphen separated vals.
delim_regex <- "[,-]"
gene_r_multi <- gene_r[grepl(delim_regex, gene_r$GENE),]
gene_r <- gene_r[!grepl(delim_regex, gene_r$GENE),]
for (i in 1:nrow(gene_r_multi)) {
  acc <- gene_r_multi$STUDY_ACCESSION[i]
  genes <- strsplit(as.character(gene_r_multi$GENE[i]), delim_regex, perl=T)[[1]]
  gene_r <- rbind(gene_r, data.frame(STUDY_ACCESSION=acc, GENE=genes))
}

# Split hyphen_or_comma separated vals.
delim_regex <- "[,-]"
gene_m_multi <- gene_m[grepl(delim_regex, gene_m$GENE),]
gene_m <- gene_m[!grepl(delim_regex, gene_m$GENE),]
for (i in 1:nrow(gene_m_multi)) {
  acc <- gene_m_multi$STUDY_ACCESSION[i]
  genes <- strsplit(as.character(gene_m_multi$GENE[i]), delim_regex, perl=T)[[1]]
  gene_m <- rbind(gene_m, data.frame(STUDY_ACCESSION=acc, GENE=genes))
}

gene_r[['MAPPED_OR_REPORTED']] <- 'R'
gene_m[['MAPPED_OR_REPORTED']] <- 'M'

gene <- rbind(gene_r, gene_m)
gene <- gene[order(gene$STUDY_ACCESSION),]
#
write_delim(gene, ofile, delim="\t")

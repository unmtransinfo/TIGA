#!/usr/bin/env Rscript
#############################################################################
### gwax_gt_stats.R - Produce gt_stats.csv, for GWAS Explorer (GWAX) app.
### BETA not comparable for different units, thus currently only OR used.
#############################################################################
library(readr)
library(data.table, quietly=T)
#
Sys.time()
t0 <- proc.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==5) {
  (ifile_assn <- args[1])
  (ifile_snp2gene <- args[2])
  (ifile_trait <- args[3])
  (ifile_tcrd <- args[4])
  (ofile <- args[5])
} else if (length(args)==0) {
  ifile_assn <- "data/gwascat_assn.tsv"
  ifile_snp2gene <- "data/gwascat_snp2gene.tsv"
  ifile_trait <- "data/gwascat_trait.tsv"
  ifile_tcrd <- "data/tcrd_targets.csv"
  ofile <- "data/gt_stats.tsv"
} else {
  message("ERROR: Syntax: gwax_gt_stats.R ASSNFILE SNP2GENEFILE TRAITFILE TCRDFILE OFILE\n...or... no args for defaults")
  quit()
}
writeLines(sprintf("Input assn file: %s", ifile_assn))
writeLines(sprintf("Input snp2gene file: %s", ifile_snp2gene))
writeLines(sprintf("Input trait file: %s", ifile_trait))
writeLines(sprintf("Input TCRD file: %s", ifile_tcrd))
writeLines(sprintf("Output file: %s", ofile))
#
###
assn <- read_delim(ifile_assn, "\t", 
	col_types=cols(.default=col_character(),
	DATE=col_date(format="%Y-%m-%d"),
	DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"),
	SNP_ID_CURRENT=col_character()))
setDT(assn)
#
snp2gene <- read_delim(ifile_snp2gene, "\t", 
	col_types=cols(.default=col_character(),
	REPORTED_OR_MAPPED=col_factor(c("r","m","md","mu"))))
setDT(snp2gene)
#
trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
#
tcrd <- read_csv(ifile_tcrd, col_types=cols(.default=col_character()))
setDT(tcrd)
###
#Clean & transform:
names(trait) <- c("STUDY_ACCESSION","TRAIT","TRAIT_URI")
trait <- trait[!is.na(trait$TRAIT_URI)]
trait$TRAIT <- iconv(trait$TRAIT, from="latin1", to="UTF-8")
###
# Counts:
writeLines(sprintf("Studies: %d", uniqueN(assn$STUDY_ACCESSION)))
# MAPPED_GENE field may include chromosomal locations
writeLines(sprintf("MAPPED_GENE values: %d", uniqueN(assn$MAPPED_GENE)))
###
gsyms_tcrd <- unique(tcrd$protein_sym)
writeLines(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))

gsyms_gwax <- unique(snp2gene$GSYMB)
gsyms_common <- intersect(gsyms_gwax, gsyms_tcrd)
writeLines(sprintf("GSYMBs mapped to TCRD: %d", length(gsyms_common)))

tcrd <- merge(tcrd, data.table(gsym=gsyms_gwax, in_gwascat=rep(T, length(gsyms_gwax))),
	by.x="protein_sym", by.y="gsym", all.x=T, all.y=F)
tcrd$in_gwascat[is.na(tcrd$in_gwascat)] <- F
tcrd$idg2 <- as.logical(tcrd$idg2)
t2 <- table(tcrd$tdl[tcrd$in_gwascat])
writeLines(sprintf("%s: %d", names(t2), t2))
###
### g2t should have one row for each gene-snp-study-trait association.
g2t <- unique(snp2gene[, c("GSYMB", "SNP", "STUDY_ACCESSION")])
g2t <- merge(g2t, assn[, c("SNPS", "STUDY_ACCESSION","PVALUE_MLOG","OR_or_BETA","oddsratio","beta")], 
	all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"))

g2t <- merge(g2t, trait, all.x=F, all.y=F, by="STUDY_ACCESSION", allow.cartesian=T)
g2t <- g2t[!is.na(GSYMB)]
g2t <- g2t[!is.na(OR_or_BETA)]
g2t <- g2t[!grepl("(^LOC|^intergenic)", GSYMB)] #non-coding RNA, etc.

message(sprintf("DEBUG: with pvalue_mlog, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$PVALUE_MLOG),]),
	 uniqueN(g2t$GSYMB[!is.na(g2t$PVALUE_MLOG)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$PVALUE_MLOG)])))
message(sprintf("DEBUG: with or_or_beta, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$OR_or_BETA),]),
	 uniqueN(g2t$GSYMB[!is.na(g2t$OR_or_BETA)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$OR_or_BETA)])))
message(sprintf("DEBUG: with oddsratio, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$oddsratio),]),
	 uniqueN(g2t$GSYMB[!is.na(g2t$oddsratio)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$oddsratio)])))
message(sprintf("DEBUG: with beta, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$beta),]),
	 uniqueN(g2t$GSYMB[!is.na(g2t$beta)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$beta)])))
###
### GENE-TRAIT stats
### From g2t, create gt_stats table for TSV export.
### Too slow (~4h). Vectorize/optimize!
NROW <- 0
for (gsymb in unique(g2t$GSYMB)) {
  NROW <- NROW + uniqueN(g2t[GSYMB==gsymb, TRAIT_URI])
}
gt_stats <- data.table(gsymb=rep(NA, NROW), trait_uri=rep(NA, NROW),
	trait=rep(NA, NROW), n_study=rep(NA, NROW), n_snp=rep(NA, NROW),
	n_traits_g=rep(NA, NROW), 
	n_genes_t=rep(NA, NROW),
	pvalue_mlog_median=rep(NA, NROW),
	or_median=rep(NA, NROW))
#
writeLines(sprintf("Initialized rows to be populated: nrow(gt_stats) = %d\n", nrow(gt_stats)));
i_row <- 0 #gt_stats populated row count
for (gsymb in unique(g2t$GSYMB)) { #gene-loop
  n_traits_g <- uniqueN(g2t[GSYMB==gsymb, TRAIT_URI]) #n_traits for gene
  gt_stats$gsymb[i_row+1:i_row+n_traits_g] <- gsymb
  gt_stats$n_traits_g[i_row+1:i_row+n_traits_g] <- n_traits_g
  for (trait_uri in unique(g2t[GSYMB==gsymb, TRAIT_URI])) { #trait-loop
    if ((i_row%%10000)==0) {
      message(sprintf("i_row: %d / %d (%.1f%%) ; %s, elapsed: %.1fs", i_row, NROW, 100*i_row/NROW, Sys.time(), (proc.time()-t0)[3]))
    }
    i_row <- i_row + 1
    gt_stats$trait_uri[i_row] <- trait_uri
    gt_stats$trait[i_row] <- g2t[GSYMB==gsymb & TRAIT_URI==trait_uri, TRAIT][1]
    gt_stats$n_study[i_row] <- uniqueN(g2t[GSYMB==gsymb & TRAIT_URI==trait_uri, STUDY_ACCESSION])
    gt_stats$n_snp[i_row] <- uniqueN(g2t[GSYMB==gsymb & TRAIT_URI==trait_uri, SNP])
    gt_stats$pvalue_mlog_median[i_row] <- median(g2t[TRAIT_URI==trait_uri, PVALUE_MLOG], na.rm=T)
    gt_stats$or_median[i_row] <- median(g2t[GSYMB==gsymb & TRAIT_URI==trait_uri, oddsratio], na.rm=T)
  }
}
Sys.time()
for (uri in unique(g2t$TRAIT_URI)) {
  n_genes_t <- uniqueN(g2t[TRAIT_URI==uri, GSYMB])
  gt_stats$n_genes_t[trait_uri==uri] <- n_genes_t
}
writeLines(sprintf("Final: nrow(gt_stats) = %d\n", nrow(gt_stats)));
writeLines(sprintf("DEBUG: sum(is.na(gt_stats$gsymb)) = %d\n", sum(is.na(gt_stats$gsymb))));
writeLines(sprintf("DEBUG: sum(is.na(gt_stats$n_genes_t)) = %d\n", sum(is.na(gt_stats$n_genes_t))));
writeLines(sprintf("DEBUG: sum(is.na(gt_stats$n_traits_g)) = %d\n", sum(is.na(gt_stats$n_traits_g))));
gt_stats <- merge(gt_stats, tcrd[, c("protein_sym", "tdl", "fam", "idg2", "name")], by.x="gsymb", by.y="protein_sym", all.x=T, all.y=F)
write_delim(gt_stats, ofile, delim="\t")
###
message(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

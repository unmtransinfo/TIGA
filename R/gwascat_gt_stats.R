#!/usr/bin/env Rscript
#############################################################################
### gwascat_gt_stats.R - Produce gt_stats.csv
### 
### OR or BETA: Reported odds ratio or beta-coefficient associated with strongest 
### SNP risk allele. Note that if an OR <1 is reported this is inverted, along with 
### the reported allele, so that all ORs included in the Catalog are >1. Appropriate 
### unit and increase/decrease are included for beta coefficients.
### Ref: https://www.ebi.ac.uk/gwas/docs/methods
###
### Jeremy Yang
#############################################################################
library(readr)
library(dplyr, quietly=T)

t0 <- proc.time()

###
assn <- read_delim("data/gwascat_assn.tsv", "\t", 
	col_types=cols(DATE=col_date(format="%Y-%m-%d"),
		DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"),
		SNP_ID_CURRENT=col_character()))
snp2gene <- read_delim("data/gwascat_snp2gene.tsv", "\t", 
	col_types=cols(reported_or_mapped=col_factor(c("r","m","md","mu"))))
trait <- read_delim("data/gwascat_trait.tsv", "\t")

#Clean & transform:
colnames(trait) <- c("STUDY_ACCESSION","TRAIT","TRAIT_URI")
trait <- trait[!is.na(trait$TRAIT_URI),]
trait$TRAIT <- iconv(trait$TRAIT, from="latin1", to="UTF-8")

###
tcrd <- read_csv("data/tcrd_targets.csv")

gsyms_tcrd <- unique(tcrd$protein_sym)
print(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))

gsyms_gwascat <- unique(snp2gene$GSYMB)
gsyms_common <- intersect(gsyms_gwascat, gsyms_tcrd)
print(sprintf("GWASCat/TCRD geneSymbols in common: %d", length(gsyms_common)))

tcrd <- merge(tcrd, data.frame(gsym=gsyms_gwascat, in_gwascat=rep(T, length(gsyms_gwascat))),
              by.x="protein_sym", by.y="gsym", all.x=T, all.y=F)
tcrd$in_gwascat[is.na(tcrd$in_gwascat)] <- F
tcrd$idg2 <- as.logical(tcrd$idg2)
t2 <- table(tcrd$tdl[tcrd$in_gwascat])
print(sprintf("%s: %d\n", names(t2), t2))

###
### gene2trait should have one row for each gene-snp-study-trait association.
gene2trait <- unique(snp2gene[,c("GSYMB", "SNP", "STUDY_ACCESSION")])
gene2trait <- merge(gene2trait, assn[,c("SNPS", "STUDY_ACCESSION","PVALUE_MLOG","OR_or_BETA","oddsratio","beta")], 
  all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"))

gene2trait <- merge(gene2trait, trait, all.x=F, all.y=F, by="STUDY_ACCESSION")
gene2trait <- gene2trait[!is.na(gene2trait$GSYMB),]

writeLines(sprintf("DEBUG: with pvalue_mlog, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$PVALUE_MLOG),]),
	 length(unique(gene2trait$GSYMB[!is.na(gene2trait$PVALUE_MLOG)])),
	 length(unique(gene2trait$TRAIT[!is.na(gene2trait$PVALUE_MLOG)]))))
writeLines(sprintf("DEBUG: with or_or_beta, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$OR_or_BETA),]),
	 length(unique(gene2trait$GSYMB[!is.na(gene2trait$OR_or_BETA)])),
	 length(unique(gene2trait$TRAIT[!is.na(gene2trait$OR_or_BETA)]))))
writeLines(sprintf("DEBUG: with oddsratio, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$oddsratio),]),
	 length(unique(gene2trait$GSYMB[!is.na(gene2trait$oddsratio)])),
	 length(unique(gene2trait$TRAIT[!is.na(gene2trait$oddsratio)]))))
writeLines(sprintf("DEBUG: with beta, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$beta),]),
	 length(unique(gene2trait$GSYMB[!is.na(gene2trait$beta)])),
	 length(unique(gene2trait$TRAIT[!is.na(gene2trait$beta)]))))


###
### GENE-TRAIT stats
### From gene2trait, create gt_stats table for TSV export (~25min).
NCHUNK <- 10000
gt_stats_empty <- data.frame(gsymb=rep(NA,NCHUNK), trait_uri=rep(NA, NCHUNK),
	trait=rep(NA, NCHUNK), n_study=rep(NA, NCHUNK), n_snp=rep(NA, NCHUNK),
	n_traits_g=rep(NA, NCHUNK), 
	n_genes_t=rep(NA, NCHUNK),
	pvalue_mlog_median=rep(NA,NCHUNK),
	or_median=rep(NA,NCHUNK))
gt_stats <- NA
i <- 0
for (gsymb in unique(gene2trait$GSYMB))
{
  if (is.na(gt_stats)) { gt_stats <- gt_stats_empty; }
  gene2trait_g <- gene2trait[gene2trait$GSYMB==gsymb,]
  n_traits_g <- length(unique(gene2trait_g$TRAIT_URI)) #n_traits for gene
  for (trait_uri in unique(gene2trait_g$TRAIT_URI))
  {
    i <- i + 1
    if (i > nrow(gt_stats))
    {
      print(sprintf("nrow(gt_stats) = %d\n", nrow(gt_stats)));
      gt_stats <- rbind(gt_stats, gt_stats_empty)
    }
    gt_stats$gsymb[i] <- gsymb
    gt_stats$n_traits_g[i] <- n_traits_g
    gt_stats$trait_uri[i] <- trait_uri

    gene2trait_gt <- gene2trait_g[gene2trait_g$TRAIT_URI==trait_uri, ]

    gt_stats$n_study[i] <- length(unique(gene2trait_gt$STUDY_ACCESSION))
    gt_stats$trait[i] <- gene2trait_gt$TRAIT[gene2trait_gt$TRAIT_URI==trait_uri][1]
    gt_stats$n_snp[i] <- length(unique(gene2trait_gt$SNP))
    gt_stats$pvalue_mlog_median[i] <- median(gene2trait_gt$PVALUE_MLOG, na.rm=T)
    #gt_stats$or_median[i] <- median(gene2trait_gt$OR_or_BETA, na.rm=T)
    gt_stats$or_median[i] <- median(gene2trait_gt$oddsratio, na.rm=T)

# print(sprintf("%d. %7s - [%s] \"%s\" ; n_study = %d ; n_snp = %d ; p_median_nlog = %.1f\n", i, gsymb, sub("^.*/", "", trait_uri), gt_stats$trait[i], gt_stats$n_study[i], gt_stats$n_snp[i], gt_stats$p_median_nlog[i]))
  }
}
gt_stats <- gt_stats[!is.na(gt_stats$gsymb), ]
for (trait_uri in unique(gene2trait$TRAIT_URI))
{
  gene2trait_t <- gene2trait[gene2trait$TRAIT_URI==trait_uri, ]
  n_genes_t <- length(unique(gene2trait_t$GSYMB)) #n_genes for trait
  gt_stats$n_genes_t[gt_stats$trait_uri==trait_uri] <- n_genes_t
}
print(sprintf("Final: nrow(gt_stats) = %d\n", nrow(gt_stats)));

#Why so many unmapped symbols?
gt_stats <- merge(gt_stats, tcrd[,c("protein_sym", "tdl", "fam", "idg2", "name")], by.x="gsymb", by.y="protein_sym", all.x=T, all.y=F)

gt_stats <- gt_stats[!is.na(gt_stats$gsymb), ] #Where do these 7 bad rows come from?
write_delim(gt_stats, "data/gt_stats.tsv", delim="\t")

###
print(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

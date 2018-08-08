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
### 13 Dec 2017
#############################################################################
#library(RMySQL, quietly=T)
library(dplyr, quietly=T)

t0 <- proc.time()

###

library(readr)
gwascat_assn <- read_delim("~/projects/idg/gwas/data/gwascat_assn.csv", ",", escape_double=FALSE, trim_ws=TRUE,
	col_types=cols(DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"), SNP_ID_CURRENT=col_character()))
colnames(gwascat_assn) <- tolower(colnames(gwascat_assn))
gwascat_snp2gene <- read_delim("~/projects/idg/gwas/data/gwascat_snp2gene.tsv", "\t", escape_double=FALSE,  trim_ws=TRUE,
	col_types=cols(reported_or_mapped=col_factor(c("r","m","md","mu"))))
gwascat_trait <- read_delim("~/projects/idg/gwas/data/gwascat_trait.tsv",
	"\t", escape_double=FALSE, trim_ws=TRUE)
colnames(gwascat_trait) <- tolower(colnames(gwascat_trait))

#Clean & transform:
gwascat_trait <- gwascat_trait[!is.na(gwascat_trait$trait_uri),]


###
#dbcon <- dbConnect(MySQL(), host="juniper.health.unm.edu", dbname="tcrd")
#sql <- "SELECT
#  t.id AS \"tcrd_tid\",  t.name,  t.tdl,  t.fam,  t.idg2,
#  p.id AS \"tcrd_pid\",  p.sym AS \"protein_sym\",  p.geneid AS \"protein_geneid\"
#FROM
#  target t
#JOIN t2tc ON t.id = t2tc.target_id
#JOIN protein p ON t2tc.protein_id = p.id
#WHERE
#  p.sym IS NOT NULL"
#tcrd <- dbGetQuery(dbcon,sql)
#ok <- dbDisconnect(dbcon)
#write.csv(tcrd, file = "data/tcrd_targets.csv", row.names = F)
###
#Or use stored file:
tcrd <- read_csv("data/tcrd_targets.csv")

gsyms_tcrd <- unique(tcrd$protein_sym)
print(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))

gsyms_gwascat <- unique(gwascat_snp2gene$gsymb)
gsyms_common <- intersect(gsyms_gwascat, gsyms_tcrd)
print(sprintf("GWASCat/TCRD geneSymbols in common: %d", length(gsyms_common)))

tcrd <- merge(tcrd, data.frame(gsym=gsyms_gwascat, in_gwascat=rep(T, length(gsyms_gwascat))),
              by.x="protein_sym", by.y="gsym", all.x=T, all.y=F)
tcrd$in_gwascat[is.na(tcrd$in_gwascat)] <- F
tcrd$idg2 <- as.logical(tcrd$idg2)
t2 <- table(tcrd$tdl[tcrd$in_gwascat])
print(sprintf("%s: %d\n", names(t2), t2))

###
#gwascat_gene <- unique(gwascat_snp2gene[,c("study_accession", "gsymb")])
#gwascat_gene <- merge(gwascat_gene, tcrd[,c("protein_sym","name","idg2")], all.x=T, all.y=F, by.x="gsymb", by.y="protein_sym")

### gene2trait should have one row for each gene-snp-study-trait association.
gene2trait <- unique(gwascat_snp2gene[,c("gsymb", "snp", "study_accession")])
gene2trait <- merge(gene2trait, gwascat_assn[,c("snps", "study_accession","pvalue_mlog","or_or_beta","oddsratio","beta")], 
  all.x=T, all.y=F, by.x=c("snp", "study_accession"), by.y=c("snps", "study_accession"))

gene2trait <- merge(gene2trait, gwascat_trait, all.x=F, all.y=F, by="study_accession")
gene2trait <- gene2trait[!is.na(gene2trait$gsymb),]

writeLines(sprintf("DEBUG: with pvalue_mlog, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$pvalue_mlog),]),
	 length(unique(gene2trait$gsymb[!is.na(gene2trait$pvalue_mlog)])),
	 length(unique(gene2trait$trait[!is.na(gene2trait$pvalue_mlog)]))))
writeLines(sprintf("DEBUG: with or_or_beta, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$or_or_beta),]),
	 length(unique(gene2trait$gsymb[!is.na(gene2trait$or_or_beta)])),
	 length(unique(gene2trait$trait[!is.na(gene2trait$or_or_beta)]))))
writeLines(sprintf("DEBUG: with oddsratio, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$oddsratio),]),
	 length(unique(gene2trait$gsymb[!is.na(gene2trait$oddsratio)])),
	 length(unique(gene2trait$trait[!is.na(gene2trait$oddsratio)]))))
writeLines(sprintf("DEBUG: with beta, gene2trait: %d ; genes: %d ; traits: %d",
	 nrow(gene2trait[!is.na(gene2trait$beta),]),
	 length(unique(gene2trait$gsymb[!is.na(gene2trait$beta)])),
	 length(unique(gene2trait$trait[!is.na(gene2trait$beta)]))))

#gene2trait <- gene2trait[!is.na(gene2trait$pvalue_mlog),]
#gene2trait <- gene2trait[!is.na(gene2trait$or_or_beta),]
#gene2trait <- gene2trait[!is.na(gene2trait$oddsratio),]

### GENE-TRAIT stats
### From gene2trait, create gt_stats table for CSV export (~25min).
NCHUNK <- 10000
gt_stats_empty <- data.frame(gsymb=rep(NA,NCHUNK), trait_uri=rep(NA, NCHUNK),
	trait=rep(NA, NCHUNK), n_study=rep(NA, NCHUNK), n_snp=rep(NA, NCHUNK),
	n_traits_g=rep(NA, NCHUNK), 
	n_genes_t=rep(NA, NCHUNK),
	pvalue_mlog_median=rep(NA,NCHUNK),
	or_median=rep(NA,NCHUNK))
gt_stats <- NA
i <- 0
for (gsymb in unique(gene2trait$gsymb))
{
  if (is.na(gt_stats)) { gt_stats <- gt_stats_empty; }
  gene2trait_g <- gene2trait[gene2trait$gsymb==gsymb,]
  n_traits_g <- length(unique(gene2trait_g$trait_uri)) #n_traits for gene
  for (trait_uri in unique(gene2trait_g$trait_uri))
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

    gene2trait_gt <- gene2trait_g[gene2trait_g$trait_uri==trait_uri, ]

    gt_stats$n_study[i] <- length(unique(gene2trait_gt$study_accession))
    gt_stats$trait[i] <- gene2trait_gt$trait[gene2trait_gt$trait_uri==trait_uri][1]
    gt_stats$n_snp[i] <- length(unique(gene2trait_gt$snp))
    gt_stats$pvalue_mlog_median[i] <- median(gene2trait_gt$pvalue_mlog, na.rm=T)
    #gt_stats$or_median[i] <- median(gene2trait_gt$or_or_beta, na.rm=T)
    gt_stats$or_median[i] <- median(gene2trait_gt$oddsratio, na.rm=T)

# print(sprintf("%d. %7s - [%s] \"%s\" ; n_study = %d ; n_snp = %d ; p_median_nlog = %.1f\n", i, gsymb, sub("^.*/", "", trait_uri), gt_stats$trait[i], gt_stats$n_study[i], gt_stats$n_snp[i], gt_stats$p_median_nlog[i]))
  }
}
gt_stats <- gt_stats[!is.na(gt_stats$gsymb), ]
for (trait_uri in unique(gene2trait$trait_uri))
{
  gene2trait_t <- gene2trait[gene2trait$trait_uri==trait_uri, ]
  n_genes_t <- length(unique(gene2trait_t$gsymb)) #n_genes for trait
  gt_stats$n_genes_t[gt_stats$trait_uri==trait_uri] <- n_genes_t
}
print(sprintf("Final: nrow(gt_stats) = %d\n", nrow(gt_stats)));

#Why so many unmapped symbols?
gt_stats <- merge(gt_stats, tcrd[,c("protein_sym", "tdl", "fam", "idg2", "name")], by.x="gsymb", by.y="protein_sym", all.x=T, all.y=F)

gt_stats <- gt_stats[!is.na(gt_stats$gsymb), ] #Where do these 7 bad rows come from?
write.csv(gt_stats, file="data/gt_stats.csv", row.names=F)

###
print(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

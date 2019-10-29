#!/usr/bin/env Rscript
#############################################################################
### gwax_gt_stats.R - Produce gt_stats.csv, for GWAS Explorer (GWAX) app.
### ~50min
#############################################################################
### To do: For the gene-to-trait mode/view, there should be new columns
### for gene-specific statistics, and it should be very clear whether 
### stats are about (1) traits, (2) genes, or (3) gene-trait pairs.
### New:
###   * traitNstudy (#studies for trait; computed from gwascat_trait.tsv.)
###   * geneNstudy (#studies for gene; computed from gwascat_assn.tsv.)
### Renamed:
###   * n_traits_g to geneNtrait (#traits for gene)
###   * n_genes_t to traitNgene (#genes for trait)
###   * tcrdGeneSymbol to geneSymbol
###   * tcrdTargetName to geneName
###   * tcrdTargetFamily to geneFamily
###   * TDL to geneIdgTdl
###   * idgList to geneIdgList
###   * trait_uri to traitUri
###   * mu_score to muScore
###   * mu_rank to muRank
### Delete:
###   * nAbove
###   * nBelow
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
library(muStat)
#
Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==5) {
  (ifile_gwas	<- args[1])
  (ifile_counts	<- args[2])
  (ifile_assn	<- args[3])
  (ifile_snp2gene	<- args[4])
  (ifile_trait	<- args[5])
  (ifile_icite	<- args[6])
  (ifile_snps	<- args[7])
  (ifile_ensembl	<- args[8])
  (ifile_tcrd	<- args[9])
  (ofile	<- args[10])
} else if (length(args)==0) {
  ifile_gwas <- "data/gwascat_gwas.tsv"
  ifile_counts <- "data/gwascat_counts.tsv"
  ifile_assn <- "data/gwascat_assn.tsv"
  ifile_snp2gene <- "data/gwascat_snp2gene.tsv"
  ifile_trait <- "data/gwascat_trait.tsv"
  ifile_icite <- "data/gwascat_icite.tsv"
  ifile_snps <- "data/gwascat_Snps.tsv.gz" #API
  ifile_ensembl <- "data/gwascat_Snps_EnsemblInfo.tsv.gz"
  ifile_tcrd <- "data/tcrd_targets.tsv"
  ofile <- "data/gt_stats.tsv"
} else {
  message("ERROR: Syntax: gwax_gt_stats.R GWASFILE COUNTSFILE ASSNFILE SNP2GENEFILE TRAITFILE ICITEFILE TCRDFILE ENSEMBLFILE OFILE\n...or... no args for defaults")
  quit()
}
writeLines(sprintf("Input gwas file: %s", ifile_gwas))
writeLines(sprintf("Input counts file: %s", ifile_counts))
writeLines(sprintf("Input assn file: %s", ifile_assn))
writeLines(sprintf("Input snp2gene file: %s", ifile_snp2gene))
writeLines(sprintf("Input trait file: %s", ifile_trait))
writeLines(sprintf("Input iCite file: %s", ifile_icite))
writeLines(sprintf("Input TCRD file: %s", ifile_tcrd))
writeLines(sprintf("Input Ensembl file: %s", ifile_ensembl))
writeLines(sprintf("Output file: %s", ofile))
#
###
gwas <- read_delim(ifile_gwas, "\t", col_types=cols(.default=col_character(), 
	DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"),
	ASSOCIATION_COUNT=col_integer(), study_N=col_integer()))
setDT(gwas)
gwas_counts <- read_delim(ifile_counts, "\t", col_types=cols(.default=col_integer(), 
	study_accession=col_character()))
setDT(gwas_counts)
#
assn <- read_delim(ifile_assn, "\t", 
	col_types=cols(.default=col_character(),
	               UPSTREAM_GENE_DISTANCE=col_integer(),
	               DOWNSTREAM_GENE_DISTANCE=col_integer(),
	DATE=col_date(format="%Y-%m-%d"),
	DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"),
	SNP_ID_CURRENT=col_character()))
setDT(assn)
#
snp2gene <- read_delim(ifile_snp2gene, "\t", 
	col_types=cols(.default=col_character(),
	REPORTED_OR_MAPPED=col_factor(c("r","m","md","mu"))))
setDT(snp2gene)
snp2gene <- unique(snp2gene[REPORTED_OR_MAPPED!="r"]) #Ignore reported.
#
ensemblInfo <- read_delim(ifile_ensembl, "\t", col_types = cols(.default=col_character(), version=col_integer(), strand=col_integer(), start=col_integer(), end=col_integer()))
setDT(ensemblInfo)
#
trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
#
trait[, traitNstudy := uniqueN(STUDY_ACCESSION), by="MAPPED_TRAIT_URI"]
#
###
# Estimate RCR prior for new publications as median.
icite <- read_delim(ifile_icite, "\t", col_types=cols(.default=col_character(),
	relative_citation_ratio=col_double(), field_citation_rate=col_double(), citation_count=col_integer(),
	nih_percentile=col_double(), expected_citations_per_year=col_double(), citations_per_year=col_double(), year=col_integer()), na=c("", "NA", "None"))
setDT(icite)
rcr_median <- median(icite$relative_citation_ratio, na.rm=T)
icite[is.na(relative_citation_ratio) & (as.integer(format(Sys.time(), "%Y"))-year<2) , relative_citation_ratio := rcr_median]
#
icite_gwas <- merge(icite[, .(pmid, relative_citation_ratio, year)], gwas[, .(PUBMEDID, STUDY_ACCESSION)], by.x="pmid", by.y="PUBMEDID", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, gwas_counts[, .(study_accession, trait_count, gene_r_count, gene_m_count)], by.x="STUDY_ACCESSION", by.y="study_accession", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, icite_gwas[, .(study_perpmid_count = uniqueN(STUDY_ACCESSION)), by="pmid"], by="pmid")
# RCRAS = RCR-Aggregated-Score
icite_gwas[, rcras_pmid := (log2(relative_citation_ratio)+1)/study_perpmid_count]
icite_gwas <- icite_gwas[gene_r_count>0 | gene_m_count>0] #Need genes to be useful
icite_gwas[gene_r_count==0, gene_r_count := NA]
icite_gwas[gene_m_count==0, gene_m_count := NA]
icite_gwas[, rcras_study := 1/gene_r_count * rcras_pmid]
icite_gwas[is.na(rcras_study), rcras_study := 0]
icite_gwas <- icite_gwas[, .(pmid, STUDY_ACCESSION, year, relative_citation_ratio, rcras_pmid, rcras_study, trait_count, gene_r_count, gene_m_count, study_perpmid_count)]
setkey(icite_gwas, pmid, STUDY_ACCESSION)
###
# Link to Ensembl via IDs from Catalog API.
# (Switching to EnsemblIDs from gene symbols.)
studySnps <- read_delim(ifile_snps, "\t", col_types=cols(.default=col_character(), merged=col_logical(), lastUpdateDate=col_datetime(),   genomicContext_isIntergenic = col_logical(), genomicContext_isUpstream = col_logical(), genomicContext_isDownstream = col_logical(), genomicContext_distance = col_double(), genomicContext_isClosestGene = col_logical(), loc_chromosomePosition = col_double()))
setDT(studySnps)
ensemblInfo <- ensemblInfo[biotype=="protein_coding" & description!="novel transcript", .(ensemblId=id, ensemblName=display_name)][, protein_coding := T]
proteinSnps <- merge(studySnps[, .(rsId, ensemblId=gene_ensemblGeneIds, geneName=gene_geneName)], ensemblInfo, by="ensemblId", all.x=F)
snp2gene <- unique(merge(snp2gene, proteinSnps, by.x="SNP", by.y="rsId", all.x=F, all.y=F, allow.cartesian=T))
snp2gene[, `:=`(GSYMB=NULL, geneName=NULL, REPORTED_OR_MAPPED=NULL, protein_coding=NULL)]
snp2gene <- unique(snp2gene)
#
###
tcrd <- read_delim(ifile_tcrd, "\t", na=c("", "NA", "NULL"), col_types=cols(.default=col_character(), idgList=col_logical()))
setDT(tcrd)
###
#Clean & transform:
setnames(trait, c("STUDY_ACCESSION", "TRAIT", "TRAIT_URI", "id", "efo_label", "traitNstudy"))
trait <- trait[!is.na(trait$TRAIT_URI)]
trait$TRAIT <- iconv(trait$TRAIT, from="latin1", to="UTF-8")
###
# Counts:
writeLines(sprintf("Studies: %d", uniqueN(assn$STUDY_ACCESSION)))
# MAPPED_GENE field may include chromosomal locations
writeLines(sprintf("MAPPED_GENE values: %d", uniqueN(assn$MAPPED_GENE)))
#
assn_reported <- assn[, .(STUDY_ACCESSION, `REPORTED_GENE(S)`)]
assn_reported <- unique(assn_reported[, list(GENE=unlist(strsplit(`REPORTED_GENE(S)`, ", *"))), by=STUDY_ACCESSION])
writeLines(sprintf("REPORTED_GENE values: %d", uniqueN(assn_reported$GENE)))
###
#gsyms_tcrd <- unique(tcrd$tcrdGeneSymbol)
#writeLines(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))
ensgs_tcrd <- unique(tcrd$ensemblGeneId)
writeLines(sprintf("TCRD targets: %d ; ensemblGeneIds: %d", nrow(tcrd), length(ensgs_tcrd)))
#
ensgs_gwax <- unique(snp2gene$ensemblId)
ensgs_common <- intersect(ensgs_gwax, ensgs_tcrd)
writeLines(sprintf("EnsemblIds mapped to TCRD: %d", length(ensgs_common)))
#
tcrd <- merge(tcrd, data.table(ensg=ensgs_gwax, in_gwascat=rep(T, length(ensgs_gwax))),
	by.x="ensemblGeneId", by.y="ensg", all.x=T, all.y=F)
tcrd[, in_gwascat := !is.na(in_gwascat)]
writeLines(sprintf("IDG-List targets mapped by GWAS: %d", tcrd[(idgList & in_gwascat), .N]))
###
# Currently, we require OR and ignore BETA. (TO BE ADDRESSED.)
assn <- assn[!is.na(oddsratio)]
### g2t should have one row for each gene-snp-study-trait association.
#g2t <- unique(snp2gene[, c("GSYMB", "SNP", "STUDY_ACCESSION")])
g2t <- unique(snp2gene[, c("ensemblId", "ensemblName", "SNP", "STUDY_ACCESSION")])
g2t <- merge(g2t, gwas[, c("STUDY_ACCESSION", "study_N")], by="STUDY_ACCESSION", all.x=T, all.y=F)
g2t <- merge(g2t, assn[, c("SNPS", "STUDY_ACCESSION", "PVALUE_MLOG", "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE", "oddsratio", "beta")], all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"))
#
g2t <- merge(g2t, trait, all.x=F, all.y=F, by="STUDY_ACCESSION", allow.cartesian=T)
#
#g2t <- g2t[!is.na(GSYMB)]
#g2t <- g2t[!grepl("(^LOC|^intergenic)", GSYMB)] #non-coding RNA, etc.
g2t <- g2t[!is.na(ensemblId)] #None to delete.
#
message(sprintf("Deleting assocations lacking OR: %d", sum(is.na(g2t$oddsratio))))
g2t <- g2t[!is.na(oddsratio)] #Many to delete.
###
# Gene-distance weighting function.
g2t[, GDistWt := 2^(-pmin(g2t$UPSTREAM_GENE_DISTANCE, g2t$DOWNSTREAM_GENE_DISTANCE, na.rm=T)/5e4)]
message(sprintf("DEBUG: GDistWt: count: %d / %d (%.1f%%)", 
                sum(!is.na(g2t$GDistWt)), nrow(g2t), 100*sum(!is.na(g2t$GDistWt))/nrow(g2t)))
#
message(sprintf("DEBUG: with pvalue_mlog, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$PVALUE_MLOG),]),
	 uniqueN(g2t$ensemblId[!is.na(g2t$PVALUE_MLOG)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$PVALUE_MLOG)])))
message(sprintf("DEBUG: with OR, g2t: %d ; genes: %d ; traits: %d",
	 nrow(g2t[!is.na(g2t$oddsratio),]),
	 uniqueN(g2t$ensemblId[!is.na(g2t$oddsratio)]),
	 uniqueN(g2t$TRAIT[!is.na(g2t$oddsratio)])))
###
### GENE-TRAIT stats
### From g2t, create gt_stats table for TSV export.
### Too slow. Vectorize/optimize!
NROW <- 0
for (ensg in unique(g2t$ensemblId)) {
  NROW <- NROW + uniqueN(g2t[ensemblId==ensg, TRAIT_URI])
}
message(sprintf("Building gt_stats with NROW: %s", NROW))
gt_stats <- data.table(ensemblId=rep(NA, NROW), 
	trait_uri=rep(NA, NROW), trait=rep(NA, NROW), 
	n_study=as.integer(rep(NA, NROW)), 
	n_snp=as.integer(rep(NA, NROW)),
	n_snpw=as.numeric(rep(NA, NROW)),
	geneNtrait=as.integer(rep(NA, NROW)),
	geneNstudy=as.integer(rep(NA, NROW)),
	traitNgene=as.integer(rep(NA, NROW)),
	traitNstudy=as.integer(rep(NA, NROW)),
	pvalue_mlog_median=as.numeric(rep(NA, NROW)),
	or_median=as.numeric(rep(NA, NROW)),
	study_N_mean=as.numeric(rep(NA, NROW)),
	rcras=rep(NA, NROW)
	)
#
message(sprintf("Initialized rows to be populated: nrow(gt_stats) = %d", nrow(gt_stats)))
i_row <- 0 #gt_stats populated row count
t0 <- proc.time()
# gene-loop:
for (ensg in unique(g2t$ensemblId)) {
  # trait-loop:
  for (trait_uri in unique(g2t[ensemblId==ensg, TRAIT_URI])) {
    if ((i_row%%10000)==0) {
      message(sprintf("i_row: %d / %d (%.1f%%) ; %s, elapsed: %.1fs", i_row, NROW, 100*i_row/NROW, Sys.time(), (proc.time()-t0)[3]))
    }
    i_row <- i_row + 1
    #
    gt_stats$ensemblId[i_row] <- ensg
    gt_stats$trait_uri[i_row] <- trait_uri
    gt_stats$trait[i_row] <- g2t[ensemblId==ensg & TRAIT_URI==trait_uri, TRAIT][1]
    gt_stats$traitNstudy[i_row] <- g2t[TRAIT_URI==trait_uri, traitNstudy][1]
    gt_stats$n_study[i_row] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])
    gt_stats$pvalue_mlog_median[i_row] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, PVALUE_MLOG], na.rm=T)
    gt_stats$or_median[i_row] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, oddsratio], na.rm=T)
    gt_stats$study_N_mean[i_row] <- round(mean(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, study_N], na.rm=T), 1)
    gt_stats$n_snp[i_row] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, SNP])
    # Deduplicate (group-by) SNPs for `n_snpw` computation. 
    gt_stats$n_snpw[i_row] <- sum(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, .(GDistWt = median(GDistWt)), by="SNP"][, GDistWt], na.rm=T)
    #
    rcras <- 0.0
    for (stacc in unique(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])) {
      grc <- gwas_counts[study_accession==stacc, gene_r_count]
      if (is.na(grc) | length(grc)==0 | grc==0) { next; }
      rcras_study <- 0.0
      for (pmid_this in icite_gwas[STUDY_ACCESSION==stacc, pmid]) {
        spp <- icite_gwas[STUDY_ACCESSION==stacc & pmid==pmid_this, study_perpmid_count]
        if (is.na(spp) | spp==0) { next; }
        rcr <- icite_gwas[STUDY_ACCESSION==stacc & pmid==pmid_this, relative_citation_ratio]
        if (is.na(rcr) | rcr==0.0) { next; }
        rcras_study <- rcras_study + log2(rcr+1)/spp
      }
      rcras <- rcras + rcras_study
    }
    gt_stats$rcras[i_row] <- rcras
  }
}
gt_stats[, traitNgene := .N, by="trait_uri"]
gt_stats[, geneNtrait := .N, by="ensemblId"]
gt_stats[, geneNstudy := uniqueN(STUDY_ACCESSION), by="ensemblId"]
#
gt_stats$or_median <- round(as.double(gt_stats$or_median), 3)
gt_stats$pvalue_mlog_median <- round(as.double(gt_stats$pvalue_mlog_median), 3)
gt_stats$rcras <- round(as.double(gt_stats$rcras), 3)
gt_stats$n_snpw <- round(as.double(gt_stats$n_snpw), 3)
#
message(sprintf("%s, elapsed: %.1fs", Sys.time(), (proc.time()-t0)[3]))
message(sprintf("Final: nrow(gt_stats) = %d", nrow(gt_stats)))
message(sprintf("gene (ensemblId) count: %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("trait (trait_uri) count: %d", uniqueN(gt_stats$trait_uri)))
message(sprintf("DEBUG: traitNgene: [%d,%d]", min(gt_stats$traitNgene), max(gt_stats$traitNgene)))
message(sprintf("DEBUG: traitNstudy: [%d,%d]", min(gt_stats$traitNstudy), max(gt_stats$traitNstudy)))
message(sprintf("DEBUG: geneNtrait: [%d,%d]", min(gt_stats$geneNtrait), max(gt_stats$geneNtrait)))
message(sprintf("DEBUG: pvalue_mlog_median: [%.2f,%.2f]", min(gt_stats$pvalue_mlog_median, na.rm=T), max(gt_stats$pvalue_mlog_median, na.rm=T)))
message(sprintf("DEBUG: or_median: [%.2f,%.2f]", min(gt_stats$or_median, na.rm=T), max(gt_stats$or_median, na.rm=T)))
message(sprintf("DEBUG: study_N_mean: [%.1f,%.1f]", min(gt_stats$study_N_mean, na.rm=T), max(gt_stats$study_N_mean, na.rm=T)))
message(sprintf("DEBUG: rcras: [%.2f,%.2f]", min(gt_stats$rcras, na.rm=T), max(gt_stats$rcras, na.rm=T)))
message(sprintf("DEBUG: n_snpw: [%.2f,%.2f]", min(gt_stats$n_snpw, na.rm=T), max(gt_stats$n_snpw, na.rm=T)))
#
gt_stats <- merge(gt_stats, tcrd[, c("ensemblGeneId", "tcrdGeneSymbol", "TDL", "tcrdTargetFamily", "idgList", "tcrdTargetName")], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F)
setnames(gt_stats,
	old=c("tcrdGeneSymbol", "tcrdTargetName", "tcrdTargetFamily", "TDL", "idgList"),
	new=c("geneSymbol", "geneName", "geneFamily", "geneIdgTdl", "geneIdgList"))
#
#############################################################################
# Beginning of former gwax_gt_stats_mu.R

###
# Multivariable non-parametric ranking via &mu; scores.
# &mu; scores
###
# Non-dominated solutions are not inferior to any other case at any variable.
# A &mu; score is defined as the number of lower cases minus the number of higher.
# The resulting ranking is the useful result, not so much the score itself.
# Using muStat package.
###
# Issue: cases globally superior in one variable and inferior in one variable
# have nAbove=0 and nBelow=0 and muScore=0. What should be the rank?
###
#
gt_stats <- gt_stats[!is.na(ensemblId)] #Should be no-op.
gt_stats <- gt_stats[!is.na(geneIdgTdl)] #Removes non-protein-coding.
gt_stats <- gt_stats[!is.na(or_median)] #OR required until beta -able.
###
#
message(sprintf("Genes (ensemblIDs): %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("Genes (symbols): %d", uniqueN(gt_stats$geneSymbol)))
message(sprintf("Traits in dataset: %d", uniqueN(gt_stats$trait_uri)))
message(sprintf("G-T associations in dataset: %d", nrow(gt_stats)))

# Use inverse of n\_traits\_g and n\_genes\_t so bigger is better as needed for muStat.
gt_stats[, geneNtrait_inv := 1 / geneNtrait]
gt_stats[, traitNgene_inv := 1 / traitNgene]

# We are interested in rankings __for a given trait__. So for each trait,
# convert to matrix for muStat::mu.GE().
# The (i,j) entry of GE matrix is 1 if \code{x_i >= x_j}, 0 otherwise. The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.

gt_stats[, `:=`(muScore=as.integer(NA), muRank=as.integer(NA))]

ii <- 0
for (trait_this in unique(gt_stats$trait)) {
  gtmat <- as.matrix(gt_stats[trait==trait_this, .(n_study, n_snp, n_snpw, geneNtrait_inv, traitNgene_inv, pvalue_mlog_median, or_median, rcras)])
  ii <- ii + 1
  message(sprintf("[%d / %d] (N_gene: %3d) \"%s\"", ii, uniqueN(gt_stats$trait), dim(gtmat)[1], trait_this))
  if (dim(gtmat)[1]<2) { #No ranking for singleton.
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt_stats[trait==trait_this, .(geneName)]
  # Bug in muStat?
  if (sum(is.na(sums$score))>0) {
    message(sprintf("DEBUG: missing score count: %d", sum(is.na(sums$score))))
    badrows <- is.na(sums$score)
    print(sums[badrows]) #DEBUG
    sums[is.na(score), score := nAbove - nBelow]
    message(sprintf("DEBUG: missing score count: %d (FIXED?)", sum(is.na(sums$score))))
    print(sums[badrows]) #DEBUG
  }
  sums <- setorder(sums, -score, na.last=T)
  sums[, rank := 1:nrow(sums)]
  gt_stats[trait==trait_this]$muScore <- sums$score
  gt_stats[trait==trait_this]$muRank <- sums$rank
  #gt_stats[trait==trait_this]$nAbove <- sums$nAbove #unneeded
  #gt_stats[trait==trait_this]$nBelow <- sums$nBelow #unneeded
}
#
gt_stats[, `:=`(geneNtrait_inv = NULL, traitNgene_inv = NULL)]
#
write_delim(gt_stats, ofile, delim="\t")
#
message(sprintf("NOTE: elapsed time: %.2fs",(proc.time()-t0)[3]))

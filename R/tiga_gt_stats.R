#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT stats
### tiga_gt_stats.R - Produce gt_stats.csv, for TIGA Shiny app.
### ~20hr
#############################################################################
### SEE FOR ALL INPUT WORKFLOWS: Go_gwascat_GetData.sh
#############################################################################
# Multivariable non-parametric ranking via &mu; scores.
###
# Non-dominated solutions are not inferior to any other case at any variable.
# A &mu; score is defined as the number of lower cases minus the number of higher.
# The resulting ranking is the useful result, not so much the score itself.
# Using muStat package.
###
# mu.GE: greater than or equal to
# mu.AND: Logical AND GEs for all variables
###
# Issue: cases globally superior in one variable and inferior in one variable
# have nAbove=0 and nBelow=0 and muScore=0. What should be the rank?
# If a case has nAbove=0 and nBelow=0 then weight=0 and mu.Sums() returns score = NA.
# This is arbitrary so we assign score = nBelow = nAbove = 0 - 0 = 0.
# ?mu.Sums:
# score  = (nB-nA) * ifelse(weight==0, NA, 1)
#############################################################################
# Writes gt_provenance file with TRAIT_URI, ensemblId, STUDY_ACCESSION and PUBMEDID.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
library(muStat)
#
t_start <- Sys.time()
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
  (ofile_prov	<- args[11])
} else if (length(args)==0) {
  ifile_gwas <- "data/gwascat_gwas.tsv.gz"	#gwascat_gwas.R
  ifile_counts <- "data/gwascat_counts.tsv"	#gwascat_counts.R
  ifile_assn <- "data/gwascat_assn.tsv"	#gwascat_assn.R
  ifile_snp2gene <- "data/gwascat_snp2gene.tsv" #snp2gene_mapped.pl, snp2gene_reported.pl
  ifile_trait <- "data/gwascat_trait.tsv"	#gwascat_trait.R
  ifile_icite <- "data/gwascat_icite.tsv" #BioClients.icite API
  ifile_snps <- "data/gwascat_Snps.tsv.gz" #BioClients.gwascatalog API (for addl data)
  ifile_ensembl <- "data/gwascat_Snps_EnsemblInfo.tsv.gz" #BioClients.ensembl API
  ifile_tcrd <- "data/tcrd_targets.tsv" #BioClients.idg API
  ofile <- "data/gt_stats.tsv.gz"
  ofile_prov <- "data/gt_provenance.tsv.gz"
} else {
  message("ERROR: Syntax: tiga_gt_stats.R GWASFILE COUNTSFILE ASSNFILE SNP2GENEFILE TRAITFILE ICITEFILE TCRDFILE ENSEMBLFILE OFILE OFILE_PROVENANCE\n...or... no args for defaults")
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
writeLines(sprintf("Output provenance file: %s", ofile_prov))
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
trait[, traitNstudy := uniqueN(STUDY_ACCESSION), by="MAPPED_TRAIT_URI"]
setnames(trait, old=c("MAPPED_TRAIT_URI", "MAPPED_TRAIT"), new=c("TRAIT_URI", "TRAIT"))
trait <- trait[!is.na(trait$TRAIT_URI)]
trait$TRAIT <- iconv(trait$TRAIT, from="latin1", to="UTF-8")
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
# (Switched to EnsemblIDs from gene symbols.)
studySnps <- read_delim(ifile_snps, "\t", col_types=cols(.default=col_character(), merged=col_logical(), lastUpdateDate=col_datetime(),   genomicContext_isIntergenic = col_logical(), genomicContext_isUpstream = col_logical(), genomicContext_isDownstream = col_logical(), genomicContext_distance = col_double(), genomicContext_isClosestGene = col_logical(), loc_chromosomePosition = col_double()))
setDT(studySnps)
ensemblInfo <- ensemblInfo[biotype=="protein_coding" & description!="novel transcript", .(ensemblId=id, ensemblSymb=display_name)][, protein_coding := T]
proteinSnps <- merge(studySnps[, .(rsId, ensemblId=gene_ensemblGeneIds, geneName=gene_geneName)], ensemblInfo, by="ensemblId", all.x=F)
snp2gene <- unique(merge(snp2gene, proteinSnps, by.x="SNP", by.y="rsId", all.x=F, all.y=F, allow.cartesian=T))
snp2gene[, `:=`(GSYMB=NULL, geneName=NULL, REPORTED_OR_MAPPED=NULL, protein_coding=NULL)]
snp2gene <- unique(snp2gene)
#
###
tcrd <- read_delim(ifile_tcrd, "\t", na=c("", "NA", "NULL"), col_types=cols(.default=col_character(), idgList=col_logical()))
setDT(tcrd)
###

###
# Counts:
writeLines(sprintf("Studies: %d", uniqueN(assn$STUDY_ACCESSION)))
#
###
# Reported genes ignored by TIGA.
#assn_reported <- assn[, .(STUDY_ACCESSION, `REPORTED_GENE(S)`)]
#assn_reported <- unique(assn_reported[, list(GENE=unlist(strsplit(`REPORTED_GENE(S)`, ", *"))), by=STUDY_ACCESSION])
###
ensgs_tcrd <- unique(tcrd$ensemblGeneId)
writeLines(sprintf("TCRD targets: %d ; ensemblGeneIds: %d", nrow(tcrd), length(ensgs_tcrd)))
#
ensgs_tiga <- unique(snp2gene$ensemblId)
ensgs_common <- intersect(ensgs_tiga, ensgs_tcrd)
writeLines(sprintf("EnsemblIds mapped to TCRD: %d", length(ensgs_common)))
#
tcrd <- merge(tcrd, data.table(ensg=ensgs_tiga, in_gwascat=rep(T, length(ensgs_tiga))),
	by.x="ensemblGeneId", by.y="ensg", all.x=T, all.y=F)
tcrd[, in_gwascat := !is.na(in_gwascat)]
writeLines(sprintf("IDG-List targets mapped by GWAS: %d", tcrd[(idgList & in_gwascat), .N]))
#print(tcrd[tcrdGeneSymbol == "CDK1" | ensemblGeneId == snp2gene[ensemblSymb == "CDK1"]$ensemblId[1]]) #DEBUG
#
### g2t should have one row for each gene-snp-study-trait association. Only 1 PUBMEDID per study.
g2t <- unique(snp2gene[, .(ensemblId, ensemblSymb, SNP, STUDY_ACCESSION)])
g2t <- merge(g2t, gwas[, .(STUDY_ACCESSION, PUBMEDID, study_N)], by="STUDY_ACCESSION", all.x=T, all.y=F)
g2t <- merge(g2t, assn[, .(SNPS, STUDY_ACCESSION, PVALUE_MLOG, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, oddsratio, beta)], all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"))
g2t <- merge(g2t, trait, all.x=F, all.y=F, by="STUDY_ACCESSION", allow.cartesian=T)
#
###
# Initially, we required OR and ignored BETA.
# This filtered many studies, traits, and genes.
print(sprintf("Studies without and with OR filter: %d -> %d (-%d; -%.1f%%)", 
  g2t[, uniqueN(STUDY_ACCESSION)], g2t[!is.na(oddsratio), uniqueN(STUDY_ACCESSION)],
  g2t[, uniqueN(STUDY_ACCESSION)] - g2t[!is.na(oddsratio), uniqueN(STUDY_ACCESSION)],
  100*(g2t[, uniqueN(STUDY_ACCESSION)] - g2t[!is.na(oddsratio), uniqueN(STUDY_ACCESSION)])/g2t[, uniqueN(STUDY_ACCESSION)]
))

print(sprintf("Traits without and with OR filter: %d -> %d (-%d; -%.1f%%)", 
  g2t[, uniqueN(TRAIT_URI)], g2t[!is.na(oddsratio), uniqueN(TRAIT_URI)],
  g2t[, uniqueN(TRAIT_URI)] - g2t[!is.na(oddsratio), uniqueN(TRAIT_URI)],
  100*(g2t[, uniqueN(TRAIT_URI)] - g2t[!is.na(oddsratio), uniqueN(TRAIT_URI)])/g2t[, uniqueN(TRAIT_URI)]
  ))

print(sprintf("Genes without and with OR filter: %d -> %d (-%d; -%.1f%%)", 
  g2t[, uniqueN(ensemblId)], g2t[!is.na(oddsratio), uniqueN(ensemblId)],
  g2t[, uniqueN(ensemblId)] - g2t[!is.na(oddsratio), uniqueN(ensemblId)],
  100*(g2t[, uniqueN(ensemblId)] - g2t[!is.na(oddsratio), uniqueN(ensemblId)])/g2t[, uniqueN(ensemblId)]
  ))
#
#badrows <- is.na(g2t$oddsratio)
#reason_txt <- "Missing OR"
badrows <- (is.na(g2t$oddsratio) & is.na(g2t$beta))
reason_txt <- "Missing both OR and beta"
print(sprintf("badrows: %d", sum(badrows)))
###
# Write files accounting for filtered studies, traits and genes.
filtered_studies <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows]$STUDY_ACCESSION, g2t[!badrows]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F))
filtered_studies[, reason := reason_txt]
print(sprintf("Studies removed by filter (%s): %d", reason_txt, filtered_studies[, uniqueN(STUDY_ACCESSION)]))
write_delim(filtered_studies, "data/filtered_studies.tsv.gz", delim="\t")
#
filtered_traits <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows]$TRAIT_URI, g2t[!badrows]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F))
filtered_traits[, reason := reason_txt]
print(sprintf("Traits removed by filter (%s): %d", reason_txt, filtered_traits[, uniqueN(TRAIT_URI)]))
write_delim(filtered_traits, "data/filtered_traits.tsv.gz", delim="\t")
#
filtered_genes <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows]$ensemblId, g2t[!badrows]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F))
filtered_genes[, reason := reason_txt]
print(sprintf("Genes removed by filter (%s): %d", reason_txt, filtered_genes[, uniqueN(ensemblId)]))
filtered_genes <- merge(filtered_genes, tcrd[, .(ensemblGeneId, geneName=tcrdTargetName, geneFamily=tcrdTargetFamily, TDL)], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F)
filtered_genes <- filtered_genes[, .(ensemblId, ensemblSymb, geneName, geneFamily, TDL, reason)]
write_delim(filtered_genes, "data/filtered_genes.tsv.gz", delim="\t")
#
message(sprintf("Filtering assocations (%s): %d", reason_txt, sum(badrows)))
g2t <- g2t[!badrows] #Many filtered.
#
#stop("DEBUG: STOP.")
#
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
# GENE-TRAIT provenance
gt_prov <- NULL
i_row_prov <- 0
#
# gene-loop:
for (ensg in unique(g2t$ensemblId)) {
  geneNstudy <- g2t[ensemblId==ensg, uniqueN(STUDY_ACCESSION)]
  # trait-loop:
  for (trait_uri in unique(g2t[ensemblId==ensg, TRAIT_URI])) {
    studies_this <- unique(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, .(STUDY_ACCESSION, PUBMEDID)])
      gt_prov_this <- data.table(ensemblId=rep(ensg, nrow(studies_this)), TRAIT_URI=rep(trait_uri, nrow(studies_this)),
		  STUDY_ACCESSION=studies_this[, STUDY_ACCESSION], PUBMEDID=studies_this[, PUBMEDID])
    if (is.null(gt_prov)) {
      gt_prov <- gt_prov_this
    } else {
      gt_prov <- rbindlist(list(gt_prov, gt_prov_this))
    }
  }
}
#
gt_prov[, efoId := sub("^.*/", "", TRAIT_URI)]
write_delim(gt_prov, ofile_prov, delim="\t")
writeLines(sprintf("Output provenance file written: %s", ofile_prov))
#
###
### GENE-TRAIT stats
### One row per unique gene-trait pair.
### From g2t, create gt_stats table for TSV export.
### Slow. Vectorize/optimize!?
NROW <- nrow(unique(g2t[, .(ensemblId, TRAIT_URI)]))
message(sprintf("Building gt_stats with NROW: %s", NROW))
gt_stats <- data.table(ensemblId=rep(NA, NROW), 
	efoId=rep(NA, NROW), trait=rep(NA, NROW), 
	n_study=as.integer(rep(NA, NROW)), 
	n_snp=as.integer(rep(NA, NROW)),
	n_snpw=as.numeric(rep(NA, NROW)),
	geneNtrait=as.integer(rep(NA, NROW)),
	geneNstudy=as.integer(rep(NA, NROW)),
	traitNgene=as.integer(rep(NA, NROW)),
	traitNstudy=as.integer(rep(NA, NROW)),
	pvalue_mlog_median=as.numeric(rep(NA, NROW)),
	or_median=as.numeric(rep(NA, NROW)),
	n_beta=as.integer(rep(NA, NROW)), #simple count of beta values
	study_N_mean=as.numeric(rep(NA, NROW)),
	rcras=rep(NA, NROW)
	)
#
message(sprintf("Initialized rows to be populated: nrow(gt_stats) = %d", nrow(gt_stats)))
#
#stop("DEBUG: STOP.")
#
i_row <- 0 #gt_stats populated row count
t0 <- proc.time()
# gene-loop:
for (ensg in unique(g2t$ensemblId)) {
  geneNstudy <- g2t[ensemblId==ensg, uniqueN(STUDY_ACCESSION)]
  # trait-loop:
  for (trait_uri in unique(g2t[ensemblId==ensg, TRAIT_URI])) {
    if ((i_row%%10000)==0) {
      message(sprintf("i_row: %d / %d (%.1f%%) ; %s, elapsed: %.1fs", i_row, NROW, 100*i_row/NROW, Sys.time(), (proc.time()-t0)[3]))
    }
    i_row <- i_row + 1
    #
    gt_stats$ensemblId[i_row] <- ensg
    gt_stats$efoId[i_row] <- sub("^.*/", "", trait_uri)
    gt_stats$geneNstudy[i_row] <- geneNstudy
    gt_stats$trait[i_row] <- g2t[ensemblId==ensg & TRAIT_URI==trait_uri, TRAIT][1]
    gt_stats$traitNstudy[i_row] <- g2t[TRAIT_URI==trait_uri, traitNstudy][1]
    gt_stats$n_study[i_row] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])
    gt_stats$pvalue_mlog_median[i_row] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, PVALUE_MLOG], na.rm=T)
    gt_stats$or_median[i_row] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, oddsratio], na.rm=T) #NA if no ORs
    gt_stats$n_beta[i_row] <- g2t[ensemblId==ensg & TRAIT_URI==trait_uri & !is.na(beta), .N] #0 if no betas
    gt_stats$study_N_mean[i_row] <- round(mean(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, study_N], na.rm=T), 1)
    gt_stats$n_snp[i_row] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, SNP])
    # Deduplicate (group-by) SNPs for `n_snpw` computation. 
    gt_stats$n_snpw[i_row] <- sum(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, .(GDistWt = median(GDistWt, na.rm=T)), by="SNP"][, GDistWt], na.rm=T)
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
gt_stats[, traitNgene := .N, by="efoId"]
gt_stats[, geneNtrait := .N, by="ensemblId"]
#
gt_stats$or_median <- round(as.double(gt_stats$or_median), 3)
gt_stats$pvalue_mlog_median <- round(as.double(gt_stats$pvalue_mlog_median), 3)
gt_stats$rcras <- round(as.double(gt_stats$rcras), 3)
gt_stats$n_snpw <- round(as.double(gt_stats$n_snpw), 3)
#
message(sprintf("%s, elapsed: %.1fs", Sys.time(), (proc.time()-t0)[3]))
message(sprintf("Final: nrow(gt_stats) = %d", nrow(gt_stats)))
message(sprintf("gene (ensemblId) count: %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("trait (efoId) count: %d", uniqueN(gt_stats$efoId)))
message(sprintf("traitNgene: [%d,%d]", min(gt_stats$traitNgene), max(gt_stats$traitNgene)))
message(sprintf("traitNstudy: [%d,%d]", min(gt_stats$traitNstudy), max(gt_stats$traitNstudy)))
message(sprintf("geneNtrait: [%d,%d]", min(gt_stats$geneNtrait), max(gt_stats$geneNtrait)))
message(sprintf("pvalue_mlog_median: [%.2f,%.2f]", min(gt_stats$pvalue_mlog_median, na.rm=T), max(gt_stats$pvalue_mlog_median, na.rm=T)))
message(sprintf("or_median: [%.2f,%.2f]", min(gt_stats$or_median, na.rm=T), max(gt_stats$or_median, na.rm=T)))
message(sprintf("n_beta: [%d,%d]", min(gt_stats$n_beta, na.rm=T), max(gt_stats$n_beta, na.rm=T)))
message(sprintf("study_N_mean: [%.1f,%.1f]", min(gt_stats$study_N_mean, na.rm=T), max(gt_stats$study_N_mean, na.rm=T)))
message(sprintf("rcras: [%.2f,%.2f]", min(gt_stats$rcras, na.rm=T), max(gt_stats$rcras, na.rm=T)))
message(sprintf("n_snpw: [%.2f,%.2f]", min(gt_stats$n_snpw, na.rm=T), max(gt_stats$n_snpw, na.rm=T)))
#
gt_stats <- merge(gt_stats, tcrd[, c("ensemblGeneId", "tcrdGeneSymbol", "TDL", "tcrdTargetFamily", "idgList", "tcrdTargetName")], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F)
setnames(gt_stats,
	old=c("tcrdGeneSymbol", "tcrdTargetName", "tcrdTargetFamily", "TDL", "idgList"),
	new=c("geneSymbol", "geneName", "geneFamily", "geneIdgTdl", "geneIdgList"))
#
gt_stats <- gt_stats[!is.na(ensemblId)] #Should be no-op.
gt_stats <- gt_stats[!is.na(geneIdgTdl)] #Non-protein-coding removed by IDG TDL requirement.
#gt_stats <- gt_stats[!is.na(or_median)] #OR required until beta -able.
###
#
message(sprintf("Genes (ensemblIDs): %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("Genes (symbols): %d", uniqueN(gt_stats$geneSymbol)))
message(sprintf("Traits in dataset: %d", uniqueN(gt_stats$efoId)))
message(sprintf("G-T associations in dataset: %d", nrow(gt_stats)))

###
# &mu; scores/rankings
#   * __for a given trait__
#   * __for a given gene__
###
#
# Use inverse of n\_traits\_g and n\_genes\_t so bigger is better as needed for muStat.
gt_stats[, geneNtrait_inv := 1 / geneNtrait]
gt_stats[, traitNgene_inv := 1 / traitNgene]

# Save pre-mu to file for debugging.
write_delim(gt_stats, "data/tmp.tsv.gz", delim="\t")

# For each (trait|gene), convert to matrix for muStat::mu.GE().
# The (i,j) entry of GE matrix is 1 if \code{x_i >= x_j}, 0 otherwise.
# The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.
###
# Gene &mu; scores:
# Some or_median will be NA, since beta included.
gt_stats[, `:=`(geneMuScore=as.integer(NA), geneMuRank=as.integer(NA))]
ii <- 0
for (efoId_this in unique(gt_stats$efoId)) {
  gtmat <- as.matrix(gt_stats[efoId==efoId_this, .(n_study, n_snp, n_snpw, geneNtrait_inv, traitNgene_inv, pvalue_mlog_median, or_median, n_beta, rcras)])
  ii <- ii + 1
  trait_this <- gt_stats[efoId==efoId_this, trait][1]
  message(sprintf("[%d / %d] (N_gene: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$efoId), dim(gtmat)[1], efoId_this, trait_this))
  if (dim(gtmat)[1]<2) { #Singletons
    gt_stats[efoId==efoId_this]$geneMuScore <- 0
    gt_stats[efoId==efoId_this]$geneMuRank <- 1
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt_stats[efoId==efoId_this, .(geneName)]
  sums[is.na(score) & weight==0, score := 0] # Override NA scores
  sums[order(-score, -weight), rank := 1:nrow(sums)]
  gt_stats[efoId==efoId_this]$geneMuScore <- sums$score
  gt_stats[efoId==efoId_this]$geneMuRank <- sums$rank
}
#
###
# Trait &mu; scores:
gt_stats[, `:=`(traitMuScore=as.integer(NA), traitMuRank=as.integer(NA))]
ii <- 0
for (ensemblId_this in unique(gt_stats$ensemblId)) {
  gtmat <- as.matrix(gt_stats[ensemblId==ensemblId_this, .(n_study, n_snp, n_snpw, geneNtrait_inv, traitNgene_inv, pvalue_mlog_median, or_median, n_beta, rcras)])
  ii <- ii + 1
  geneSymbol_this <- gt_stats[ensemblId==ensemblId_this, geneSymbol][1]
  message(sprintf("[%d / %d] (N_trait: %3d) %s:\"%s\"", ii, uniqueN(gt_stats$ensemblId), dim(gtmat)[1], ensemblId_this, geneSymbol_this))
  if (dim(gtmat)[1]<2) { #Singletons
    gt_stats[ensemblId==ensemblId_this]$traitMuScore <- 0
    gt_stats[ensemblId==ensemblId_this]$traitMuRank <- 1
    next;
  }
  ge <- mu.GE(gtmat)
  sums <- mu.Sums(mu.AND(ge)) # Logical AND GEs for all variables 
  setDT(sums)
  sums$name <- gt_stats[ensemblId==ensemblId_this, .(geneName)]
  sums[is.na(score) & weight==0, score := 0] # Override NA scores
  sums[order(-score, -weight), rank := 1:nrow(sums)]
  gt_stats[ensemblId==ensemblId_this]$traitMuScore <- sums$score
  gt_stats[ensemblId==ensemblId_this]$traitMuRank <- sums$rank
}
###
gt_stats[, `:=`(geneNtrait_inv = NULL, traitNgene_inv = NULL)]
#
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2fs (%.2f %s)", (proc.time()-t0)[3], t_elapsed, attr(t_elapsed, "units")))

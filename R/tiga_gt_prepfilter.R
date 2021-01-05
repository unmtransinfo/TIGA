#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT PREPFILTER - preprocess and filter
### 
### (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
### (2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
### (2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
### (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR FULL WORKFLOW: Go_gwascat_GetData.sh
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
message(t_start)
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==10) {
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
  ifile_gwas <- "data/gwascat_gwas.tsv.gz"	#gwascat_gwas.R
  ifile_counts <- "data/gwascat_counts.tsv"	#Go_gwascat_DbCreate.sh
  ifile_assn <- "data/gwascat_assn.tsv"	#gwascat_assn.R
  ifile_snp2gene <- "data/gwascat_snp2gene.tsv" #snp2gene_mapped.pl, snp2gene_reported.pl
  ifile_trait <- "data/gwascat_trait.tsv"	#gwascat_trait.R
  ifile_icite <- "data/gwascat_icite.tsv" #BioClients.icite API
  ifile_snps <- "data/gwascat_Snps.tsv.gz" #BioClients.gwascatalog API (for addl data)
  ifile_ensembl <- "data/gwascat_EnsemblInfo.tsv.gz" #BioClients.ensembl API
  ifile_tcrd <- "data/tcrd_targets.tsv" #BioClients.idg API
  ofile <- "data/gt_prepfilter.Rdata"
} else {
  message("ERROR: Syntax: tiga_gt_prepfilter.R GWASFILE COUNTSFILE ASSNFILE SNP2GENEFILE TRAITFILE ICITEFILE TCRDFILE ENSEMBLFILE OFILE \n...or... no args for defaults")
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
writeLines(sprintf("Output prepfilter file: %s", ofile))
#
###
gwas <- read_delim(ifile_gwas, "\t", col_types=cols(.default=col_character(), DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"), ASSOCIATION_COUNT=col_integer(), study_N=col_integer()))
setDT(gwas)
gwas_counts <- read_delim(ifile_counts, "\t", col_types=cols(.default=col_integer(), study_accession=col_character()))
setDT(gwas_counts)
#
assn <- read_delim(ifile_assn, "\t", col_types=cols(.default=col_character(), 
	DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"), 
	UPSTREAM_GENE_DISTANCE=col_integer(), DOWNSTREAM_GENE_DISTANCE=col_integer(), 
	oddsratio=col_double(), beta=col_double(), OR_or_BETA=col_double(), 
	PVALUE_MLOG=col_double(), `P-VALUE`=col_double()))
setDT(assn)
###
# Missing (UP|DOWN)STREAM_GENE_DISTANCE implies within MAPPED_GENE, thus distances zero for gene-distance weighting function, GDistWt.
assn[(!is.na(MAPPED_GENE) & is.na(UPSTREAM_GENE_DISTANCE) & is.na(DOWNSTREAM_GENE_DISTANCE)), `:=`(UPSTREAM_GENE_DISTANCE=0, DOWNSTREAM_GENE_DISTANCE=0)]
###
#
snp2gene <- read_delim(ifile_snp2gene, "\t", col_types=cols(.default=col_character(), REPORTED_OR_MAPPED=col_factor(c("r", "m", "md", "mu"))))
setDT(snp2gene)
snp2gene <- unique(snp2gene[REPORTED_OR_MAPPED!="r"]) #Ignore reported.
#
ensemblInfo <- read_delim(ifile_ensembl, "\t", col_types = cols(.default=col_character(), version=col_integer(), strand=col_integer(), start=col_integer(), end=col_integer()))
setDT(ensemblInfo)
#
trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, old=c("MAPPED_TRAIT_URI", "MAPPED_TRAIT"), new=c("TRAIT_URI", "TRAIT"))
trait[, TRAIT := iconv(TRAIT, from="latin1", to="UTF-8")]
#
###
# Estimate RCR for new publications as global median. Also any undefined RCR.
icite <- read_delim(ifile_icite, "\t", col_types=cols(.default=col_character(), relative_citation_ratio=col_double(), field_citation_rate=col_double(), citation_count=col_integer(), nih_percentile=col_double(), expected_citations_per_year=col_double(), citations_per_year=col_double(), year=col_integer()), na=c("", "NA", "None"))
setDT(icite)
rcr_median <- median(icite$relative_citation_ratio, na.rm=T) #global median ignoring NAs.
year_this <- as.integer(format(Sys.time(), "%Y"))
message(sprintf("Estimating undefined RCR for new publications [%d-%d] (and any undefined) as global median.", year_this-1, year_this))
icite[is.na(relative_citation_ratio) & (year>=year_this-1) , relative_citation_ratio := rcr_median]
icite[is.na(relative_citation_ratio), relative_citation_ratio := rcr_median]
#
icite_gwas <- merge(icite[, .(pmid, relative_citation_ratio, year)], gwas[, .(PUBMEDID, STUDY_ACCESSION)], by.x="pmid", by.y="PUBMEDID", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, gwas_counts[, .(study_accession, trait_count, gene_r_count, gene_m_count)], by.x="STUDY_ACCESSION", by.y="study_accession", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, icite_gwas[, .(study_perpmid_count = uniqueN(STUDY_ACCESSION)), by="pmid"], by="pmid")
###
# RCRAS = RCR-Aggregated-Score
icite_gwas[, rcras_pmid := log2(relative_citation_ratio+1)/study_perpmid_count]
icite_gwas <- icite_gwas[gene_r_count>0 | gene_m_count>0] #Need genes to be useful
icite_gwas[gene_r_count==0, gene_r_count := NA]
icite_gwas[gene_m_count==0, gene_m_count := NA]
icite_gwas[, rcras_study := 1/gene_r_count * rcras_pmid]
icite_gwas[is.na(rcras_study), rcras_study := 0]
icite_gwas <- icite_gwas[, .(pmid, STUDY_ACCESSION, year, relative_citation_ratio, rcras_pmid, rcras_study, trait_count, gene_r_count, gene_m_count, study_perpmid_count)]
setkey(icite_gwas, pmid, STUDY_ACCESSION) #Ensures unique study-pmid pairs each row.
###
# Link to Ensembl via IDs from Catalog.
# (In 2020 (not 2018), EnsemblIDs available in downloads, so API not required.)
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
writeLines(sprintf("Association studies: %d", uniqueN(assn$STUDY_ACCESSION)))
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
tcrd <- merge(tcrd, data.table(ensg=ensgs_tiga, in_gwascat=T), by.x="ensemblGeneId", by.y="ensg", all.x=T, all.y=F)
tcrd[, in_gwascat := !is.na(in_gwascat)]
writeLines(sprintf("IDG-List targets mapped by GWAS: %d", tcrd[(idgList & in_gwascat), .N]))
#
### g2t should have one row for each gene-snp-study-trait association. Only 1 PUBMEDID per study.
g2t <- unique(snp2gene[, .(ensemblId, ensemblSymb, SNP, STUDY_ACCESSION)])
g2t <- merge(g2t, gwas[, .(STUDY_ACCESSION, PUBMEDID, study_N)], by="STUDY_ACCESSION", all.x=T, all.y=F)
g2t <- merge(g2t, assn[, .(SNPS, STUDY_ACCESSION, PVALUE_MLOG, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, oddsratio, beta)], all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"))
#
trait[, traitNstudy := fifelse(is.na(TRAIT_URI), as.integer(NA), uniqueN(STUDY_ACCESSION)), by="TRAIT_URI"]
g2t <- merge(g2t, trait, all.x=F, all.y=F, by="STUDY_ACCESSION", allow.cartesian=T)
#
###
# FILTERS:
###
# (1) TRAIT_URI filter: require MAPPED_TRAIT_URI (EFO)
reason_txt <- "Missing MAPPED_TRAIT_URI"
badrows_traituri <- (is.na(g2t$TRAIT_URI))
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_traituri, na.rm=T)))
filtered_studies <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_traituri]$STUDY_ACCESSION, g2t[!badrows_traituri]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F))
filtered_studies[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies[, uniqueN(STUDY_ACCESSION)], filtered_studies[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
#
filtered_traits <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_traituri]$TRAIT_URI, g2t[!badrows_traituri]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F))
filtered_traits[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits[, uniqueN(TRAIT_URI)], filtered_traits[, uniqueN(TRAIT_URI)], 100*filtered_traits[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
#
filtered_genes <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_traituri]$ensemblId, g2t[!badrows_traituri]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F))
filtered_genes[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes[, uniqueN(ensemblId)], filtered_genes[, uniqueN(ensemblId)], 100*filtered_genes[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
###
# (2) P-value filter: must exceed standard threshold for significance (pvalue <= 1e-8, pvalue_mlog >= 8)
pval_threshold <- 5e-8
pval_mlog_threshold <- -log10(5e-8)
reason_txt <- sprintf("P-value fails significance threshold (%g)", pval_threshold)
badrows_pval <- (g2t$PVALUE_MLOG<pval_mlog_threshold)
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_pval, na.rm=T)))
filtered_studies_pval <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_pval]$STUDY_ACCESSION, g2t[!badrows_pval]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F))
filtered_studies_pval[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies_pval[, uniqueN(STUDY_ACCESSION)], filtered_studies_pval[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies_pval[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_pval))
#
filtered_traits_pval <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_pval]$TRAIT_URI, g2t[!badrows_pval]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F))
filtered_traits_pval[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits_pval[, uniqueN(TRAIT_URI)],
filtered_traits_pval[, uniqueN(TRAIT_URI)], 100*filtered_traits_pval[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
filtered_traits <- rbindlist(list(filtered_traits, filtered_traits_pval))
#
filtered_genes_pval <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_pval]$ensemblId, g2t[!badrows_pval]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F))
filtered_genes_pval[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes_pval[, uniqueN(ensemblId)],
filtered_genes_pval[, uniqueN(ensemblId)], 100*filtered_genes_pval[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
filtered_genes <- rbindlist(list(filtered_genes, filtered_genes_pval))
###
# (3) Effect size filter: either OR or beta required.
badrows_effect <- (is.na(g2t$oddsratio) & is.na(g2t$beta))
reason_txt <- "Missing both OR and beta"
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_effect, na.rm=T)))
###
filtered_studies_effect <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_effect]$STUDY_ACCESSION, g2t[!badrows_effect]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F))
filtered_studies_effect[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies_effect[, uniqueN(STUDY_ACCESSION)], filtered_studies_effect[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies_effect[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_effect))
#
filtered_traits_effect <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_effect]$TRAIT_URI, g2t[!badrows_effect]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F))
filtered_traits_effect[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits_effect[, uniqueN(TRAIT_URI)],
filtered_traits_effect[, uniqueN(TRAIT_URI)], 100*filtered_traits_effect[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
filtered_traits <- rbindlist(list(filtered_traits, filtered_traits_effect))
#
filtered_genes_effect <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_effect]$ensemblId, g2t[!badrows_effect]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F))
filtered_genes_effect[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes_effect[, uniqueN(ensemblId)],
filtered_genes_effect[, uniqueN(ensemblId)], 100*filtered_genes_effect[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
filtered_genes <- rbindlist(list(filtered_genes, filtered_genes_effect))
#
###
filtered_studies <- filtered_studies[, reason := paste0(reason, collapse=" ; "), by=c("STUDY_ACCESSION", "STUDY")]
filtered_traits <- filtered_traits[, reason := paste0(reason, collapse=" ; "), by=c("TRAIT_URI", "TRAIT")]
filtered_genes <- filtered_genes[, reason := paste0(reason, collapse=" ; "), by=c("ensemblId", "ensemblSymb")]
#
filtered_genes <- merge(filtered_genes, tcrd[, .(ensemblGeneId, geneName=tcrdTargetName, geneFamily=tcrdTargetFamily, TDL)], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F)
filtered_genes <- filtered_genes[, .(ensemblId, ensemblSymb, geneName, geneFamily, TDL, reason)]
#
# Write files accounting for filtered studies, traits and genes.
write_delim(filtered_studies, "data/filtered_studies.tsv", delim="\t")
write_delim(filtered_traits, "data/filtered_traits.tsv", delim="\t")
write_delim(filtered_genes, "data/filtered_genes.tsv", delim="\t")
#
badrows <- (badrows_traituri | badrows_pval | badrows_effect)
message(sprintf("Filtered associations (total): %d -> %d (-%d; -%.1f%%)", nrow(g2t), nrow(g2t)-sum(badrows), sum(badrows), 100*sum(badrows)/nrow(g2t)))
g2t <- g2t[!badrows] #Many filtered.
#
###
#
message(sprintf("Final: nrow(g2t) = %d", nrow(g2t)))
message(sprintf("G-T associations in dataset: %d", nrow(unique(g2t[, .(ensemblId, efoId)]))))
message(sprintf("Study (STUDY_ACCESSION) count: %d", uniqueN(g2t$STUDY_ACCESSION)))
message(sprintf("Gene (ensemblId) count: %d", uniqueN(g2t$ensemblId)))
message(sprintf("Trait (efoId) count: %d", uniqueN(g2t$efoId)))
message(sprintf("PVALUE_MLOG (instance) count: %d", nrow(g2t[!is.na(PVALUE_MLOG)])))
message(sprintf("oddsratio (instance) count: %d", nrow(g2t[!is.na(oddsratio)])))
message(sprintf("beta (instance) count: %d", nrow(g2t[!is.na(beta)])))
#
###
# Save to file.
save(g2t, tcrd, icite_gwas, file=ofile)
writeLines(sprintf("Output Rdata file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT PREPFILTER - preprocess and filter
### 
### (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
### (2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
### (2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
### (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR FULL WORKFLOW: Go_TIGA_Workflow.sh
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
message(t_start)
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
#
GC_REL <- trimws(read_file("LATEST_RELEASE.txt"))
message(sprintf("LATEST_RELEASE: %s", GC_REL))
ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))
#
ifile_gwas <-	ifelse((length(args)>0), args[1], paste0(ODIR, "/gwascat_gwas.tsv")) #gwascat_gwas.R
ifile_counts <-	ifelse((length(args)>1), args[2], paste0(ODIR, "/gwascat_gwas_counts.tsv")) # NOW tiga_gwas_counts.py (OLD Go_gwascat_DbCreate.sh)
ifile_assn <-	ifelse((length(args)>2), args[3], paste0(ODIR, "/gwascat_assn.tsv")) #gwascat_assn.R
ifile_snp2gene <-ifelse((length(args)>3), args[4], paste0(ODIR, "/gwascat_snp2gene_MERGED.tsv")) #snp2gene.py
ifile_trait <-	ifelse((length(args)>4), args[5], paste0(ODIR, "/gwascat_trait.tsv")) #gwascat_trait.R
ifile_icite <-	ifelse((length(args)>5), args[6], paste0(ODIR, "/gwascat_icite.tsv")) #BioClients.icite API
ifile_ensembl <-ifelse((length(args)>6), args[7], paste0(ODIR, "/gwascat_EnsemblInfo.tsv")) #BioClients.ensembl API
ifile_tcrd <-	ifelse((length(args)>7), args[8], paste0(ODIR, "/tcrd_targets.tsv")) #BioClients.idg API
ofile <-	ifelse((length(args)>8), args[9], paste0(ODIR, "/gt_prepfilter.Rdata"))
ofile_filtered_studies <- ifelse((length(args)>9), args[10], paste0(ODIR, "/filtered_studies.tsv"))
ofile_filtered_traits <- ifelse((length(args)>10), args[11], paste0(ODIR, "/filtered_traits.tsv"))
ofile_filtered_genes <-	ifelse((length(args)>11), args[12], paste0(ODIR, "/filtered_genes.tsv"))
#
if (length(args)>12) {
  message("ERROR: Syntax: tiga_gt_prepfilter.R GWASFILE [COUNTSFILE [ASSNFILE [SNP2GENEFILE [TRAITFILE [ICITEFILE [TCRDFILE [ENSEMBLFILE [OFILE [OFILE_FILTERED_STUDIES [OFILE_FILTERED_TRAITS [OFILE_FILTERED_GENES]]]]]]]]]]]]\n...or... no args for defaults")
  quit()
}
message(sprintf("Input gwas file: %s", ifile_gwas))
message(sprintf("Input counts file: %s", ifile_counts))
message(sprintf("Input assn file: %s", ifile_assn))
message(sprintf("Input snp2gene file: %s", ifile_snp2gene))
message(sprintf("Input trait file: %s", ifile_trait))
message(sprintf("Input iCite file: %s", ifile_icite))
message(sprintf("Input TCRD file: %s", ifile_tcrd))
message(sprintf("Input Ensembl file: %s", ifile_ensembl))
message(sprintf("Output prepfilter file: %s", ofile))
message(sprintf("Output filtered studies file: %s", ofile_filtered_studies))
message(sprintf("Output filtered traits file: %s", ofile_filtered_traits))
message(sprintf("Output filtered genes file: %s", ofile_filtered_genes))
#
###
# Filter cutoffs:
PVAL_THRESHOLD <- 5e-8
PVAL_MLOG_THRESHOLD <- -log10(PVAL_THRESHOLD)
###
gwas <- read_delim(ifile_gwas, "\t", col_types=cols(.default=col_character(), DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"), ASSOCIATION_COUNT=col_integer(), study_N=col_integer()))
setDT(gwas)
gwas_counts <- read_delim(ifile_counts, "\t", col_types=cols(.default=col_number(), study_accession=col_character()))
setDT(gwas_counts)
message(sprintf("GWAS with mapped genes: %d", gwas_counts[gene_m_count>0, uniqueN(study_accession)]))
#
assn <- read_delim(ifile_assn, "\t", col_types=cols(.default=col_character(), 
	DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d"), 
	UPSTREAM_GENE_DISTANCE=col_integer(), DOWNSTREAM_GENE_DISTANCE=col_integer(), 
	oddsratio=col_double(), beta=col_double(), OR_or_BETA=col_double(), 
	PVALUE_MLOG=col_double(), `P-VALUE`=col_double()))
setDT(assn)
#
###
ensg_debug <- "ENSG00000042304"
assn_debug <- assn[grepl(ensg_debug, paste0(UPSTREAM_GENE_ID, DOWNSTREAM_GENE_ID, SNP_GENE_IDS)),]
message(sprintf("%s: assn count: %d: trait count: %d; SNPS count: %d; SNPS count (pval_mlog>%g): %d", ensg_debug, nrow(assn_debug), 
                assn_debug[, uniqueN(MAPPED_TRAIT)], assn_debug[, uniqueN(SNPS)], PVAL_MLOG_THRESHOLD, assn_debug[PVALUE_MLOG>=PVAL_MLOG_THRESHOLD, uniqueN(SNPS)]))
message(sprintf("%s: SNPS (pval_mlog>%g): %s", ensg_debug, PVAL_MLOG_THRESHOLD, paste(assn_debug[PVALUE_MLOG>=PVAL_MLOG_THRESHOLD, unique(SNPS)], collapse=", ")))
###
# Missing (UP|DOWN)STREAM_GENE_DISTANCE implies within MAPPED_GENE, thus distances zero for gene-distance weighting function, GDistWt.
assn[(!is.na(MAPPED_GENE) & is.na(UPSTREAM_GENE_DISTANCE) & is.na(DOWNSTREAM_GENE_DISTANCE)), `:=`(UPSTREAM_GENE_DISTANCE=0, DOWNSTREAM_GENE_DISTANCE=0)]
###
#
snp2gene <- read_delim(ifile_snp2gene, "\t", col_types=cols(.default=col_character(), MAPPED_OR_REPORTED=col_factor(c("r", "m", "md", "mu"))))
setDT(snp2gene)
message(sprintf("%s: snp2gene count: %d; snp count: %d", ensg_debug, nrow(snp2gene[ENSG==ensg_debug,]), snp2gene[ENSG==ensg_debug, uniqueN(SNP)]))
snp2gene <- unique(snp2gene[MAPPED_OR_REPORTED!="r"]) #Ignore reported.
message(sprintf("%s: snp2gene count: %d; snp count: %d", ensg_debug, nrow(snp2gene[ENSG==ensg_debug,]), snp2gene[ENSG==ensg_debug, uniqueN(SNP)]))
message(sprintf("%s: snps: %s", ensg_debug, paste(snp2gene[ENSG==ensg_debug, unique(SNP)], collapse=", ")))
#
ensemblInfo <- read_delim(ifile_ensembl, "\t", col_types = cols(.default=col_character(), version=col_integer(), strand=col_integer(), start=col_integer(), end=col_integer()))
setDT(ensemblInfo)
#
trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, old=c("MAPPED_TRAIT_URI", "MAPPED_TRAIT", "TRAIT"), new=c("TRAIT_URI", "TRAIT", "TRAIT_ORIG"))
trait[, TRAIT := iconv(TRAIT, from="latin1", to="UTF-8")]
#
###
# Estimate RCR for new publications as global median. Also any undefined RCR.
icite <- read_delim(ifile_icite, "\t", col_types=cols(.default=col_character(), relative_citation_ratio=col_double(), field_citation_rate=col_double(), citation_count=col_integer(), nih_percentile=col_double(), expected_citations_per_year=col_double(), citations_per_year=col_double(), year=col_integer()), na=c("", "NA", "None"))
setDT(icite)
rcr_median <- median(icite$relative_citation_ratio, na.rm=T) #global median ignoring NAs.
year_this <- as.integer(format(Sys.time(), "%Y"))
message(sprintf("Estimated RCR for new publications (and any undefined) (global median): %.3f", rcr_median))
#icite[is.na(relative_citation_ratio) & (year>=year_this-1) , relative_citation_ratio := rcr_median] #redundant
icite[is.na(relative_citation_ratio), relative_citation_ratio := rcr_median]
#
icite_gwas <- merge(icite[, .(pmid, relative_citation_ratio, year)], gwas[, .(PUBMEDID, STUDY_ACCESSION)], by.x="pmid", by.y="PUBMEDID", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, gwas_counts[, .(study_accession, trait_count, gene_r_count, gene_m_count)], by.x="STUDY_ACCESSION", by.y="study_accession", all.x=T, all.y=T)
icite_gwas <- merge(icite_gwas, icite_gwas[, .(study_perpmid_count = uniqueN(STUDY_ACCESSION)), by="pmid"], by="pmid")
message(sprintf("icite_gwas rows: %d", nrow(icite_gwas)))
message(sprintf("icite_gwas studies (STUDY_ACCESSION): %d", icite_gwas[, uniqueN(STUDY_ACCESSION)]))
###
# RCRAS = RCR-Aggregated-Score
# rcras_pmid is per-publication RCRAS normalized by study count in publication.
# rcras_study is per-study RCRAS, summing rcras_pmid, normalized by gene count in study.
# rcras_gt is per gene-trait association, computed in tiga_gt_variables.R, from rcras_study.
#
icite_gwas[, rcras_pmid := log2(relative_citation_ratio+1)/study_perpmid_count]
icite_gwas <- icite_gwas[gene_r_count>0 | gene_m_count>0] #Need genes to be useful
icite_gwas[gene_r_count==0, gene_r_count := NA]
icite_gwas[gene_m_count==0, gene_m_count := NA]
icite_gwas[, rcras_study := 1/gene_m_count * rcras_pmid] #Use mapped genes, not reported.
icite_gwas[is.na(rcras_study), rcras_study := 0]
icite_gwas <- icite_gwas[, .(pmid, STUDY_ACCESSION, year, relative_citation_ratio, rcras_pmid, rcras_study, trait_count, gene_r_count, gene_m_count, study_perpmid_count)]
setkey(icite_gwas, pmid, STUDY_ACCESSION) #Ensures unique study-pmid pairs each row.
message(sprintf("rcras_pmid: [%.3f, %.3f]; rcras_study: [%.3f, %.3f]", icite_gwas[, min(rcras_pmid)], icite_gwas[, max(rcras_pmid)], icite_gwas[, min(rcras_study)], icite_gwas[, max(rcras_study)]))
###
# Link to Ensembl via IDs from Catalog.
# (In 2020 (not 2018), EnsemblIDs available in catalog downloads, so API not required.)
setnames(ensemblInfo, old=c("id", "display_name"), new=c("ensemblId", "ensemblSymb"))
message(sprintf("%s: ensemblInfo symbols: %d", ensg_debug, ensemblInfo[ensemblId==ensg_debug, uniqueN(ensemblSymb)]))
###
# To avoid removing genes, rely on TCRD to ensure protein-coding.
# But, report how many genes would have been filtered.
# ensemblInfo annotations used: (1) ensemblSymb, (2) biotype (only logging, no impact on workflow).
#ensemblInfo <- ensemblInfo[biotype=="protein_coding" & description!="novel transcript"][, protein_coding := T]
###
snp2gene <- merge(snp2gene, ensemblInfo, by.x="ENSG", by.y="ensemblId", all.x=T, all.y=F, allow.cartesian=T)
message(sprintf("Genes NOT in Ensembl Info: %d", snp2gene[is.na(biotype), uniqueN(ENSG)]))
message(sprintf("Genes NOT defined protein_coding by Ensembl: %d / %d (%.1f%%)", 
                snp2gene[biotype != "protein_coding", uniqueN(ENSG)],
                snp2gene[, uniqueN(ENSG)],
                100 * snp2gene[biotype != "protein_coding", uniqueN(ENSG)] / snp2gene[, uniqueN(ENSG)]))
print(unique(snp2gene[, .(ENSG, biotype)])[, .(.N), by=biotype][order(-N)]) #DEBUG
snp2gene[, `:=`(GSYMB=NULL, MAPPED_OR_REPORTED=NULL)]
#snp2gene[, `:=`(protein_coding=NULL)]
setnames(snp2gene, "ENSG", "ensemblId")
snp2gene <- unique(snp2gene[, .(STUDY_ACCESSION, SNP, ensemblId, ensemblSymb)])
#
###
tcrd <- read_delim(ifile_tcrd, "\t", na=c("", "NA", "NULL"), col_types=cols(.default=col_character(), idgList=col_logical()))
setDT(tcrd)
###

###
# Counts:
message(sprintf("Association studies: %d", uniqueN(assn$STUDY_ACCESSION)))
#
###
# Reported genes ignored by TIGA.
#assn_reported <- assn[, .(STUDY_ACCESSION, `REPORTED_GENE(S)`)]
#assn_reported <- unique(assn_reported[, list(GENE=unlist(strsplit(`REPORTED_GENE(S)`, ", *"))), by=STUDY_ACCESSION])
###
ensgs_tcrd <- unique(tcrd$ensemblGeneId)
message(sprintf("TCRD targets: %d ; ensemblGeneIds: %d", nrow(tcrd), length(ensgs_tcrd)))
#
ensgs_tiga <- unique(snp2gene$ensemblId)
ensgs_common <- intersect(ensgs_tiga, ensgs_tcrd)
message(sprintf("EnsemblIds mapped to TCRD: %d / %d (%.1f%%)", length(ensgs_common), length(ensgs_tiga), 100 * length(ensgs_common) / length(ensgs_tiga)))
message(sprintf("%s: in TIGA?: %s; TCRD?: %s; common?: %s", ensg_debug, (ensg_debug %in% ensgs_tiga), (ensg_debug %in% ensgs_tcrd), (ensg_debug %in% ensgs_common)))
#
tcrd <- merge(tcrd, data.table(ensg=ensgs_tiga, in_gwascat=T), by.x="ensemblGeneId", by.y="ensg", all.x=T, all.y=F)
tcrd[, in_gwascat := !is.na(in_gwascat)]
message(sprintf("IDG-List targets mapped by GWAS: %d", tcrd[(idgList & in_gwascat), .N]))
#
### g2t should have one row for each gene-snp-study-trait association. Only 1 PUBMEDID per study.
g2t <- unique(snp2gene[, .(ensemblId, ensemblSymb, SNP, STUDY_ACCESSION)])
g2t <- merge(g2t, gwas[, .(STUDY_ACCESSION, PUBMEDID, study_N)], by="STUDY_ACCESSION", all.x=T, all.y=F, allow.cartesian=T)
g2t <- merge(g2t, assn[, .(SNPS, STUDY_ACCESSION, PVALUE_MLOG, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, oddsratio, beta)], 
             all.x=T, all.y=F, by.x=c("SNP", "STUDY_ACCESSION"), by.y=c("SNPS", "STUDY_ACCESSION"), allow.cartesian=T)
#
trait[, traitNstudy := fifelse(is.na(TRAIT_URI), as.integer(NA), uniqueN(STUDY_ACCESSION)), by="TRAIT_URI"]
g2t <- merge(g2t, trait, all.x=F, all.y=F, by="STUDY_ACCESSION", allow.cartesian=T)
message("DEBUG: merge 3 done.")
###
# For PC filter:
g2t <- merge(g2t, tcrd[!is.na(ensemblGeneId), c("ensemblGeneId", "TDL")], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F, allow.cartesian=T)
message("DEBUG: merge 4 done.")
#
debug_test <- function(ensg_test, g2t) {
  g2t_test <- g2t[ensemblId==ensg_test]
  message(sprintf("DEBUG: %s (%s) rows: %d", ensg_test, g2t_test$ensemblSymb[1], nrow(g2t_test)))
  print(g2t_test[, .(STUDY_ACCESSION, ensemblSymb, oddsratio, beta, PVALUE_MLOG, TRAIT, efoId)])
}
debug_test("ENSG00000170312", g2t)
###
# FILTERS:
###
# (0) PROTEIN-CODING filter: require mapping to TCRD/TDL.
reason_txt <- "Not protein-coding with TCRD TDL."
badrows_notpc <- g2t[, is.na(TDL)]
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_notpc, na.rm=T)))
filtered_studies <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_notpc]$STUDY_ACCESSION, g2t[!badrows_notpc]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F, allow.cartesian=T))
message("DEBUG: merge 5 done.")
filtered_studies[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies[, uniqueN(STUDY_ACCESSION)], filtered_studies[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
#
filtered_traits <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_notpc]$TRAIT_URI, g2t[!badrows_notpc]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F, allow.cartesian=T))
filtered_traits[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits[, uniqueN(TRAIT_URI)], filtered_traits[, uniqueN(TRAIT_URI)], 100*filtered_traits[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
#
filtered_genes <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_notpc]$ensemblId, g2t[!badrows_notpc]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F, allow.cartesian=T))
filtered_genes[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes[, uniqueN(ensemblId)], filtered_genes[, uniqueN(ensemblId)], 100*filtered_genes[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
###
# (1) TRAIT_URI filter: require MAPPED_TRAIT_URI (EFO)
reason_txt <- "Missing MAPPED_TRAIT_URI"
badrows_traituri <- (is.na(g2t$TRAIT_URI))
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_traituri, na.rm=T)))
filtered_studies <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_traituri]$STUDY_ACCESSION, g2t[!badrows_traituri]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F, allow.cartesian=T))
filtered_studies[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies[, uniqueN(STUDY_ACCESSION)], filtered_studies[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
#
filtered_traits <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_traituri]$TRAIT_URI, g2t[!badrows_traituri]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F, allow.cartesian=T))
filtered_traits[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits[, uniqueN(TRAIT_URI)], filtered_traits[, uniqueN(TRAIT_URI)], 100*filtered_traits[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
#
filtered_genes <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_traituri]$ensemblId, g2t[!badrows_traituri]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F, allow.cartesian=T))
filtered_genes[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes[, uniqueN(ensemblId)], filtered_genes[, uniqueN(ensemblId)], 100*filtered_genes[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
debug_test("ENSG00000170312", g2t[!badrows_traituri])
###
# (2) P-value filter: must exceed standard threshold for significance (pvalue <= 5e-8, pvalue_mlog >= 7.3).
# Usual genome-wide significance in GWAS, per LJJ 8/10/20 email.
message(sprintf("P-value significance threshold: %g; P-value minus-log threshold: %.3f)", PVAL_THRESHOLD, PVAL_MLOG_THRESHOLD))
reason_txt <- sprintf("P-value fails significance threshold (%g) (mlog=%.1f)", PVAL_THRESHOLD, PVAL_MLOG_THRESHOLD)
badrows_pval <- (g2t$PVALUE_MLOG<PVAL_MLOG_THRESHOLD)
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_pval, na.rm=T)))
filtered_studies_pval <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_pval]$STUDY_ACCESSION, g2t[!badrows_pval]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F, allow.cartesian=T))
filtered_studies_pval[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies_pval[, uniqueN(STUDY_ACCESSION)], filtered_studies_pval[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies_pval[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_pval))
#
filtered_traits_pval <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_pval]$TRAIT_URI, g2t[!badrows_pval]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F, allow.cartesian=T))
filtered_traits_pval[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits_pval[, uniqueN(TRAIT_URI)],
filtered_traits_pval[, uniqueN(TRAIT_URI)], 100*filtered_traits_pval[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
filtered_traits <- rbindlist(list(filtered_traits, filtered_traits_pval))
#
filtered_genes_pval <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_pval]$ensemblId, g2t[!badrows_pval]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F, allow.cartesian=T))
filtered_genes_pval[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes_pval[, uniqueN(ensemblId)],
filtered_genes_pval[, uniqueN(ensemblId)], 100*filtered_genes_pval[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
filtered_genes <- rbindlist(list(filtered_genes, filtered_genes_pval))
debug_test("ENSG00000170312", g2t[!badrows_pval])
###
# (3) Effect size filter: either OR or beta required.
badrows_effect <- (is.na(g2t$oddsratio) & is.na(g2t$beta))
reason_txt <- "Missing both OR and beta"
print(sprintf("%s: rows: %d", reason_txt, sum(badrows_effect, na.rm=T)))
###
filtered_studies_effect <- unique(merge(data.table(STUDY_ACCESSION = setdiff(g2t[badrows_effect]$STUDY_ACCESSION, g2t[!badrows_effect]$STUDY_ACCESSION)), gwas[, .(STUDY_ACCESSION, STUDY)], by="STUDY_ACCESSION", all.x=T, all.y=F, allow.cartesian=T))
filtered_studies_effect[, reason := reason_txt]
message(sprintf("Filtered studies (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(STUDY_ACCESSION)], g2t[, uniqueN(STUDY_ACCESSION)]-filtered_studies_effect[, uniqueN(STUDY_ACCESSION)], filtered_studies_effect[, uniqueN(STUDY_ACCESSION)], 100*filtered_studies_effect[, uniqueN(STUDY_ACCESSION)]/g2t[, uniqueN(STUDY_ACCESSION)]))
filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_effect))
#
filtered_traits_effect <- unique(merge(data.table(TRAIT_URI = setdiff(g2t[badrows_effect]$TRAIT_URI, g2t[!badrows_effect]$TRAIT_URI)), trait[, .(TRAIT_URI, TRAIT)], by="TRAIT_URI", all.x=T, all.y=F, allow.cartesian=T))
filtered_traits_effect[, reason := reason_txt]
message(sprintf("Filtered traits (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(TRAIT_URI)], g2t[, uniqueN(TRAIT_URI)]-filtered_traits_effect[, uniqueN(TRAIT_URI)],
filtered_traits_effect[, uniqueN(TRAIT_URI)], 100*filtered_traits_effect[, uniqueN(TRAIT_URI)]/g2t[, uniqueN(TRAIT_URI)]))
filtered_traits <- rbindlist(list(filtered_traits, filtered_traits_effect))
#
filtered_genes_effect <- unique(merge(data.table(ensemblId = setdiff(g2t[badrows_effect]$ensemblId, g2t[!badrows_effect]$ensemblId)), unique(g2t[, .(ensemblId, ensemblSymb)]), by="ensemblId", all.x=T, all.y=F, allow.cartesian=T))
filtered_genes_effect[, reason := reason_txt]
message(sprintf("Filtered genes (%s): %d -> %d (-%d; -%.1f%%)", reason_txt, g2t[, uniqueN(ensemblId)], g2t[, uniqueN(ensemblId)]-filtered_genes_effect[, uniqueN(ensemblId)],
filtered_genes_effect[, uniqueN(ensemblId)], 100*filtered_genes_effect[, uniqueN(ensemblId)]/g2t[, uniqueN(ensemblId)]))
filtered_genes <- rbindlist(list(filtered_genes, filtered_genes_effect))
#
###
# Informative: Any studies with BOTH OR and beta?
message(sprintf("Studies with BOTH OR and beta: %d", g2t[(!is.na(g2t$oddsratio) & !is.na(g2t$beta)), uniqueN(STUDY_ACCESSION)]))
###
filtered_studies <- filtered_studies[, reason := paste0(reason, collapse=" ; "), by=c("STUDY_ACCESSION", "STUDY")]
filtered_traits <- filtered_traits[, reason := paste0(reason, collapse=" ; "), by=c("TRAIT_URI", "TRAIT")]
filtered_genes <- filtered_genes[, reason := paste0(reason, collapse=" ; "), by=c("ensemblId", "ensemblSymb")]
#
filtered_genes <- merge(filtered_genes, tcrd[, .(ensemblGeneId, geneName=tcrdTargetName, geneFamily=tcrdTargetFamily, TDL)], by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F, allow.cartesian=T)
filtered_genes <- filtered_genes[, .(ensemblId, ensemblSymb, geneName, geneFamily, TDL, reason)]
message(sprintf("Filtered: Gene (ensemblId) count: %d", filtered_genes[, uniqueN(ensemblId)]))
message(sprintf("Filtered: Missing gene symbols: %d", filtered_genes[is.na(ensemblSymb), uniqueN(ensemblId)]))
#
# Write files accounting for filtered studies, traits and genes.
write_delim(filtered_studies, ofile_filtered_studies, delim="\t")
write_delim(filtered_traits, ofile_filtered_traits, delim="\t")
write_delim(filtered_genes, ofile_filtered_genes, delim="\t")
#
badrows <- (badrows_notpc | badrows_traituri | badrows_pval | badrows_effect)
message(sprintf("Filtered associations (total): %d -> %d (-%d; -%.1f%%)", nrow(g2t), nrow(g2t)-sum(badrows), sum(badrows), 100*sum(badrows)/nrow(g2t)))
g2t <- g2t[!badrows] #Many filtered.
#
###
#
message(sprintf("GT: Final nrow(g2t) = %d", nrow(g2t)))
message(sprintf("GT: associations in dataset: %d", nrow(unique(g2t[, .(ensemblId, efoId)]))))
message(sprintf("GT: Study (STUDY_ACCESSION) count: %d", g2t[, uniqueN(STUDY_ACCESSION)]))
message(sprintf("GT: Gene (ensemblId) count: %d", g2t[, uniqueN(ensemblId)]))
message(sprintf("GT: Missing gene symbols: %d", g2t[is.na(ensemblSymb), uniqueN(ensemblId)])) #needed for autocomplete
message(sprintf("GT: Trait (efoId) count: %d", g2t[, uniqueN(efoId)]))
message(sprintf("GT: PVALUE_MLOG (instance) count: %d", nrow(g2t[!is.na(PVALUE_MLOG)])))
message(sprintf("GT: oddsratio (instance) count: %d", nrow(g2t[!is.na(oddsratio)])))
message(sprintf("GT: beta (instance) count: %d", nrow(g2t[!is.na(beta)])))
#
###
# Save to file.
save(g2t, tcrd, icite_gwas, file=ofile)
writeLines(sprintf("Output Rdata file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))

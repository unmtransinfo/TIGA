#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT VARIABLES
###   (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
###  (2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
### *(2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
###   (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR ALL INPUT WORKFLOWS: Go_TIGA_Workflow.sh
#############################################################################
# Now using mean-rank instead of mu_scores.
# Based on DISEASES benchmark results, using only variables:
#   * n_study
#   * pvalue_mlog_median (MAYBE CHANGING TO pvalue_mlog_max)
#   * rcras
#############################################################################
# TCRD-TDL-MAPPING/PROTEIN-CODING FILTERING MOVED TO tiga_gt_prepfilter.R
# FOR ACCOUNTING AND REDUCED COMPUTATION.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
#
#
if (interactive()) {
  rel_y <- as.integer(readline(prompt="Enter RELEASE_YEAR: "))
  rel_m <- as.integer(readline(prompt="Enter RELEASE_MONTH: "))
  rel_d <- as.integer(readline(prompt="Enter RELEASE_DAY: "))
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (length(args)==3) {
  rel_y <- as.integer(args[1])
  rel_m <- as.integer(args[2])
  rel_d <- as.integer(args[3])
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (file.exists("LATEST_RELEASE_GWC.txt")) {
  GC_REL <- trimws(read_file("LATEST_RELEASE_GWC.txt"))
  rel_y <- as.integer(sub("\\-.*$", "", GC_REL))
  rel_m <- as.integer(sub("\\d+\\-(\\d+)\\-.*$", "\\1", GC_REL))
  rel_d <- as.integer(sub("\\d+\\-\\d+\\-(\\d+).*$", "\\1", GC_REL))
  message(sprintf("LATEST_RELEASE_GWC: %s", GC_REL))
  ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))
} else {
  message("ERROR: Syntax: tiga_gt_variables.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
#
ifile	<- paste0(ODIR, "/gt_prepfilter.Rdata")
ofile	<- paste0(ODIR, "/gt_variables.tsv.gz")
#
message(sprintf("Input prepfilter file: %s", ifile))
message(sprintf("Output variables file: %s", ofile))
#
###
load(ifile)
#
###
# Check file contents.
message(sprintf("g2t rows: %d", nrow(g2t)))
message(sprintf("g2t genes (ensemblId): %d", g2t[, uniqueN(ensemblId)]))
message(sprintf("g2t traits (TRAIT_URI): %d", g2t[, uniqueN(TRAIT_URI)]))
message(sprintf("g2t gene-trait pairs: %d", nrow(unique(g2t[, .(ensemblId, TRAIT_URI)]))))
message(sprintf("g2t studies (STUDY_ACCESSION): %d", g2t[, uniqueN(STUDY_ACCESSION)]))
message(sprintf("g2t papers (PUBMEDID): %d", g2t[, uniqueN(PUBMEDID)]))
message(sprintf("g2t study_N (instances): %d (%.1f%%)", nrow(g2t[!is.na(study_N)]),
                100*nrow(g2t[!is.na(study_N)]) /nrow(g2t)))
message(sprintf("g2t traitNstudy (instances): %d (%.1f%%)", nrow(g2t[!is.na(traitNstudy)]),
                100*nrow(g2t[!is.na(traitNstudy)]) /nrow(g2t)))
message(sprintf("g2t PVALUE_MLOG (instances): %d (%.1f%%)", nrow(g2t[!is.na(PVALUE_MLOG)]),
                100*nrow(g2t[!is.na(PVALUE_MLOG)]) /nrow(g2t)))
message(sprintf("g2t oddsratio (instances): %d (%.1f%%)", nrow(g2t[!is.na(oddsratio)]),
                100*nrow(g2t[!is.na(oddsratio)]) /nrow(g2t)))
message(sprintf("g2t beta (instances): %d (%.1f%%)", nrow(g2t[!is.na(beta)]),
                100*nrow(g2t[!is.na(beta)]) /nrow(g2t)))
#
message(sprintf("icite_gwas rows: %d", nrow(icite_gwas)))
message(sprintf("icite_gwas studies (STUDY_ACCESSION): %d", icite_gwas[, uniqueN(STUDY_ACCESSION)]))
message(sprintf("icite_gwas relative_citation_ratio (instances): %d (%.1f%%)", nrow(icite_gwas[!is.na(relative_citation_ratio)]),
                100*nrow(icite_gwas[!is.na(relative_citation_ratio)]) /nrow(icite_gwas)))
message(sprintf("icite_gwas study_perpmid_count (instances): %d (%.1f%%)", nrow(icite_gwas[!is.na(study_perpmid_count)]),
                100*nrow(icite_gwas[!is.na(study_perpmid_count)]) /nrow(icite_gwas)))
message(sprintf("icite_gwas rcras_pmid (instances): %d (%.1f%%)", nrow(icite_gwas[!is.na(rcras_pmid)]),
                100*nrow(icite_gwas[!is.na(rcras_pmid)]) /nrow(icite_gwas)))
message(sprintf("icite_gwas rcras_study (instances): %d (%.1f%%)", nrow(icite_gwas[!is.na(rcras_study)]),
                100*nrow(icite_gwas[!is.na(rcras_study)]) /nrow(icite_gwas)))
message(sprintf("icite_gwas gene_m_count (instances): %d (%.1f%%)", nrow(icite_gwas[!is.na(gene_m_count)]),
                100*nrow(icite_gwas[!is.na(gene_m_count)]) /nrow(icite_gwas)))
#
message(sprintf("tcrd rows: %d", nrow(tcrd)))
###
# Gene-distance weighting function.
g2t[, GDistWt := 2^(-pmin(g2t$UPSTREAM_GENE_DISTANCE, g2t$DOWNSTREAM_GENE_DISTANCE, na.rm=T)/5e4)]
message(sprintf("DEBUG: GDistWt: count: %d / %d (%.1f%%)", sum(!is.na(g2t$GDistWt)), nrow(g2t), 100*sum(!is.na(g2t$GDistWt))/nrow(g2t)))
#
###
### GENE-TRAIT stats
### One row per unique gene-trait pair.
### From g2t, create gt_stats table for TSV export.
### Slow. Vectorize/optimize!?
#
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
	pvalue_mlog_max=as.numeric(rep(NA, NROW)),
	or_median=as.numeric(rep(NA, NROW)),
	n_beta=as.integer(rep(NA, NROW)), #simple count of beta values
	study_N_mean=as.numeric(rep(NA, NROW)),
	rcras=rep(NA, NROW)
	)
#
message(sprintf("Initialized rows to be populated: nrow(gt_stats) = %d", nrow(gt_stats)))
#
i_gt <- 0 #gt_stats populated row count (gene-trait pairs)
# gene-loop:
for (ensg in unique(g2t$ensemblId)) {
  geneNstudy <- g2t[ensemblId==ensg, uniqueN(STUDY_ACCESSION)]
  # trait-loop:
  for (trait_uri in unique(g2t[ensemblId==ensg, TRAIT_URI])) {
    i_gt <- i_gt + 1
    if ((i_gt%%10000)==0) {
      t_elapsed <- (Sys.time()-t_start)
      message(sprintf("%d / %d (%.1f%%) GTs; %s, elapsed: %.1f %s", i_gt, NROW, 100*i_gt/NROW, Sys.time(), t_elapsed, attr(t_elapsed, "units")))
    }
    gt_stats$ensemblId[i_gt] <- ensg
    gt_stats$efoId[i_gt] <- sub("^.*/", "", trait_uri)
    gt_stats$geneNstudy[i_gt] <- geneNstudy
    gt_stats$trait[i_gt] <- g2t[ensemblId==ensg & TRAIT_URI==trait_uri, TRAIT][1]
    gt_stats$traitNstudy[i_gt] <- g2t[TRAIT_URI==trait_uri, traitNstudy][1]
    gt_stats$n_study[i_gt] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])
    gt_stats$pvalue_mlog_median[i_gt] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, PVALUE_MLOG], na.rm=T)
    gt_stats$pvalue_mlog_max[i_gt] <- max(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, PVALUE_MLOG], na.rm=T)
    gt_stats$or_median[i_gt] <- median(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, oddsratio], na.rm=T) #NA if no ORs
    gt_stats$n_beta[i_gt] <- g2t[ensemblId==ensg & TRAIT_URI==trait_uri & !is.na(beta), .N] #0 if no betas
    gt_stats$study_N_mean[i_gt] <- round(mean(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, study_N], na.rm=T), 1)
    gt_stats$n_snp[i_gt] <- uniqueN(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, SNP])
    # Deduplicate (group-by) SNPs for `n_snpw` computation. 
    gt_stats$n_snpw[i_gt] <- sum(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, .(GDistWt = median(GDistWt, na.rm=T)), by="SNP"][, GDistWt], na.rm=T)
    #
    staccs <- unique(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])
    gt_stats$rcras[i_gt] <- icite_gwas[STUDY_ACCESSION %in% staccs, sum(rcras_study)]
  }
}
message(sprintf("DEBUG: nrow(gt_stats): %d; final i_gt: %d", nrow(gt_stats), i_gt))
gt_stats[, traitNgene := .N, by="efoId"]
gt_stats[, geneNtrait := .N, by="ensemblId"]
#
gt_stats$or_median <- round(as.double(gt_stats$or_median), 3)
gt_stats$pvalue_mlog_median <- round(as.double(gt_stats$pvalue_mlog_median), 3)
gt_stats$pvalue_mlog_max <- round(as.double(gt_stats$pvalue_mlog_max), 3)
gt_stats$rcras <- round(as.double(gt_stats$rcras), 3)
gt_stats$n_snpw <- round(as.double(gt_stats$n_snpw), 3)
#
message(sprintf("%s, elapsed time: %.2f %s", Sys.time(), t_elapsed, attr(t_elapsed, "units")))
message(sprintf("Final: nrow(gt_stats) = %d", nrow(gt_stats)))
message(sprintf("gene (ensemblId) count: %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("trait (efoId) count: %d", uniqueN(gt_stats$efoId)))
message(sprintf("traitNgene: [%d,%d]", min(gt_stats$traitNgene), max(gt_stats$traitNgene)))
message(sprintf("traitNstudy: [%d,%d]", min(gt_stats$traitNstudy, na.rm=T), max(gt_stats$traitNstudy, na.rm=T)))
message(sprintf("geneNtrait: [%d,%d]", min(gt_stats$geneNtrait), max(gt_stats$geneNtrait)))
message(sprintf("pvalue_mlog_median: [%.2f,%.2f]", min(gt_stats$pvalue_mlog_median, na.rm=T), max(gt_stats$pvalue_mlog_median, na.rm=T)))
message(sprintf("pvalue_mlog_max: [%.2f,%.2f]", min(gt_stats$pvalue_mlog_max, na.rm=T), max(gt_stats$pvalue_mlog_max, na.rm=T)))
message(sprintf("or_median: [%.2f,%.2f]", min(gt_stats$or_median, na.rm=T), max(gt_stats$or_median, na.rm=T)))
message(sprintf("n_beta: [%d,%d]", min(gt_stats$n_beta, na.rm=T), max(gt_stats$n_beta, na.rm=T)))
message(sprintf("study_N_mean: [%.1f,%.1f]", min(gt_stats$study_N_mean, na.rm=T), max(gt_stats$study_N_mean, na.rm=T)))
message(sprintf("rcras: [%.2f,%.2f]", min(gt_stats$rcras, na.rm=T), max(gt_stats$rcras, na.rm=T)))
message(sprintf("n_snpw: [%.2f,%.2f]", min(gt_stats$n_snpw, na.rm=T), max(gt_stats$n_snpw, na.rm=T)))
#
message(sprintf("Rows missing ensemblId: %d", nrow(gt_stats[is.na(ensemblId)])))
message(sprintf("Rows missing efoId: %d", nrow(gt_stats[is.na(efoId)])))
message(sprintf("Rows missing ensemblId or efoId: %d", nrow(gt_stats[is.na(ensemblId)|is.na(efoId)])))
gt_stats <- gt_stats[!(is.na(ensemblId)|is.na(efoId))] #Should be no-op.
#
gt_stats <- merge(gt_stats, tcrd[!is.na(ensemblGeneId), c("ensemblGeneId", "tcrdGeneSymbol", "TDL", "tcrdTargetFamily", "idgList", "tcrdTargetName")], 
                  by.x="ensemblId", by.y="ensemblGeneId", all.x=T, all.y=F, allow.cartesian=F)
setnames(gt_stats,
	old=c("tcrdGeneSymbol", "tcrdTargetName", "tcrdTargetFamily", "TDL", "idgList"),
	new=c("geneSymbol", "geneName", "geneFamily", "geneIdgTdl", "geneIdgList"))
#
message(sprintf("Genes NOT mapped to TCRD/TDL (thus not protein-coding): %d", gt_stats[is.na(geneIdgTdl), uniqueN(ensemblId)]))
#Non-protein-coding removed by IDG TDL requirement (MOVED TO PREPFILTER).
#gt_stats <- gt_stats[!is.na(geneIdgTdl)]
###
#
message(sprintf("TOTAL Genes (ensemblIDs): %d", gt_stats[, uniqueN(ensemblId)]))
message(sprintf("TOTAL Genes (symbols): %d", gt_stats[, uniqueN(geneSymbol)]))
message(sprintf("TOTAL Traits in dataset: %d", gt_stats[, uniqueN(efoId)]))
message(sprintf("TOTAL GT associations in dataset: %d", nrow(gt_stats)))
#
message(sprintf("PROTEIN-CODING Genes (ensemblIDs): %d", gt_stats[!is.na(geneIdgTdl), uniqueN(ensemblId)]))
message(sprintf("PROTEIN-CODING Genes (symbols): %d", gt_stats[!is.na(geneIdgTdl), uniqueN(geneSymbol)]))
message(sprintf("PROTEIN-CODING Traits in dataset: %d", gt_stats[!is.na(geneIdgTdl), uniqueN(efoId)]))
message(sprintf("PROTEIN-CODING GT associations in dataset: %d", nrow(gt_stats[!is.na(geneIdgTdl)])))
###
# Save computed variables to file.
write_delim(gt_stats, ofile, delim="\t")
message(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2f %s", t_elapsed, attr(t_elapsed, "units")))

#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT VARIABLES
###   (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
###  (2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
### *(2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
###   (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR ALL INPUT WORKFLOWS: Go_gwascat_GetData.sh
#############################################################################
# Now using mean-rank instead of mu_scores.
# Based on DISEASES benchmark results, using only variables:
#   * n_study
#   * pvalue_mlog_median
#   * rcras
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t0 <- proc.time()
t_start <- Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==2) {
  (ifile	<- args[1])
  (ofile	<- args[2])
} else if (length(args)==0) {
  ifile <- "data/gt_prepfilter.Rdata"
  ofile <- "data/gt_variables.tsv.gz"
} else {
  message("ERROR: Syntax: tiga_gt_variables.R GT_PREPFILTER_FILE OFILE\n...or... no args for defaults")
  quit()
}
writeLines(sprintf("Input prepfilter file: %s", ifile))
writeLines(sprintf("Output variables file: %s", ofile))
#
###
load(ifile)
#
###
# Gene-distance weighting function.
g2t[, GDistWt := 2^(-pmin(g2t$UPSTREAM_GENE_DISTANCE, g2t$DOWNSTREAM_GENE_DISTANCE, na.rm=T)/5e4)]
message(sprintf("DEBUG: GDistWt: count: %d / %d (%.1f%%)", 
                sum(!is.na(g2t$GDistWt)), nrow(g2t), 100*sum(!is.na(g2t$GDistWt))/nrow(g2t)))
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
	or_median=as.numeric(rep(NA, NROW)),
	n_beta=as.integer(rep(NA, NROW)), #simple count of beta values
	study_N_mean=as.numeric(rep(NA, NROW)),
	rcras=rep(NA, NROW)
	)
#
message(sprintf("Initialized rows to be populated: nrow(gt_stats) = %d", nrow(gt_stats)))
#
i_row <- 0 #gt_stats populated row count
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
    rcras <- 0
    for (stacc in unique(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, STUDY_ACCESSION])) {
      rcras_study <- 0
      for (pmid_this in icite_gwas[STUDY_ACCESSION==stacc, pmid]) {
        if (nrow(icite_gwas[STUDY_ACCESSION==stacc & pmid==pmid_this])==0) { next }
        spp <- icite_gwas[STUDY_ACCESSION==stacc & pmid==pmid_this]$study_perpmid_count[1]
        rcr <- icite_gwas[STUDY_ACCESSION==stacc & pmid==pmid_this]$relative_citation_ratio[1]
        if (is.na(rcr) | spp==0) { next; }
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
###
#
message(sprintf("Genes (ensemblIDs): %d", uniqueN(gt_stats$ensemblId)))
message(sprintf("Genes (symbols): %d", uniqueN(gt_stats$geneSymbol)))
message(sprintf("Traits in dataset: %d", uniqueN(gt_stats$efoId)))
message(sprintf("G-T associations in dataset: %d", nrow(gt_stats)))
###

# Save computed variables to file.
write_delim(gt_stats, ofile, delim="\t")
writeLines(sprintf("Output file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2fs (%.2f %s)", (proc.time()-t0)[3], t_elapsed, attr(t_elapsed, "units")))

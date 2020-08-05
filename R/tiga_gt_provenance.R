#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT PROVENANCE
###   (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
### *(2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
###  (2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
###   (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR ALL INPUT WORKFLOWS: Go_gwascat_GetData.sh
#############################################################################
# Writes gt_provenance file with TRAIT_URI, ensemblId, STUDY_ACCESSION and PUBMEDID.
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
  ofile <- "data/gt_provenance.tsv.gz"
} else {
  message("ERROR: Syntax: tiga_gt_provenance.R GT_PREPFILTER_FILE OFILE\n...or... no args for defaults")
  quit()
}
writeLines(sprintf("Input prepfilter file: %s", ifile_gt))
writeLines(sprintf("Output provenance file: %s", ofile))
#
###
load(ifile)
#
###
# GENE-TRAIT provenance
message("Generating provenance file for each gene-trait pair.")
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
#
###
message(sprintf("Genes (ensemblIDs): %d", uniqueN(gt_prov$ensemblId)))
message(sprintf("Traits (efoIds): %d", uniqueN(gt_prov$efoId)))
message(sprintf("G-T associations in dataset: %d", nrow(unique(gt_prov[, .(ensemblId, efoId)]))))
#
# Save to file.
write_delim(gt_prov, ofile, delim="\t")
writeLines(sprintf("Output provenance file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2fs (%.2f %s)", (proc.time()-t0)[3], t_elapsed, attr(t_elapsed, "units")))

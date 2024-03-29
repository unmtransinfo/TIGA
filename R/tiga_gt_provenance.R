#!/usr/bin/env Rscript
#############################################################################
### GENE-TRAIT PROVENANCE
###   (1) tiga_gt_prepfilter.R - Merge input files, preprocess and filter.
### *(2a) tiga_gt_provenance.R - Produce gt_provenance.tsv.gz, for TIGA app.
###  (2b) tiga_gt_variables.R - Produce gt_variables.tsv.gz
###   (3) tiga_gt_stats.R, to produce gt_stats.tsv.gz, for TIGA app.
#############################################################################
### SEE FOR ALL INPUT WORKFLOWS: Go_TIGA_Workflow.sh
#############################################################################
# Writes gt_provenance file with TRAIT_URI, ensemblId, STUDY_ACCESSION and PUBMEDID.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
t_start <- Sys.time()
#
message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
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
  message("ERROR: Syntax: tiga_gt_provenance.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
#
ifile	<- paste0(ODIR, "/gt_prepfilter.Rdata")
ofile	<- paste0(ODIR, "/gt_provenance.tsv.gz")
#
writeLines(sprintf("Input prepfilter file: %s", ifile))
writeLines(sprintf("Output provenance file: %s", ofile))
#
###
load(ifile)
#
n_gt_pairs <- nrow(unique(g2t[(!is.na(STUDY_ACCESSION) & !is.na(PUBMEDID)), .(ensemblId, TRAIT_URI)]))
message(sprintf("Input genes: %d", uniqueN(g2t$ensemblId)))
message(sprintf("Input traits: %d", uniqueN(g2t$TRAIT_URI)))
message(sprintf("Input GTs (gene-trait pairs) with provenance: %d", n_gt_pairs))
###
# GENE-TRAIT provenance
message("Enumerating study & publication provenance for each gene-trait pair.")
gt_prov <- NULL
i <- 0
# gene-loop:
for (ensg in unique(g2t$ensemblId)) {
  geneNstudy <- g2t[ensemblId==ensg, uniqueN(STUDY_ACCESSION)]
  # trait-loop:
  for (trait_uri in unique(g2t[ensemblId==ensg, TRAIT_URI])) {
    i <- i + 1
    studies_this <- unique(g2t[ensemblId==ensg & TRAIT_URI==trait_uri, .(STUDY_ACCESSION, PUBMEDID)])
    if ((i%%10000)==0) {
      t_elapsed <- (Sys.time()-t_start)
      message(sprintf("%d / %d (%.1f%%) GTs; %s, elapsed: %.1f %s", i, n_gt_pairs, 100*i/n_gt_pairs, Sys.time(), t_elapsed, attr(t_elapsed, "units")))
    }
    gt_prov_this <- data.table(ensemblId=rep(ensg, nrow(studies_this)), TRAIT_URI=rep(trait_uri, nrow(studies_this)), STUDY_ACCESSION=studies_this[, STUDY_ACCESSION], PUBMEDID=studies_this[, PUBMEDID])
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
message(sprintf("Provenance genes (ensemblIDs): %d", uniqueN(gt_prov$ensemblId)))
message(sprintf("Provenance traits (efoIds): %d", uniqueN(gt_prov$efoId)))
message(sprintf("Provenance GT associations: %d", nrow(unique(gt_prov[, .(ensemblId, efoId)]))))
#
# Save to file.
write_delim(gt_prov, ofile, delim="\t")
writeLines(sprintf("Output provenance file written: %s", ofile))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2f %s", t_elapsed, attr(t_elapsed, "units")))

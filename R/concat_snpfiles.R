#!/usr/bin/env Rscript
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
#
message(paste(commandArgs(), collapse=" "))
#
ifile_snps_a <- "data/gwascat_Snps_a.tsv"
ifile_snps_b <- "data/gwascat_Snps_skip125000.tsv"
ofile <- "data/gwascat_Snps.tsv"
message(sprintf("ifile_snps_a: %s", ifile_snps_a))
message(sprintf("ifile_snps_b: %s", ifile_snps_b))
message(sprintf("ofile: %s", ofile))
#
###
studySnps_a <- read_delim(ifile_snps_a, "\t", col_types=cols(.default=col_character(), merged=col_logical(), lastUpdateDate=col_datetime(),   isIntergenic = col_logical(), isUpstream = col_logical(), isDownstream = col_logical(), distance = col_double(), isClosestGene = col_logical(), chromosomePosition = col_double()))
setDT(studySnps_a)

studySnps_b <- read_delim(ifile_snps_b, "\t", col_types=cols(.default=col_character(), merged=col_logical(), lastUpdateDate=col_datetime(),   isIntergenic = col_logical(), isUpstream = col_logical(), isDownstream = col_logical(), distance = col_double(), isClosestGene = col_logical(), chromosomePosition = col_double()))
setDT(studySnps_b)

message(sprintf("nrow(studySnps_a): %d", nrow(studySnps_a)))
message(sprintf("nrow(studySnps_b): %d", nrow(studySnps_b)))
studySnps <- rbindlist(list(studySnps_a, studySnps_b))
message(sprintf("nrow(studySnps): %d", nrow(studySnps)))
studySnps <- unique(studySnps)
message(sprintf("nrow(studySnps): %d UNIQUE", nrow(studySnps)))
#
#
write_delim(studySnps, ofile, delim="\t")
#

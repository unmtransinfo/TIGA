#!/usr/bin/env Rscript
##########################################################################################
### Create efo_sub_gwas.tsv with subclasses for GWAS/EFO traits.
##########################################################################################
library(readr)
library(data.table, quietly=T)

message(paste(commandArgs(), collapse=" "))

if (length(commandArgs(trailingOnly=T))>0) {
  args <- commandArgs(trailingOnly=T)
}

if (length(args)==3) {
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
  message("ERROR: Syntax: gwascat_trait_efokg.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}

ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
#
ifile_trait <- paste0(ODIR, "/gwascat_trait.tsv")
ifile_efo <- paste0(ODIR, "/efo.tsv")
ofile_subclass <- paste0(ODIR, "/efo_sub_gwas.tsv.gz") #Large, thus gzip.

message(sprintf("Input traits: %s", ifile_trait))
message(sprintf("Input EFO: %s", ifile_efo))
message(sprintf("Output subclass: %s", ofile_subclass))

###
# Trait file produced by gwascat_trait.R
#
trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
###
# EFO (full ontology)
#
efo <- read_delim(ifile_efo, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_node <- efo[node_or_edge == "node"]
efo_node[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]
efo_node[['EFO_prefix']] <- as.factor(sub("_.*$", "", efo_node$id))
efo_counts <- efo_node[, .(N_class = .N), by="EFO_prefix"][order(-N_class)]
print(efo_counts[1:10])
message(sprintf("Other prefix (total=%d): %d", uniqueN(efo_counts$EFO_prefix), efo_counts[11:nrow(efo_counts), sum(N_class)]))
efo_node[, EFO_prefix := NULL]


###
# Subclass file, of study pairs with subclass-related traits:
# The file produced is huge (112G for 20260720 release), and not used for the TIGA app.
# Memory use can crash R. Thus, maybe should be skipped.
#
efo_sub <- efo[node_or_edge == "edge" & label=="has_subclass"]
efo_sub[, `:=`(node_or_edge = NULL, id = NULL, label = NULL, uri = NULL, comment = NULL)]
setnames(efo_sub, old=c("source", "target"), new=c("trait_uri", "subclass_uri"))
efo_sub <- efo_sub[trait_uri %in% trait$MAPPED_TRAIT_URI & subclass_uri %in% trait$MAPPED_TRAIT_URI]
efo_sub <- merge(efo_sub, efo_node[, .(trait_uri = uri, trait = label)], by="trait_uri")
efo_sub <- merge(efo_sub, efo_node[, .(subclass_uri = uri, subclass_trait = label)], by="subclass_uri")
efo_sub <- merge(efo_sub, unique(trait[, .(study_accession = STUDY_ACCESSION, trait_uri = MAPPED_TRAIT_URI)]), by="trait_uri", allow.cartesian = T)
efo_sub <- merge(efo_sub, unique(trait[, .(study_accession_subclass = STUDY_ACCESSION, subclass_uri = MAPPED_TRAIT_URI)]), by="subclass_uri", allow.cartesian = T)

efo_sub <- unique(efo_sub[, .(study_accession, trait, trait_uri, study_accession_subclass, subclass_trait, subclass_uri)])

#Example
unique(efo_sub[grepl("(mood|bipolar)", trait), .(efoId = sub("^.*/", "", trait_uri), subclass_efoId = sub("^.*/", "", subclass_uri), N_study = uniqueN(study_accession), N_study_subclass = uniqueN(study_accession_subclass)), by=c("trait", "subclass_trait")])
#
write_delim(efo_sub, ofile_subclass, "\t")


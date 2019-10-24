#!/usr/bin/env Rscript
##########################################################################################
### Analyze/describe assn file traits.
### Note that EFO includes entities and IDs from other ontologies.
### Output gwascat_trait.tsv for use by gwax_gt_stats.R, but only for mapping from
### STUDY_ACCESSION to MAPPED_TRAIT_URI.
##########################################################################################
library(readr)
require(data.table, quietly=T)

ifile_default <- "/home/data/gwascatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv"
ofile_default <- "data/gwascat_trait.tsv"

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  (ifile <- args[1])
} else if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- ifile_default
  ofile <- ofile_default
} else {
  message("ERROR: Syntax: gwascat_trait.R [GWASFILE [OFILE]]\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

ifile_efo <- "data/efo.tsv"

trait <- read_delim(ifile, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, old="STUDY ACCESSION", new="STUDY_ACCESSION")
trait <- trait[, .(STUDY_ACCESSION, MAPPED_TRAIT, MAPPED_TRAIT_URI)]

trait <- unique(trait[complete.cases(trait),])

# Split comma separated vals.
trait_multi <- trait[grepl(",", trait$MAPPED_TRAIT_URI)]
trait <- trait[!grepl(",", trait$MAPPED_TRAIT_URI)]
for (i in 1:nrow(trait_multi)) {
  uris <- strsplit(trait_multi$MAPPED_TRAIT_URI[i], ', ', perl=T)[[1]]
  traits <- strsplit(trait_multi$MAPPED_TRAIT[i], ', ', perl=T)[[1]]
  accs <- rep(trait_multi$STUDY_ACCESSION[i], length(uris))
  if (length(uris)!=length(traits)) {
    message(sprintf("ERROR: length(uris)!=length(traits) (%d!=%d) \"%s\"", length(uris), length(traits), trait_multi$MAPPED_TRAIT[i]))
    traits <- rep(trait_multi$MAPPED_TRAIT[i], length(uris)) #Commas in trait names, so must be curated manually.
  }
  trait <- rbind(trait, data.frame(STUDY_ACCESSION=accs, MAPPED_TRAIT=traits, MAPPED_TRAIT_URI=uris))
}

trait[['id']] <- as.factor(sub("^.*/","", trait$MAPPED_TRAIT_URI))
trait[['ontology']] <- as.factor(sub("_.*$","", trait$id))
trait <- unique(trait)

###
#v1.0.2 counts:
#EFO: 1794
#GO: 69
#HP: 49
#Orphanet: 43
#CHEBI: 2
#PATO: 1
###

trait_ont <- trait[, .N, by="ontology"][order(-N)]
message(sprintf("%12s: %4d / %4d (%4.1f%%)\n", trait_ont$ontology, trait_ont$N, sum(trait_ont$N), 100*trait_ont$N/sum(trait_ont$N)))

trait_study_counts <- trait[, .(N_study = uniqueN(STUDY_ACCESSION)), by="MAPPED_TRAIT_URI"][order(-N_study)]
for (i in 1:max(trait_study_counts$N_study)) {
  if (i %in% trait_study_counts$N_study) {
    message(sprintf("N_study=%d: N_trait=%d", i, trait_study_counts[N_study==i, uniqueN(MAPPED_TRAIT_URI)]))
  }
}

###
efo <- read_delim(ifile_efo, "\t")
setDT(efo)
efo <- efo[node_or_edge == "node"]
efo[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]

# EFO includes other ontologies.
efo[['ontology']] <- as.factor(sub("_.*$","", efo$id))
efo_counts <- efo[, .N, by="ontology"][order(-N)]
message(sprintf("%12s: %4d / %4d (%4.1f%%)\n", efo_counts$ontology, efo_counts$N, sum(efo_counts$N), 100*efo_counts$N/sum(efo_counts$N)))
efo[, ontology := NULL]

trait <- merge(trait, efo[, .(id, efo_label = label)], by="id", all.x=T, all.y=F)

trait_unmapped <- trait[is.na(efo_label), .(id, MAPPED_TRAIT)]
n_trait_mapped <- nrow(trait[!is.na(efo_label)])
message(sprintf("EFO+ IDs mapped to efo.owl: %d / %d (%.1f%%)", n_trait_mapped, 
                n_trait_mapped + uniqueN(trait_unmapped$id),
                100 * (n_trait_mapped / (n_trait_mapped + uniqueN(trait_unmapped$id)))))
#
write_delim(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT, MAPPED_TRAIT_URI, id, efo_label)], ofile, "\t")

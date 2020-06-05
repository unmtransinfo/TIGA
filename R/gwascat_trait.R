#!/usr/bin/env Rscript
##########################################################################################
### Analyze/describe assn file traits.
### Note that EFO includes entities and IDs from other ontologies.
### Output gwascat_trait.tsv for use by tiga_gt_stats.R, but only for mapping from
### STUDY_ACCESSION to MAPPED_TRAIT_URI.
### efoId parsed from MAPPED_TRAIT_URI. Note that "Orphanet_2445" is a valid efoId,
### since EFO includes other ontologies.
##########################################################################################
library(readr)
require(data.table, quietly=T)

ifile_default <- paste0(Sys.getenv("HOME"), "/../data/gwascatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv")
efofile_default <- "data/efo.tsv"

ofile_default <- "data/gwascat_trait.tsv"

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (efofile <- args[2])
} else if (length(args)==3) {
  (ifile <- args[1])
  (efofile <- args[2])
  (ofile <- args[3])
} else if (length(args)==0) {
  ifile <- ifile_default
  efofile <- efofile_default
  ofile <- ofile_default
} else {
  message("ERROR: Syntax: gwascat_trait.R GWASFILE EFOFILE [OFILE]")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Input EFO: %s", efofile))
writeLines(sprintf("Output: %s", ofile))

trait <- read_delim(ifile, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, old=c("STUDY ACCESSION"), new=c("STUDY_ACCESSION"))
study <- unique(trait[, .(STUDY_ACCESSION, STUDY)])
trait <- unique(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT)])

writeLines(sprintf("Studies missing MAPPED_TRAIT_URI: %d", trait[is.na(MAPPED_TRAIT_URI), uniqueN(STUDY_ACCESSION)]))

filtered_studies <- merge(data.table(STUDY_ACCESSION = trait[is.na(MAPPED_TRAIT_URI), unique(STUDY_ACCESSION)]),
                          study, by="STUDY_ACCESSION", all.x=T, all.y=F)
filtered_studies[, comment := "Filtered due to missing MAPPED_TRAIT_URI."]
write_delim(filtered_studies, "data/filtered_studies_trait.tsv", "\t")

trait <- unique(trait[!is.na(MAPPED_TRAIT_URI)])

# Split comma separated vals.
trait_multi <- trait[grepl(",", MAPPED_TRAIT_URI)]
trait <- trait[!is.na(MAPPED_TRAIT_URI) & !grepl(",", MAPPED_TRAIT_URI)]
for (i in 1:nrow(trait_multi)) {
  uris <- strsplit(trait_multi$MAPPED_TRAIT_URI[i], ', ', perl=T)[[1]]
  traits <- strsplit(trait_multi$MAPPED_TRAIT[i], ', ', perl=T)[[1]]
  accs <- rep(trait_multi$STUDY_ACCESSION[i], length(uris))
  if (length(uris)!=length(traits)) {
    message(sprintf("ERROR: length(uris)!=length(traits) (%d!=%d) (probably due to commas in trait names) \"%s\"", length(uris), length(traits), trait_multi$MAPPED_TRAIT[i]))
    traits <- rep(trait_multi$MAPPED_TRAIT[i], length(uris)) #Commas in trait names, so must be curated manually.
  }
  trait <- rbind(trait, data.frame(STUDY_ACCESSION=accs, MAPPED_TRAIT_URI=uris, MAPPED_TRAIT=traits))
}

trait[['efoId']] <- as.factor(sub("^.*/", "", trait$MAPPED_TRAIT_URI)) 
trait[['ontology']] <- as.factor(sub("_.*$", "", trait$efoId))
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
# Here we could consider parentage. For a given query trait, perhaps all child traits should be
# included in the associated genes. But then which query traits are allowed? How high level?
# Need statistics mapping all EFO to GWAS Catalog with parentage considered.
#
efo <- read_delim(efofile, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo <- efo[node_or_edge == "node"]
efo[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]
efo[['ontology']] <- as.factor(sub("_.*$", "", efo$id))
efo_counts <- efo[, .N, by="ontology"][order(-N)]
message(sprintf("%12s: %4d / %4d (%4.1f%%)\n", efo_counts$ontology, efo_counts$N, sum(efo_counts$N), 100*efo_counts$N/sum(efo_counts$N)))
efo[, ontology := NULL]

trait <- merge(trait, efo[, .(id, efo_label = label)], by.x="efoId", by.y="id", all.x=T, all.y=F)

trait_unmapped <- trait[is.na(efo_label), .(efoId, MAPPED_TRAIT)]
n_trait_mapped <- nrow(trait[!is.na(efo_label)])
message(sprintf("EFO+ IDs mapped to efo.owl: %d / %d (%.1f%%)", n_trait_mapped, 
                n_trait_mapped + uniqueN(trait_unmapped$efoId),
                100 * (n_trait_mapped / (n_trait_mapped + uniqueN(trait_unmapped$efoId)))))
#
#
write_delim(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT, efoId, efo_label)], ofile, "\t")

#!/usr/bin/env Rscript
##########################################################################################
### Analyze/describe assn file traits.
### Note that EFO includes entities and IDs from other ontologies.
### Output gwascat_trait.tsv for use by tiga_gt_stats.R, but only for mapping from
### STUDY_ACCESSION to MAPPED_TRAIT_URI.
### efoId parsed from MAPPED_TRAIT_URI. Note that "Orphanet_2445" is a valid efoId,
### since EFO includes other ontologies.
##########################################################################################
### Also create efo_sub_gwas.tsv with subclasses for GWAS/EFO traits.
##########################################################################################
library(readr)
library(data.table, quietly=T)

#ifile_default <- paste0(Sys.getenv("HOME"), "/../data/GWASCatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv")
#ifile_default <- paste0(Sys.getenv("HOME"), "/../data/GWASCatalog/data/gwas_catalog_v1.0.2-studies_r2020-07-14.tsv")
ifile_default <- paste0(Sys.getenv("HOME"), "/../data/GWASCatalog/releases/2020/07/15/gwas-catalog-studies_ontology-annotated.tsv")

efofile_default <- "data/efo.tsv"
ofile_default <- "data/gwascat_trait.tsv"
ofile_subclass_default <- "data/efo_sub_gwas.tsv"

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==4) {
  (ifile <- args[1])
  (efofile <- args[2])
  (ofile <- args[3])
  (ofile_subclass <- args[4])
} else if (length(args)==0) {
  ifile <- ifile_default
  efofile <- efofile_default
  ofile <- ofile_default
  ofile_subclass <- ofile_subclass_default
} else {
  message("ERROR: Syntax: gwascat_trait.R GWASFILE EFOFILE OFILE OFILE_SUBCLASS")
  quit()
}
message(sprintf("Input: %s", ifile))
message(sprintf("Input EFO: %s", efofile))
message(sprintf("Output: %s", ofile))

trait <- read_delim(ifile, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, old=c("STUDY ACCESSION"), new=c("STUDY_ACCESSION"))
study <- unique(trait[, .(STUDY_ACCESSION, STUDY)])
trait <- unique(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT)])

writeLines(sprintf("Studies missing MAPPED_TRAIT_URI: %d", trait[is.na(MAPPED_TRAIT_URI), uniqueN(STUDY_ACCESSION)]))

filtered_studies <- merge(data.table(STUDY_ACCESSION = trait[is.na(MAPPED_TRAIT_URI), unique(STUDY_ACCESSION)]),
                          study, by="STUDY_ACCESSION", all.x=T, all.y=F)
filtered_studies[, reason := "missing MAPPED_TRAIT_URI"]
write_delim(filtered_studies, "data/filtered_studies_trait.tsv", "\t")
system("gzip -f data/filtered_studies_trait.tsv")

trait <- unique(trait[!is.na(MAPPED_TRAIT_URI)])

###
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
trait[['EFO_prefix']] <- as.factor(sub("_.*$", "", trait$efoId))
trait <- unique(trait)

message(sprintf("Total trait unique ID count (MAPPED_TRAIT_URI): %d", trait[, uniqueN(MAPPED_TRAIT_URI)]))
message(sprintf("Total trait instance count (MAPPED_TRAIT_URI): %d", trait[!is.na(MAPPED_TRAIT_URI), .N]))

###
#
###
trait_prefix <- trait[, .(N_trait = .N), by="EFO_prefix"][order(-N_trait)]
message(sprintf("%12s: %4d / %4d (%4.1f%%)\n", trait_prefix$EFO_prefix, trait_prefix$N_trait, sum(trait_prefix$N_trait), 100*trait_prefix$N_trait/sum(trait_prefix$N_trait)))
#
trait_study_counts <- trait[, .(N_study = uniqueN(STUDY_ACCESSION)), by="MAPPED_TRAIT_URI"][order(-N_study)]
for (i in 1:10) {
  if (i %in% trait_study_counts$N_study) {
    message(sprintf("%4d traits involve %3d studies", trait_study_counts[N_study==i, uniqueN(MAPPED_TRAIT_URI)], i))
  }
}
message(sprintf("%4d traits involve [%d,%d] studies", trait_study_counts[N_study>10, uniqueN(MAPPED_TRAIT_URI)], 11, max(trait_study_counts$N_study)))


###
# EFO (full ontology)
#
efo <- read_delim(efofile, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_node <- efo[node_or_edge == "node"]
efo_node[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]
efo_node[['EFO_prefix']] <- as.factor(sub("_.*$", "", efo_node$id))
efo_counts <- efo_node[, .(N_class = .N), by="EFO_prefix"][order(-N_class)]
print(efo_counts[1:10])
message(sprintf("Other prefix (total=%d): %d", uniqueN(efo_counts$EFO_prefix), efo_counts[11:nrow(efo_counts), sum(N_class)]))
efo_node[, EFO_prefix := NULL]

trait <- merge(trait, efo_node[, .(id, efo_label = label)], by.x="efoId", by.y="id", all.x=T, all.y=F)

trait_unmapped <- trait[is.na(efo_label), .(efoId, MAPPED_TRAIT)]
n_trait_mapped <- nrow(trait[!is.na(efo_label)])
message(sprintf("Trait IDs mapped to EFO: %d / %d (%.1f%%)", n_trait_mapped, 
                n_trait_mapped + uniqueN(trait_unmapped$efoId),
                100 * (n_trait_mapped / (n_trait_mapped + uniqueN(trait_unmapped$efoId)))))
#
write_delim(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT, efoId, efo_label)], ofile, "\t")
#
###
# Subclass file, of study pairs with subclass-related traits:
efo_sub <- efo[node_or_edge == "edge" & label=="has_subclass"]
efo_sub[, `:=`(node_or_edge = NULL, id = NULL, label = NULL, uri = NULL, comment = NULL)]
setnames(efo_sub, old=c("source", "target"), new=c("trait_uri", "subclass_uri"))
efo_sub <- efo_sub[trait_uri %in% trait$MAPPED_TRAIT_URI & subclass_uri %in% trait$MAPPED_TRAIT_URI]
efo_sub <- merge(efo_sub, efo_node[, .(trait_uri = uri, trait = label)], by="trait_uri")
efo_sub <- merge(efo_sub, efo_node[, .(subclass_uri = uri, subclass_trait = label)], by="subclass_uri")
efo_sub <- merge(efo_sub, unique(trait[, .(study_accession = STUDY_ACCESSION, trait_uri = MAPPED_TRAIT_URI)]), by="trait_uri", allow.cartesian = T)
efo_sub <- merge(efo_sub, unique(trait[, .(study_accession_subclass = STUDY_ACCESSION, subclass_uri = MAPPED_TRAIT_URI)]), by="subclass_uri", allow.cartesian = T)

efo_sub <- unique(efo_sub[, .(study_accession, trait, trait_uri, study_accession_subclass, subclass_trait, subclass_uri)])

write_delim(efo_sub, ofile_subclass, "\t")

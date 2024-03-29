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
  message("ERROR: Syntax: gwascat_trait.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}

ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
#
# FILESET CHANGED 2024-01-22 (see Elliot Sollis email, blog post https://ebispot.github.io/gwas-blog/cohorts-now-available)
#ifile <- paste0(Sys.getenv("HOME"), sprintf("/../data/GWASCatalog/releases/%d/%02d/%02d/gwas-catalog-studies_ontology-annotated.tsv", rel_y, rel_m, rel_d))
ifile <- paste0(Sys.getenv("HOME"), sprintf("/../data/GWASCatalog/releases/%d/%02d/%02d/gwas-catalog-studies-download-alternative-v1.0.2.1.txt", rel_y, rel_m, rel_d))

ifile_efo <- paste0(ODIR, "/efo.tsv")
ofile <- paste0(ODIR, "/gwascat_trait.tsv")
ofile_subclass <- paste0(ODIR, "/efo_sub_gwas.tsv")

message(sprintf("Input: %s", ifile))
message(sprintf("Input EFO: %s", ifile_efo))
message(sprintf("Output: %s", ofile))
message(sprintf("Output subclass: %s", ofile_subclass))

# escape_double=F needed to parse all lines!
trait <- read_delim(ifile, "\t", col_types=cols(.default=col_character()), escape_double=F)
setDT(trait)
setnames(trait, old=c("STUDY ACCESSION", "DISEASE/TRAIT"), new=c("STUDY_ACCESSION", "TRAIT"))
study <- unique(trait[, .(STUDY_ACCESSION, STUDY)])
trait <- unique(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT, TRAIT)])

writeLines(sprintf("Studies missing MAPPED_TRAIT_URI: %d", trait[is.na(MAPPED_TRAIT_URI), uniqueN(STUDY_ACCESSION)]))

###
# Split comma separated vals.
trait_multi <- trait[grepl(",", MAPPED_TRAIT_URI, TRAIT)]
trait <- trait[!grepl(",", MAPPED_TRAIT_URI, TRAIT)]
for (i in 1:nrow(trait_multi)) {
  uris <- strsplit(trait_multi$MAPPED_TRAIT_URI[i], ', ', perl=T)[[1]]
  mapped_traits <- strsplit(trait_multi$MAPPED_TRAIT[i], ', ', perl=T)[[1]]
  accs <- rep(trait_multi$STUDY_ACCESSION[i], length(uris))
  traits <- rep(trait_multi$TRAIT[i], length(uris))
  if (length(uris)!=length(mapped_traits)) {
    message(sprintf("ERROR: length(uris)!=length(mapped_traits) (%d!=%d) (probably due to commas in trait names) \"%s\"", length(uris), length(mapped_traits), trait_multi$MAPPED_TRAIT[i]))
    mapped_traits <- rep(trait_multi$MAPPED_TRAIT[i], length(uris)) #Commas in trait names, so must be curated manually.
  }
  trait <- rbind(trait, data.frame(STUDY_ACCESSION=accs, MAPPED_TRAIT_URI=uris, MAPPED_TRAIT=mapped_traits, TRAIT=traits))
}

trait <- unique(trait)
trait[, efoId := as.factor(sub("^.*/", "", MAPPED_TRAIT_URI))] 
trait[, EFO_prefix := as.factor(sub("_.*$", "", efoId))]

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
efo <- read_delim(ifile_efo, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_node <- efo[node_or_edge == "node"]
efo_node[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]
efo_node[['EFO_prefix']] <- as.factor(sub("_.*$", "", efo_node$id))
efo_counts <- efo_node[, .(N_class = .N), by="EFO_prefix"][order(-N_class)]
print(efo_counts[1:10])
message(sprintf("Other prefix (total=%d): %d", uniqueN(efo_counts$EFO_prefix), efo_counts[11:nrow(efo_counts), sum(N_class)]))
efo_node[, EFO_prefix := NULL]

trait <- merge(trait, efo_node[, .(id, efo_label = label)], by.x="efoId", by.y="id", all.x=T, all.y=F)

#

n_trait_unmapped <- trait[is.na(MAPPED_TRAIT_URI), uniqueN(TRAIT)]
n_trait_mapped <- trait[!is.na(MAPPED_TRAIT_URI), uniqueN(MAPPED_TRAIT_URI)]
message(sprintf("Trait IDs mapped to EFO: %d / %d (%.1f%%)", n_trait_mapped, 
                n_trait_mapped + n_trait_unmapped,
                100 * (n_trait_mapped / (n_trait_mapped + n_trait_unmapped))))
#
write_delim(trait[, .(STUDY_ACCESSION, MAPPED_TRAIT_URI, MAPPED_TRAIT, TRAIT, efoId, efo_label)], ofile, "\t")
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

#Example
unique(efo_sub[grepl("(mood|bipolar)", trait), .(efoId = sub("^.*/", "", trait_uri), subclass_efoId = sub("^.*/", "", subclass_uri), N_study = uniqueN(study_accession), N_study_subclass = uniqueN(study_accession_subclass)), by=c("trait", "subclass_trait")])
#
write_delim(efo_sub, ofile_subclass, "\t")


#!/usr/bin/env Rscript
##########################################################################################
### Analyze/describe assn file traits.
##########################################################################################
library(readr)
require(data.table, quietly=T)

ifile_efo <- "data/efo.tsv"

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  (ifile <- args[1])
} else if (length(args)==0) {
  ifile <- "/home/data/gwascatalog/data/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv"
} else {
  message("ERROR: Syntax: gwascat_trait.R GWASFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))

trait <- read_delim(ifile, "\t")
setDT(trait)

trait <- trait[, .(MAPPED_TRAIT, MAPPED_TRAIT_URI)]

trait <- unique(trait[complete.cases(trait),])

# Split comma separated vals.
trait_multi <- trait[grepl(",", trait$MAPPED_TRAIT_URI)]
trait <- trait[!grepl(",", trait$MAPPED_TRAIT_URI)]
for (i in 1:nrow(trait_multi)) {
  uris <- strsplit(as.character(trait_multi$MAPPED_TRAIT_URI[i]), ', ', perl=T)[[1]]
  traits <- strsplit(as.character(trait_multi$MAPPED_TRAIT[i]), ', ', perl=T)[[1]]
  if (length(uris)!=length(traits)) {
    message(sprintf("ERROR: length(uris)!=length(traits) (%d!=%d) \"%s\"", length(uris), length(traits), trait_multi$MAPPED_TRAIT[i]))
    traits <- NA #Commas in trait names, so must be curated manually.
  }
  trait <- rbind(trait, data.frame(MAPPED_TRAIT=traits, MAPPED_TRAIT_URI=uris))
}

trait[['ontology']] <- as.factor(sub("[^/]+$","", trait$MAPPED_TRAIT_URI))
trait[['id']] <- as.factor(sub("^.*/","", trait$MAPPED_TRAIT_URI))

trait <- unique(trait)

###
#v1.0.2 counts:
#HP: http://purl.obolibrary.org/obo/: 295
#EFO: http://www.ebi.ac.uk/efo/: 3899
#Orphanet: http://www.orpha.net/ORDO/: 78
###

tbl <- table(trait$ontology)
message(sprintf("Ontology: %32s: %4d / %4d (%4.1f%%)\n", names(tbl), tbl, sum(tbl), 100*tbl/sum(tbl)))

###
efo <- read_delim(ifile_efo, "\t")
setDT(efo)
efo <- efo[node_or_edge == "node" & grepl("^EFO", id)]
efo[, `:=`(node_or_edge = NULL, source = NULL, target = NULL)]

trait <- merge(trait, efo, by="id", all.x=T, all.y=F)

trait_unmapped <- trait[grepl("^EFO", id) & is.na(label), .(id, MAPPED_TRAIT)]
n_trait_mapped <- nrow(trait[grepl("^EFO", id) & !is.na(label)])

message(sprintf("EFO IDs mapped to efo.owl: %d / %d (%.1f%%)", n_trait_mapped, 
                n_trait_mapped + uniqueN(trait_unmapped$id),
                100 * (n_trait_mapped / (n_trait_mapped + uniqueN(trait_unmapped$id)))))
#
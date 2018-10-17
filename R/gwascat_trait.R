#!/usr/bin/env Rscript
##########################################################################################
### Analyze/describe assn file traits.
##########################################################################################
require(dplyr, quietly=T)
library(readr)


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

trait <- trait[,c('MAPPED_TRAIT', 'MAPPED_TRAIT_URI')]

trait <- unique(trait[complete.cases(trait),])

# Split comma separated vals.
trait_multi <- trait[grepl(",", trait$MAPPED_TRAIT_URI),]
trait <- trait[!grepl(",", trait$MAPPED_TRAIT_URI),]
for (i in 1:nrow(trait_multi)) {
  uris <- strsplit(as.character(trait_multi$MAPPED_TRAIT_URI[i]), ', ', perl=T)[[1]]
  traits <- strsplit(as.character(trait_multi$MAPPED_TRAIT[i]), ', ', perl=T)[[1]]
  if (length(uris)!=length(traits)) {
    writeLines(sprintf("ERROR: length(uris)!=length(traits) (%d!=%d) \"%s\"", length(uris), length(traits), trait_multi$MAPPED_TRAIT[i]))
    traits <- NA #Commas in trait names, so must be curated manually.
  }
  trait <- rbind(trait, data.frame(MAPPED_TRAIT=traits, MAPPED_TRAIT_URI=uris))
}

trait[['ontology']] <- as.factor(sub("[^/]+$","", trait$MAPPED_TRAIT_URI))
trait[['id']] <- as.factor(sub("^.*/","", trait$MAPPED_TRAIT_URI))

###
#v1.0.2 counts:
#HP: http://purl.obolibrary.org/obo/: 295
#EFO: http://www.ebi.ac.uk/efo/: 3899
#Orphanet: http://www.orpha.net/ORDO/: 78
###

tbl <- table(trait$ontology)
writeLines(sprintf("Ontology: %32s: %4d / %4d (%4.1f%%)", names(tbl), tbl, sum(tbl), 100*tbl/sum(tbl)))

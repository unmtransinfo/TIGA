#!/usr/bin/env Rscript
###
#
library(readr)
library(data.table)

ifile_trait <- "data/gwascat_trait.tsv" #trait2study
ifile_efo <- "data/efo.tsv" #from efo.owl
ifile_efo2do <- "data/oxo_efo2doid_d2.tsv" #from OxO
ifile_efo_gwas <- "data/efo_sub_gwas.tsv"

trait <- read_delim(ifile_trait, "\t", col_types=cols(.default=col_character()))
setDT(trait)
setnames(trait, c("STUDY_ACCESSION", "TRAIT_URI", "TRAIT", "trait_id", "efo_label"))
trait <- trait[!is.na(TRAIT_URI)]
trait[, TRAIT := iconv(TRAIT, from="latin1", to="UTF-8")]

# GWAS traits and EFO
# EFO = Experimental Factor Ontology. Includes GO, Orphanet, PO, Mondo and Uberon classes.
# TSV from source OWL. 

efo <- read_delim(ifile_efo, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_classes <- efo[node_or_edge=="node"][, `:=`(node_or_edge=NULL, source=NULL, target=NULL)]
efo_classes[, in_gwascat := as.logical(uri %in% trait$TRAIT_URI)]
efo_counts <- efo_classes[, Ontology := sub("_.*$", "", id)][, .(N_in_gwas = sum(in_gwascat), N_total = .N), by=Ontology][order(-N_in_gwas)]
message(sprintf("EFO total classes: %d; in GWAS: %d", uniqueN(efo_classes$uri),
sum(efo_counts$N_in_gwas)))

## GWAS trait-subclass relationships 

efo_sub <- efo[node_or_edge=="edge"][, `:=`(node_or_edge=NULL, id=NULL, comment=NULL, uri=NULL)]
message(sprintf("EFO classes: %d ; total subclass relationships: %d",
uniqueN(efo$uri), nrow(efo_sub)))
efo_sub[, source_in_gwas := (source %in% trait$TRAIT_URI)]
efo_sub[, target_in_gwas := (target %in% trait$TRAIT_URI)]
efo_sub <- efo_sub[(source_in_gwas) & (target_in_gwas)]
efo_sub[, `:=`(source_in_gwas=NULL, target_in_gwas=NULL, label=NULL)]
setnames(efo_sub, old=c("source", "target"), new=c("trait_uri", "subclass_uri"))
trait_counts <- trait[, .(N_study = uniqueN(STUDY_ACCESSION)), by=c("TRAIT_URI", "TRAIT")]
efo_sub <- merge(efo_sub, trait_counts[, .(TRAIT_URI, trait_N_study=N_study, trait_name=TRAIT)], by.x="trait_uri", by.y="TRAIT_URI", all.x=T, all.y=F)
efo_sub <- merge(efo_sub, trait_counts[, .(TRAIT_URI, subclass_N_study=N_study, subclass_name=TRAIT)], by.x="subclass_uri", by.y="TRAIT_URI", all.x=T, all.y=F)
efo_sub[, `:=`(trait_id = sub("^.*/", "", trait_uri), subclass_id = sub("^.*/", "", subclass_uri))]
efo_sub <- efo_sub[, .(trait_id, trait_name, subclass_id, subclass_name, trait_N_study, subclass_N_study)]
message(sprintf("GWAS trait-subclass pairs: %d", nrow(efo_sub)))
setorder(efo_sub, -trait_N_study, -subclass_N_study)

## GWAS studies related by EFO subclass links.

efo_gwas <- read_delim(ifile_efo_gwas, "\t", col_types=cols(.default=col_character()))
setDT(efo_gwas)
message(sprintf("GWAS studies related by EFO-subclass: %d", uniqueN(union(efo_gwas[study_accession != study_accession_subclass, study_accession], efo_gwas[study_accession != study_accession_subclass, study_accession_subclass]))))

# EFO to DOID (Disease Ontology ID)
# From EBI Ontology Xref Service (OxO).
# One-to-many and many-to-one mappings exist.
# Keep only closest mappings, maximum distance=2.

efo2do <- read_delim(ifile_efo2do, "\t", col_types=cols(.default=col_character()))
setDT(efo2do)
efo2do[, `:=`(curie_id = sub(":", "_", curie_id, fixed=T), mapped_curie = sub(":", "_", mapped_curie, fixed=T), mapping_source_prefix = NULL, mapping_target_prefix = NULL)]
setnames(efo2do, old=c("curie_id", "label",  "mapped_curie", "mapped_label"), new=c("efo_id", "efo_name", "do_id", "doid_name"))
efo2do[, do_uri := paste0("http://purl.obolibrary.org/obo/", do_id)]
efo2do <- efo2do[, .SD[distance == min(distance)], by=c("efo_id", "do_id")] #keep only closest mappings
gwas_efo_ids <- unique(sub("^.*/", "", efo_gwas$trait_uri))
efo2do_gwas <- efo2do[(efo_id %in% gwas_efo_ids)]
setorder(efo2do_gwas, efo_id, do_id)
message(sprintf("GWAS EFO_IDs (Total): %d", length(gwas_efo_ids)))
message(sprintf("GWAS EFO_ID to DO_ID mappings (distance<=2): EFO_IDs: %d, DO_IDs: %d", 
        uniqueN(efo2do_gwas[, efo_id]), uniqueN(efo2do_gwas[, do_id])))
message(sprintf("GWAS EFO_ID to DO_ID mappings (efo_name=doid_name): EFO_IDs: %d, DO_IDs: %d", 
        uniqueN(efo2do_gwas[efo_name==doid_name, efo_id]),
uniqueN(efo2do_gwas[efo_name==doid_name, do_id])))
message(sprintf("GWAS EFO_ID to DO_ID mappings (distance=1): EFO_IDs: %d, DO_IDs: %d", 
        uniqueN(efo2do_gwas[distance==1, efo_id]), uniqueN(efo2do_gwas[distance==1, do_id])))
message(sprintf("GWAS EFO_ID to DO_ID mappings (distance=2): EFO_IDs: %d, DO_IDs: %d", 
        uniqueN(efo2do_gwas[distance==2, efo_id]), uniqueN(efo2do_gwas[distance==2, do_id])))
#
write_delim(efo2do_gwas[, .(efo_id, efo_name, do_id, doid_name, distance)], "data/gwascat_efo2doid.tsv", "\t")

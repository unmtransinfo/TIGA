#!/usr/bin/env Rscript
########################################################################################
### tiga_final_files.R
### Final files generation for Shiny app and downloads.
########################################################################################
### Input files:
###	gt_stats.tsv.gz	(from tiga_gt_stats.R) 
###	gt_provenance.tsv.gz	(from tiga_gt_provenance.R)
###	filtered_studies.tsv	(from tiga_gt_prepfilter.R)
###	filtered_traits.tsv	(from tiga_gt_prepfilter.R)
###	filtered_genes.tsv	(from tiga_gt_prepfilter.R)
###	gwascat_gwas.tsv	(from gwascat_gwas.R)
###	efo_graph.graphml.gz	(from efo_graph.R)
### Output files:
###	tiga.Rdata		(for Shiny app.R)
###	tiga_gene-trait_stats.tsv	(download)
###	tiga_gene-trait_provenance.tsv	(download)
###	tiga_genes.tsv	(download)
###	tiga_traits.tsv	(download)
########################################################################################
library(readr)
library(data.table)
library(igraph, quietly=T)
#
pkgs <- names(sessionInfo()$otherPkgs)
pkgVerTxt <- paste(sprintf("%s %s", pkgs, sapply(pkgs, function(p){paste(packageVersion(p), collapse=".")})), collapse="; ")
message(pkgVerTxt)
#
efoId2Uri <- function(efoId) { #(non-vector)
  if (is.null(efoId)) { return(NA) }
  if (grepl("^EFO_", efoId)) {
    return(sprintf("http://www.ebi.ac.uk/efo/%s", efoId))
  } else if (grepl("^Orphanet_", efoId)) {
    return(sprintf("http://www.orpha.net/ORDO/%s", efoId))
  } else if (grepl("(^HP_|^GO_|^CHEBI_|^UBERON_|^NCBITaxon_|^MONDO_|^UO_|^CL_|^PATO_|^HANCESTRO_)", efoId)) {
    return(sprintf("http://purl.obolibrary.org/obo/%s", efoId))
  } else {
    return(NA)
  }
}
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)!=1) {
  message("ERROR: Syntax: tiga_final.R DATADIR\n\tSpecify dir for input and output files.")
  quit()
}
DATADIR <- args[1]
OFILE_RDATA <- paste0(DATADIR, "/tiga.Rdata")
OFILE_GTSTATS <- paste0(DATADIR, "/tiga_gene-trait_stats.tsv")
OFILE_PROVENANCE <- paste0(DATADIR, "/tiga_gene-trait_provenance.tsv")
OFILE_GENES <- paste0(DATADIR, "/tiga_genes.tsv")
OFILE_TRAITS <- paste0(DATADIR, "/tiga_traits.tsv")
#
message(sprintf("DATADIR: %s", DATADIR))
#
MIN_ASSN <- 1
#
gt <- read_delim(paste0(DATADIR, "/gt_stats.tsv.gz"), '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(), geneNtrait=col_integer(), geneNstudy=col_integer(), traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), pvalue_mlog_max=col_double(), or_median=col_double(), n_beta=col_double(), study_N_mean=col_double(), rcras=col_double() , meanRank=col_double(), meanRankScore=col_double()))
setDT(gt)
setnames(gt, old=c("geneIdgTdl"), new=c("TDL"))
#
gt_prov <- read_delim(paste0(DATADIR, "/gt_provenance.tsv.gz"), "\t", col_types=cols(.default=col_character()))
setDT(gt_prov)
#
filtered_studies <- read_delim(paste0(DATADIR, "/filtered_studies.tsv"), "\t")
setDT(filtered_studies)
filtered_studies[, type := "study"]
filtered_traits <- read_delim(paste0(DATADIR, "/filtered_traits.tsv"), "\t")
setDT(filtered_traits)
filtered_traits[, type := "trait"]
filtered_genes <- read_delim(paste0(DATADIR, "/filtered_genes.tsv"), "\t")
setDT(filtered_genes)
filtered_genes[, type := "gene"]
filtered_gene_menu <- filtered_genes$ensemblId #named vector
#names(filtered_gene_menu) <- sprintf("%s:%s", filtered_genes$ensemblSymb, filtered_genes$geneName)
names(filtered_gene_menu) <- filtered_genes$ensemblSymb
#
for (ensemblId in filtered_genes$ensemblId) {
  if (ensemblId %in% unique(gt$ensemblId)) {
    message(sprintf("ERROR: filtered gene exists in dataset: \"%s\"", ensemblId))
  }
}
filtered_genes <- filtered_genes[!(ensemblId %in% unique(gt$ensemblId))]
#
filtered <- rbindlist(list(filtered_studies[, .(type, id=STUDY_ACCESSION, reason)], filtered_traits[, .(type, id=TRAIT_URI, reason)], filtered_genes[, .(type, id=ensemblId, reason)]))
#
trait_table <- gt[, .(trait=first(trait), N_study=first(traitNstudy), traitNgene=uniqueN(ensemblId)),  by="efoId"]
message(sprintf("Traits with N_assn<%d: %d", MIN_ASSN, trait_table[traitNgene<MIN_ASSN, uniqueN(efoId)]))
trait_table <- trait_table[traitNgene>=MIN_ASSN]
trait_table[, trait_uri := sapply(efoId, efoId2Uri)]
trait_table <- trait_table[, .(trait_uri, efoId, trait, N_study, traitNgene)][order(trait)]
#
# trait_menu for menu only.
trait_menu <- trait_table
#
gene_table <- gt[, .(geneSymbol, geneName, geneFamily, TDL, N_study = geneNstudy, geneNtrait = uniqueN(efoId), filtered=F),  by=c("ensemblId")]
gene_table <- rbindlist(list(gene_table, filtered_genes[, .(ensemblId, geneSymbol=ensemblSymb, geneName, geneFamily, TDL, N_study=0, geneNtrait=0, filtered=T)]))
gene_table <- unique(gene_table)[order(geneSymbol)]
#
# Duplicated ENSGs (could be isoforms):
dup_ENSGs <- gene_table[ensemblId %in% gene_table[duplicated(ensemblId), ensemblId], .(ensemblId, geneSymbol, geneName)][order(ensemblId)]
message(sprintf("DEBUG: duplicated ENSGs: %d; removing %d rows", uniqueN(dup_ENSGs$ensemblId), sum(duplicated(gene_table$ensemblId))))
gene_table <- gene_table[!duplicated(ensemblId)]
#
# gene_menu for menu only.
gene_menu <- gene_table
gene_menu[, geneSymbol := ifelse(!is.na(geneSymbol), geneSymbol, ensemblId)] #NAs break autocomplete.
#
study_table <- read_delim(paste0(DATADIR, "/gwascat_gwas.tsv"), "\t", col_types = cols(.default = col_character(), DATE=col_date(), DATE_ADDED_TO_CATALOG=col_date()))
setDT(study_table)
# Filter studies without TIGA evidence:
study_table <- study_table[STUDY_ACCESSION %in% gt_prov$STUDY_ACCESSION]
study_table <- study_table[, .(STUDY_ACCESSION, STUDY, MAPPED_TRAIT, MAPPED_TRAIT_URI, PUBMEDID, DATE_PUBLISHED = DATE, DATE_ADDED_TO_CATALOG)][order(DATE_PUBLISHED)]
#
system(paste0("if [ -f \"", DATADIR, "/efo_graph.graphml.gz\" ]; then gunzip -f ", DATADIR, "/efo_graph.graphml.gz ; fi"))
efoGraph <- read_graph(paste0(DATADIR, "/efo_graph.graphml"), format="graphml")
system(paste0("if [ -f \"", DATADIR, "/efo_graph.graphml\" ]; then gzip -f ", DATADIR, "/efo_graph.graphml ; fi"))
#
GWASCATALOG_RELEASE <- trimws(readr::read_file(paste0(DATADIR, "/gwascat_release.txt")))
EFO_RELEASE <- trimws(readr::read_file(paste0(DATADIR, "/efo_release.txt")))
tcrd_info <- read_delim(paste0(DATADIR, "/tcrd_info.tsv"), "\t")
TCRD_RELEASE <- tcrd_info$data_ver[1]
#
###
message(sprintf("Writing: %s", OFILE_RDATA))
save(gt, trait_table, gene_table, trait_menu, gene_menu, filtered_gene_menu, gt_prov, filtered, study_table, efoGraph, GWASCATALOG_RELEASE, EFO_RELEASE, TCRD_RELEASE, file=OFILE_RDATA)
#
###
message(sprintf("Writing: %s", OFILE_GTSTATS))
write_delim(gt, OFILE_GTSTATS, delim="\t")
message(sprintf("Writing: %s", OFILE_PROVENANCE))
write_delim(gt_prov, OFILE_PROVENANCE, delim="\t")
message(sprintf("Writing: %s", OFILE_GENES))
write_delim(gene_table[(!filtered)], OFILE_GENES, delim="\t")
message(sprintf("Writing: %s", OFILE_TRAITS))
write_delim(trait_table, OFILE_TRAITS, delim="\t")
#

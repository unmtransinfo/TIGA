#!/usr/bin/env Rscript
###
# Subclass based trait-gene associations needed for expected behavior.
# https://igraph.org/r/
# We consider parentage. For a given query trait, all child traits should be
# considered for associated genes. But which query traits are allowed? How high level?
# Need statistics mapping all EFO to GWAS Catalog with parentage considered.
# efo.tsv from efo.owl, produced by iu_idsl_jena jena_utils.
###
library(readr, quietly=T)
library(data.table, quietly=T)
library(igraph, quietly=T)
###
message(paste(commandArgs(), collapse=" "))
#
args <- commandArgs(trailingOnly=TRUE)
#
if (length(args)==3) {
  rel_y <- as.integer(args[1])
  rel_m <- as.integer(args[2])
  rel_d <- as.integer(args[3])
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (file.exists("LATEST_RELEASE.txt")) {
  GC_REL <- trimws(read_file("LATEST_RELEASE.txt"))
  rel_y <- as.integer(sub("\\-.*$", "", GC_REL))
  rel_m <- as.integer(sub("\\d+\\-(\\d+)\\-.*$", "\\1", GC_REL))
  rel_d <- as.integer(sub("\\d+\\-\\d+\\-(\\d+).*$", "\\1", GC_REL))
  message(sprintf("LATEST_RELEASE: %s", GC_REL))
  ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))
} else {
  message("ERROR: Syntax: efo_graph.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
#
efofile <- paste0(ODIR, "/efo.tsv")
efosubgwasfile <- paste0(ODIR, "/efo_sub_gwas.tsv") #from gwascat_trait.R
ofile <- paste0(ODIR, "/efo_graph.graphml")
if (file.exists(paste0(ODIR, "/efo_release.txt"))) {
  EFO_REL <- trimws(read_file(paste0(ODIR, "/efo_release.txt")))
  message(sprintf("efo_release: %s", EFO_REL))
} else {
  message(sprintf("Not found: %s", paste0(ODIR, "/efo_release.txt")))
}
#
message(sprintf("efofile: %s", efofile))
message(sprintf("efosubgwasfile: %s", efosubgwasfile))
message(sprintf("ofile: %s", ofile))
#
efo <- read_delim(efofile, "\t", col_types=cols(.default=col_character()))
setDT(efo)
message(sprintf("nodes: %d; edges: %d", nrow(efo[node_or_edge == "node"]), nrow(efo[node_or_edge == "edge"])))
efo_node <- efo[node_or_edge == "node", .(id, uri, label, comment)]
efo_node[['ontology']] <- sub("_.*$", "", efo_node$id)
efo_counts <- efo_node[, .N, by="ontology"][order(-N)]
for (i in 1:min(10, nrow(efo_counts))) {
  message(sprintf("%2d. %12s: %4d / %4d (%4.1f%%)", i, efo_counts$ontology[i], efo_counts$N[i], sum(efo_counts$N), 100*efo_counts$N[i]/sum(efo_counts$N)))
}
message(sprintf("%2d-%2d. %9s: %4d / %4d (%4.1f%%)", i+1, nrow(efo_counts), "OTHER", sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T), sum(efo_counts$N), 
          100*sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T)/sum(efo_counts$N)))
#
efo_edge <- efo[node_or_edge == "edge", .(source, target, label)]
# Any nodes without edges?
orphans <- setdiff(efo_node[, uri], c(efo_edge[, source], efo_edge[, target]))
message(sprintf("nodes with edges: %d; without edges: %d", length(unique(c(efo_edge[, source], efo_edge[, target]))), length(orphans)))
#igraph better with integer vertex IDs?
efo_node[['vId']] <- 1:nrow(efo_node)
efo_edge <- merge(efo_edge, efo_node[, .(uri, sourceId=vId)], by.x="source", by.y="uri", all.x=T, all.y=F)
efo_edge <- merge(efo_edge, efo_node[, .(uri, targetId=vId)], by.x="target", by.y="uri", all.x=T, all.y=F)
efo_edge <- efo_edge[!is.na(sourceId) & !is.na(targetId)]
#
efoG <- igraph::graph_from_edgelist(as.matrix(efo_edge[, .(sourceId, targetId)]), directed=T)
message(sprintf("igraph nodes: %d; edges: %d", length(igraph::V(efoG)), length(igraph::E(efoG))))
#
igraph::graph_attr(efoG, "name") <- "EFO"
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoG, "name"), ifelse(is_directed(efoG), "", "NOT_"), vcount(efoG), ecount(efoG)))
if (length(igraph::V(efoG)) != nrow(efo_node)) {
  message(sprintf("ERROR: length(V(efoG)) != nrow(efo_node) (%d != %d)", length(igraph::V(efoG)), nrow(efo_node)))
}
#
efoG <- igraph::set_vertex_attr(efoG, "name", V(efoG), efo_node$id)
efoG <- igraph::set_vertex_attr(efoG, "description", V(efoG), efo_node$label)
efoG <- igraph::set_vertex_attr(efoG, "uri", V(efoG), efo_node$uri)
efoG <- igraph::set_vertex_attr(efoG, "efoId", V(efoG), efo_node$id)
efoG <- igraph::set_vertex_attr(efoG, "ontology", V(efoG), efo_node$ontology)
#
efo_sub_gwas <- read_delim(efosubgwasfile, "\t")
setDT(efo_sub_gwas)
efo_sub_gwas_counts <- unique(efo_sub_gwas[, .(trait, subclass_trait, N_gwas = uniqueN(study_accession), N_sub_gwas = uniqueN(study_accession_subclass)), by=c("trait_uri", "subclass_uri")])
efo_node <- merge(efo_node, efo_sub_gwas_counts[, .(trait_uri, N_gwas)], by.x="uri", by.y="trait_uri", all.x=T, all.y=F)
efo_node[is.na(N_gwas), N_gwas := 0]
efo_node <- unique(efo_node)
efoG <- igraph::set_vertex_attr(efoG, "N_gwas", V(efoG), efo_node$N_gwas)
###
# Save annotated graph to file.
igraph::write_graph(efoG, ofile, format="graphml")
#
###
# See efo_graph_test.R for test of this output.
###

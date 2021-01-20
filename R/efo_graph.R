#!/usr/bin/env Rscript
###
# Subclass based trait-gene associations needed for expected behavior.
# https://igraph.org/r/
# We consider parentage. For a given query trait, all child traits should be
# considered for associated genes. But which query traits are allowed? How high level?
# Need statistics mapping all EFO to GWAS Catalog with parentage considered.
# efo.tsv from efo.owl, produced by iu_idsl_jena jena_utils.
###
library(readr)
library(data.table, quietly=T)
library(igraph, quietly=T)
###
#
args <- commandArgs(trailingOnly=TRUE)
ifile <- ifelse((length(args)>0), args[1], "data/efo.tsv")
ofile <- ifelse((length(args)>1), args[2], "data/efo_graph.graphml")
#
efo <- read_delim(ifile, "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_node <- efo[node_or_edge == "node", .(id, uri, label, comment)]
efo_node[['ontology']] <- sub("_.*$", "", efo_node$id)
efo_counts <- efo_node[, .N, by="ontology"][order(-N)]
for (i in 1:min(10, nrow(efo_counts))) {
  message(sprintf("%2d. %12s: %4d / %4d (%4.1f%%)", i, efo_counts$ontology[i], efo_counts$N[i], sum(efo_counts$N), 100*efo_counts$N[i]/sum(efo_counts$N)))
}
message(sprintf("%2d-%2d. %9s: %4d / %4d (%4.1f%%)", i+1, nrow(efo_counts), "OTHER", sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T), sum(efo_counts$N), 
          100*sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T)/sum(efo_counts$N)))
#
efo <- efo[node_or_edge == "edge", .(source, target, label)]
#igraph better with integer vertex IDs?
efo_node[['vId']] <- 1:nrow(efo_node)
efo <- merge(efo, efo_node[, .(uri, sourceId=vId)], by.x="source", by.y="uri", all.x=T, all.y=F)
efo <- merge(efo, efo_node[, .(uri, targetId=vId)], by.x="target", by.y="uri", all.x=T, all.y=F)
efoG <- igraph::graph_from_edgelist(as.matrix(efo[(!is.na(sourceId) & !is.na(targetId)), .(sourceId, targetId)]), directed=T)
graph_attr(efoG, "name") <- "EFO"
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoG, "name"), ifelse(is_directed(efoG), "", "NOT_"), vcount(efoG), ecount(efoG)))
efoG <- set_vertex_attr(efoG, "name", V(efoG), efo_node$id)
###
#
efoG <- set_vertex_attr(efoG, "description", V(efoG), efo_node$label)
efoG <- set_vertex_attr(efoG, "uri", V(efoG), efo_node$uri)
efoG <- set_vertex_attr(efoG, "efoId", V(efoG), efo_node$id)
efoG <- set_vertex_attr(efoG, "ontology", V(efoG), efo_node$ontology)
ontocolor <- function(ont) {
  ifelse(ont=="EFO", "#CCCCFF",
  ifelse(ont=="Orphanet", "green",
  ifelse(ont=="MONDO", "orange",
  ifelse(ont=="CHEBI", "pink",
  ifelse(ont=="NCBITaxon", "yellow",
  ifelse(ont=="HP", "cyan",
  ifelse(ont=="UBERON", "aqua",
         "#CCCCCC")))))))
}
efoG <- set_vertex_attr(efoG, "color", V(efoG), ontocolor(efo_node$ontology))

###
# Save annotated graph to file for TIGA UI.
write_graph(efoG, ofile, format="graphml")
#
###
# See efo_graph_test.R for test of this output.
###

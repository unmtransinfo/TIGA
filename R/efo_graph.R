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
# Induce and plot subgraph. BFS finds all subclasses.
EFONAME <- "mood disorder" #EFO_0004247 (11 subclasses)
#EFONAME <- "epilepsy" #500+ subclasses!
v_this <- V(efoG)[V(efoG)$description == EFONAME]
bfs_this <- igraph::bfs(efoG, v_this, neimode="out", unreachable=F)
subg <- induced_subgraph(efoG, bfs_this$order[1:sum(!is.na(bfs_this$order))])
graph_attr(subg, "name") <- sprintf("EFO_SUBGRAPH:%s (subclasses)", v_this$efoId)
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(subg, "name"), ifelse(is_directed(subg), "", "NOT_"), vcount(subg), ecount(subg)))
writeLines(sprintf("%2d. %-44s %-18s \"%s\"", 1:vcount(subg), V(subg)$uri, V(subg)$efoId, V(subg)$description))
###
# Plot (if interactive)
if (interactive()) {
	tkplot(subg, canvas.width=800, canvas.height=600, layout=layout_as_tree, vertex.label=sprintf("%s\n%s", V(subg)$efoId, V(subg)$description), 
       	vertex.frame.color="#6666AA", vertex.size=30, edge.color="#6666AA", edge.width=2, edge.arrow.size=1,
       	edge.label="has_subclass")
}
#
# CYJS not supported. (Can do igraph_utils.py graph2cyjs)
subg <- set_vertex_attr(subg, "name", V(subg), sprintf("%s: %s", V(subg)$efoId, V(subg)$description))
write_graph(subg, sprintf("data/efo_subgraph_%s.graphml", efo_node[label == EFONAME]$id), format="graphml")
#
###
# Test that we can read graphml file ok.
efoG2 <- read_graph(ofile, format="graphml")
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoG2, "name"), ifelse(is_directed(efoG2), "", "NOT_"), vcount(efoG2), ecount(efoG2)))
message(sprintf("TEST: %d =? %d (%s)", vcount(efoG2), vcount(efoG), ifelse(vcount(efoG2)==vcount(efoG), "OK", "NOT OK")))
message(sprintf("TEST: %d =? %d (%s)", ecount(efoG2), ecount(efoG), ifelse(ecount(efoG2)==ecount(efoG), "OK", "NOT OK")))
#

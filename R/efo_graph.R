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
efo <- read_delim("data/efo.tsv", "\t", col_types=cols(.default=col_character()))
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
efoGraph <- igraph::graph_from_edgelist(as.matrix(efo[(!is.na(sourceId) & !is.na(targetId)), .(sourceId, targetId)]), directed=T)
graph_attr(efoGraph, "name") <- "EFO"
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoGraph, "name"), ifelse(is_directed(efoGraph), "", "NOT_"), vcount(efoGraph), ecount(efoGraph)))
efoGraph <- set_vertex_attr(efoGraph, "name", V(efoGraph), efo_node$id)
###
#
efoGraph <- set_vertex_attr(efoGraph, "description", V(efoGraph), efo_node$label)
efoGraph <- set_vertex_attr(efoGraph, "uri", V(efoGraph), efo_node$uri)
efoGraph <- set_vertex_attr(efoGraph, "efoId", V(efoGraph), efo_node$id)
efoGraph <- set_vertex_attr(efoGraph, "ontology", V(efoGraph), efo_node$ontology)
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
efoGraph <- set_vertex_attr(efoGraph, "color", V(efoGraph), ontocolor(efo_node$ontology))

###
# Save annotated graph to file for TIGA UI.
write_graph(efoGraph, "data/efo_graph.graphml", format="graphml")
system("gzip -f data/efo_graph.graphml")
#
###
# Induce and plot subgraph. BFS finds all subclasses.
bfs_this <- igraph::bfs(efoGraph, V(efoGraph)[efo_node[label == "mood disorder"]$vId], neimode="out", unreachable=F)
subg <- induced_subgraph(efoGraph, bfs_this$order[1:sum(!is.na(bfs_this$order))])
graph_attr(subg, "name") <- sprintf("EFO_SUBGRAPH:%s (subclasses)", efo_node[label == "mood disorder"]$id)
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(subg, "name"), ifelse(is_directed(subg), "", "NOT_"), vcount(subg), ecount(subg)))
writeLines(sprintf("%2d. %-44s %-18s \"%s\"", 1:vcount(subg), V(subg)$uri, V(subg)$efoId, V(subg)$description))
###
# Plot
tkplot(subg, canvas.width=800, canvas.height=600, layout=layout_as_tree, vertex.label=sprintf("%s\n%s", V(subg)$efoId, V(subg)$description), 
       vertex.frame.color="#6666AA", vertex.size=30, edge.color="#6666AA", edge.width=2, edge.arrow.size=1,
       edge.label="has_subclass")
#
# CYJS not supported. (Can do igraph_utils.py graph2cyjs)
subg <- set_vertex_attr(subg, "name", V(subg), sprintf("%s: %s", V(subg)$efoId, V(subg)$description))
write_graph(subg, sprintf("data/efo_subgraph_%s.graphml", efo_node[label == "mood disorder"]$id), format="graphml")
#
###
# Test that we can read graphml file ok.
system("gunzip -f data/efo_graph.graphml.gz")
efoGraph2 <- read_graph("data/efo_graph.graphml", format="graphml")
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoGraph2, "name"), ifelse(is_directed(efoGraph2), "", "NOT_"), vcount(efoGraph2), ecount(efoGraph2)))
message(sprintf("TEST: %d =? %d (%s)", vcount(efoGraph2), vcount(efoGraph), ifelse(vcount(efoGraph2)==vcount(efoGraph), "OK", "NOT OK")))
message(sprintf("TEST: %d =? %d (%s)", ecount(efoGraph2), ecount(efoGraph), ifelse(ecount(efoGraph2)==ecount(efoGraph), "OK", "NOT OK")))
system("gzip -f data/efo_graph.graphml")
#
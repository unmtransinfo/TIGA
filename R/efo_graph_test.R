#!/usr/bin/env Rscript
###
# GraphML file output from efo_graph.R 
# https://igraph.org/r/
###
#library(readr)
#library(data.table, quietly=T)
library(igraph, quietly=T)
###
#
args <- commandArgs(trailingOnly=TRUE)
ifile <- ifelse((length(args)>0), args[1], "data/efo_graph.graphml")
#

###
efoG <- read_graph(ifile, format="graphml")
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
#write_graph(subg, sprintf("data/efo_subgraph_%s.graphml", efo_node[label == EFONAME]$id), format="graphml")
#
###
# Test that we can read graphml file ok.
efoG2 <- read_graph(ofile, format="graphml")
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(efoG2, "name"), ifelse(is_directed(efoG2), "", "NOT_"), vcount(efoG2), ecount(efoG2)))
message(sprintf("TEST: %d =? %d (%s)", vcount(efoG2), vcount(efoG), ifelse(vcount(efoG2)==vcount(efoG), "OK", "NOT OK")))
message(sprintf("TEST: %d =? %d (%s)", ecount(efoG2), ecount(efoG), ifelse(ecount(efoG2)==ecount(efoG), "OK", "NOT OK")))
#

#!/usr/bin/env Rscript
###
###
library(readr)
library(data.table, quietly=T)
library(igraph, quietly=T)
###
#
dto <- read_delim("data/dto.tsv", "\t", col_types=cols(.default=col_character()))
setDT(dto)
dto_node <- dto[node_or_edge == "node", .(id, uri, label, comment)]
dto_node[['ontology']] <- sub("_.*$", "", dto_node$id)
dto_counts <- dto_node[, .N, by="ontology"][order(-N)]
for (i in 1:nrow(dto_counts)) {
  message(sprintf("%2d. %12s: %4d / %4d (%4.1f%%)", i, dto_counts$ontology[i], dto_counts$N[i], sum(dto_counts$N), 100*dto_counts$N[i]/sum(dto_counts$N)))
}
#
dto <- dto[node_or_edge == "edge", .(source, target, label)]
#igraph better with integer vertex IDs?
dto_node[['vId']] <- 1:nrow(dto_node)
dto <- merge(dto, dto_node[, .(uri, sourceId=vId)], by.x="source", by.y="uri", all.x=T, all.y=F)
dto <- merge(dto, dto_node[, .(uri, targetId=vId)], by.x="target", by.y="uri", all.x=T, all.y=F)
dtoG <- igraph::graph_from_edgelist(as.matrix(dto[(!is.na(sourceId) & !is.na(targetId)), .(sourceId, targetId)]), directed=T)
graph_attr(dtoG, "name") <- "DTO"
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(dtoG, "name"), ifelse(is_directed(dtoG), "", "NOT_"), vcount(dtoG), ecount(dtoG)))
dtoG <- set_vertex_attr(dtoG, "name", V(dtoG), dto_node$id)
###
#
dtoG <- set_vertex_attr(dtoG, "description", V(dtoG), dto_node$label)
dtoG <- set_vertex_attr(dtoG, "uri", V(dtoG), dto_node$uri)
dtoG <- set_vertex_attr(dtoG, "dtoId", V(dtoG), dto_node$id)
nodecolor <- function(dtoId) {
  "orange"
}
dtoG <- set_vertex_attr(dtoG, "color", V(dtoG), nodecolor(dto_node$id))
###
# Save annotated graph.
write_graph(dtoG, "data/dto_graph.graphml", format="graphml")
#
###
# From "Gene" root node classify gene nodes by level.
# Induce and plot subgraph. BFS finds all subclasses.
dtoId_geneRoot <- "DTO_00200000"
v_this <- V(dtoG)[V(dtoG)$name == dtoId_geneRoot]
bfs_this <- igraph::bfs(dtoG, v_this, neimode="out", dist=T, unreachable=F)
subg <- induced_subgraph(dtoG, bfs_this$order[1:sum(!is.na(bfs_this$order) & !is.na(bfs_this$dist) & bfs_this$dist<3)])
graph_attr(subg, "name") <- sprintf("DTO_SUBGRAPH:%s (subclasses)", v_this$dtoId)
message(sprintf("Graph \"%s\" (%sDIRECTED): vertices: %d; edges: %d", graph_attr(subg, "name"), ifelse(is_directed(subg), "", "NOT_"), vcount(subg), ecount(subg)))
writeLines(sprintf("%2d. %-44s %-18s \"%s\"", 1:vcount(subg), V(subg)$uri, V(subg)$dtoId, V(subg)$description))
###
# Plot
tkplot(subg, canvas.width=800, canvas.height=600, layout=layout_as_tree, vertex.label=sprintf("%s\n%s", V(subg)$dtoId, V(subg)$description), 
       vertex.frame.color="#6666AA", vertex.size=30, edge.color="#6666AA", edge.width=2, edge.arrow.size=1,
       edge.label="has_subclass")

#!/usr/bin/env Rscript
###
###
library(readr)
library(data.table, quietly=T)
library(igraph, quietly=T)
###
#
efo <- read_delim("data/dto.tsv", "\t", col_types=cols(.default=col_character()))
setDT(dto)
dto_node <- dto[node_or_edge == "node", .(id, uri, label, comment)]
dto_node[['ontology']] <- sub("_.*$", "", dto_node$id)
dto_counts <- dto_node[, .N, by="ontology"][order(-N)]
for (i in 1:min(10, nrow(dto_counts))) {
  message(sprintf("%2d. %12s: %4d / %4d (%4.1f%%)", i, dto_counts$ontology[i], dto_counts$N[i], sum(dto_counts$N), 100*dto_counts$N[i]/sum(dto_counts$N)))
}
message(sprintf("%2d-%2d. %9s: %4d / %4d (%4.1f%%)", i+1, nrow(dto_counts), "OTHER", sum(dto_counts$N[i+1:nrow(dto_counts)], na.rm=T), sum(dto_counts$N), 
          100*sum(dto_counts$N[i+1:nrow(dto_counts)], na.rm=T)/sum(dto_counts$N)))
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

###
# Save annotated graph.
write_graph(dtoG, "data/dto_graph.graphml", format="graphml")
#
###

#!/usr/bin/env Rscript
###
# Higher level groupings could provide UI visual channel for trait classification.
# General groupings can include child trait-gene associations as needed for expected behavior.
###
library(readr)
library(data.table, quietly=T)
library(igraph, quietly=T)

# File from nx_analysis.py.
efo_groups <- read_delim(paste0(Sys.getenv("HOME"), "/../data/gwascatalog/tiga_data/efo_groups.tsv"), "\t", col_types=cols(.default=col_character(), level=col_integer(), N_sub=col_integer(), N_sub_gwc=col_integer(), in_gwc=col_logical()))
setDT(efo_groups)

setorder(efo_groups, level, -N_sub)

#
###
# Here we could consider parentage. For a given query trait, perhaps all child traits should be
# included in the associated genes. But then which query traits are allowed? How high level?
# Need statistics mapping all EFO to GWAS Catalog with parentage considered.
#
efo <- read_delim("data/efo.tsv", "\t", col_types=cols(.default=col_character()))
setDT(efo)
efo_node <- efo[node_or_edge == "node", .(id, uri, label, comment)]
efo_node[['ontology']] <- as.factor(sub("_.*$", "", efo_node$id))
efo_counts <- efo_node[, .N, by="ontology"][order(-N)]
for (i in 1:min(10, nrow(efo_counts))) {
  message(sprintf("%2d. %12s: %4d / %4d (%4.1f%%)", i, efo_counts$ontology[i], efo_counts$N[i], sum(efo_counts$N), 100*efo_counts$N[i]/sum(efo_counts$N)))
}
message(sprintf("%2d-%2d. %9s: %4d / %4d (%4.1f%%)", i+1, nrow(efo_counts), "OTHER", sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T), sum(efo_counts$N), 
          100*sum(efo_counts$N[i+1:nrow(efo_counts)], na.rm=T)/sum(efo_counts$N)))
#
efo <- efo[node_or_edge == "edge", .(source, target, label)]
efo <- efo[, `:=`(sourceId = sub("^.*/", "", source), targetId = sub("^.*/", "", target))]
efoGraph <- igraph::graph_from_edgelist(as.matrix(efo[, .(sourceId, targetId)]), directed=T)
message(sprintf("iGraph (%sDIRECTED): vertices: %d; edges: %d", ifelse(is_directed(efoGraph), "", "NOT_"), vcount(efoGraph), ecount(efoGraph)))
efo_node <- merge(efo_node, data.table(igraphId=1:length(V(efoGraph)), id=V(efoGraph)$name), by="id")
efoGraph <- set_vertex_attr(efoGraph, "description", V(efoGraph)[efo_node[!is.null(igraphId)]$igraphId], efo_node[!is.null(igraphId)]$label)
efoGraph <- set_vertex_attr(efoGraph, "uri", V(efoGraph)[efo_node[!is.null(igraphId)]$igraphId], efo_node[!is.null(igraphId)]$uri)
efoGraph <- set_vertex_attr(efoGraph, "efoId", V(efoGraph)[efo_node[!is.null(igraphId)]$igraphId], efo_node[!is.null(igraphId)]$id)

# Plot subgraph? Need recursive, all subclasses.
sg <- make_ego_graph(efoGraph, order = 1, nodes = V(efoGraph)[1], mode = "out")[[1]]
message(sprintf("iGraph (%sDIRECTED): vertices: %d; edges: %d", ifelse(is_directed(sg), "", "NOT_"), vcount(sg), ecount(sg)))
plot(sg, layout = layout.auto)

tkplot(sg, layout = layout.auto, vertex.label=sprintf("%s\n%s", V(sg)$efoId, V(sg)$description), 
       vertex.color="#CCCCFF", vertex.frame.color="#FFFFFF", vertex.size=30, 
       edge.label="has_subclass")

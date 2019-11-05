#!/usr/bin/env python3
###
import sys,os
import pandas
import networkx
from networkx.convert_matrix import from_pandas_edgelist

efo_edges = pandas.read_csv("data/efo_edgelist.tsv", "\t")
graph = from_pandas_edgelist(efo_edges, source="source", target="target", edge_attr="edge_attr", create_using=networkx.DiGraph)
print("nodes: {0}; edges: {1}".format(len(graph.nodes), len(graph.edges)), file=sys.stderr)

efo_nodes = pandas.read_csv("data/efo_nodelist.tsv", "\t", index_col="id")
#node_attr = efo_nodes.to_dict(orient='dict')
node_attr = efo_nodes.to_dict(orient='index')

for attr in efo_nodes.columns:

  networkx.set_node_attributes(graph,
	[node_attr[nodeId] for nodeId in graph.nodes],
	key)


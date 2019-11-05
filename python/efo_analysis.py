#!/usr/bin/env python3
"""
	Load efo.owl as NetworkX graph.
	For each class, associate (cluster) with an ancestral group of defined
	size or other criteria.
"""
###
import sys,os,json,logging
import pandas as pd
import networkx as nx
from networkx.convert_matrix import from_pandas_edgelist

###
def SubclassCount(G, n):
   """Recursive all-subclass count."""
   N_sub = 0
   for nn in G.successors(n):
     N_sub += 1
     N_sub += SubclassCount(G, nn)
   return N_sub

def SubclassCount_InSet(G, n, nodeset):
   """Recursive all-subclass-in-set count."""
   N_sub = 0
   for nn in G.successors(n):
     #in_gwc = nx.get_node_attributes(G, 'in_gwc')[nn] if nn in nx.get_node_attributes(G, 'in_gwc') else False
     in_gwc = bool(n in nodeset) #faster?
     if in_gwc:
       N_sub += 1
     N_sub += SubclassCount_InSet(G, nn, nodeset)
   return N_sub

###
def Graph2CYJS(G, ofile):
  cyjs = nx.cytoscape_data(G)
  fout = open(ofile, "w")
  fout.write(json.dumps(cyjs, indent=2))
  fout.close()

###
def Groups2TSV(groups, ofile):
  fout = open(ofile, "w")
  groups.to_csv(fout, '\t', index=False)
  fout.close()

#############################################################################
if __name__=="__main__":
  logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
  ifile = "data/efo_edgelist.tsv"
  logging.info("Reading {0}".format(ifile))
  efo_edges = pd.read_csv(ifile, "\t", dtype=str)
  G = from_pandas_edgelist(efo_edges, source="source", target="target", edge_attr="edge_attr", create_using=nx.DiGraph)
  G.graph['name'] = "EFO: Experimental Factor Ontology"
  logging.info(nx.info(G))
  logging.info("nodes: {0}; edges: {1}; directed: {2}".format(G.number_of_nodes(), G.number_of_edges(), G.is_directed()))
  logging.info("connected: {0}".format(nx.is_weakly_connected(G)))
  logging.info("connected components: {0}".format(nx.number_weakly_connected_components(G)))
  logging.info("DAG (mono-hierarchy): {0}".format(nx.is_directed_acyclic_graph(G)))
  #
  ###
  ifile = "data/efo_nodelist.tsv"
  logging.info("Reading {0}".format(ifile))
  efo_nodes = pd.read_csv(ifile, "\t", index_col="id")
  node_attr = efo_nodes.to_dict(orient='index')
  nx.set_node_attributes(G, node_attr)
  #
  ###
  efoIds=set()
  ofile = "data/gwascatalog.efoid"
  logging.info("Writing {0}".format(ofile))
  with open(ofile) as fin:
    for line in fin:
      efoIds.add(line.strip())
  logging.info("efoIds in GWAS Catalog: {0}".format(len(efoIds)))
  #
  nx.set_node_attributes(G, {efoId:{'in_gwc':True} for efoId in efoIds})
  ###
  ofile = "data/efo_nxgraph.cyjs"
  logging.info("Writing {0}".format(ofile))
  Graph2CYJS(G, ofile)
  #
  ###
  # Find classes with N_subclasses >= MIN_SIZE.
  grouplist=[];
  i_node=0;
  for n in G.nodes:
    i_node += 1
    N_sub = SubclassCount(G, n)
    label = nx.get_node_attributes(G, 'label')[n]
    #in_gwc = nx.get_node_attributes(G, 'in_gwc')[n] if n in nx.get_node_attributes(G, 'in_gwc') else False
    in_gwc = bool(n in efoIds) #faster?
    if N_sub >= 100:
      #N_sub_gwc = SubclassCount_InSet(G, n)
      N_sub_gwc = SubclassCount_InSet(G, n, efoIds)
      logging.debug("{0}/{1}. {2}: {3}; N_sub={4}; N_sub_gwc={5}".format(i_node, G.number_of_nodes(), n, label, N_sub, N_sub_gwc))
      grouplist.append((n, label, in_gwc, N_sub, N_sub_gwc))
  groups = pd.DataFrame({
	'efoId': [efoId for efoId,label,in_gwc,N_sub,N_sub_gwc in grouplist],
	'label': [label for efoId,label,in_gwc,N_sub,N_sub_gwc in grouplist],
	'in_gwc': [in_gwc for efoId,label,in_gwc,N_sub,N_sub_gwc in grouplist],
	'N_sub': [N_sub for efoId,label,in_gwc,N_sub,N_sub_gwc in grouplist],
	'N_sub_gwc': [N_sub_gwc for efoId,label,in_gwc,N_sub,N_sub_gwc in grouplist]})
  	.sort_values(by=['N_sub_gwc', 'N_sub'], ascending=False)

  print(groups.head(10))
  ofile = "data/efo_groups.tsv"
  logging.info("Writing {0}".format(ofile))
  Groups2TSV(groups, ofile)

  print(groups[(groups.in_gwc.bool())].head(10))

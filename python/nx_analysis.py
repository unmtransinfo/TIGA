#!/usr/bin/env python3
"""
Load OWL as NetworkX (NX) graph, for NX analytics.
For each class, associate (cluster) with an ancestral group of defined
size or other criteria.
"""
###
import sys,os,json,logging,argparse,time,tqdm
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
###
def SubclassCount_InSet(G, n, nodeset):
   """Recursive all-subclass-in-set count. Slow."""
   N_sub = 0
   for nn in G.successors(n):
     in_set = bool(n in nodeset)
     if in_set:
       N_sub += 1
     N_sub += SubclassCount_InSet(G, nn, nodeset)
   return N_sub
###
def Level(G, n):
  """Level relative to ancestral root[s] at level 0"""
  roots = [n for n,d in G.in_degree() if d==0] 
  min_pathlen = min([nx.shortest_path_length(G, root, n) if nx.has_path(G, root, n) else sys.maxsize for root in roots])
  return min_pathlen
###
def Graph2CYJS(G, ofile):
  cyjs = nx.cytoscape_data(G)
  fout.write(json.dumps(cyjs, indent=2))
###
def Groups2TSV(groups, fout):
  if groups is None: return
  logging.info(f"Output rows: {groups.shape[0]}; cols: {groups.shape[1]}")
  groups.to_csv(fout, '\t', index=False)
###
def GraphSummary(G):
  logging.info(nx.info(G))
  logging.info(f"nodes: {G.number_of_nodes()}; edges: {G.number_of_edges()}; directed: {G.is_directed()}")
  logging.info(f"connected: {nx.is_weakly_connected(G)}")
  logging.info(f"connected components: {nx.number_weakly_connected_components(G)}")
  logging.info(f"DAG (mono-hierarchy): {nx.is_directed_acyclic_graph(G)}")
  logging.info(f"Tree: {nx.is_tree(G)}")
  logging.info(f"Forest: {nx.is_forest(G)}")
  #
  roots = [n for n,d in G.in_degree() if d==0] 
  logging.info(f"Roots: {len(roots)}")
  for root in roots:
    logging.info("root: {}: {}".format(root, nx.get_node_attributes(G, 'label')[root] if root in nx.get_node_attributes(G, 'label') else ''))
  #
  leafs = [n for n,d in G.out_degree() if d==0] 
  logging.info(f"Leafs: {len(leafs)}")
  #
  singles = set(roots) & set(leafs)
  logging.info(f"Singles: {len(singles)}")
###
def GraphCount_InSet(G, nodeset):
  N_in = 0
  for n in G.nodes:
    if n in nodeset:
      N_in+=1
  return N_in
###
def Cluster(G, nodeSet, min_groupsize, max_level, setname):
  """Find ancestors with N_subclasses >= min_groupsize."""
  i_node=0; grouplist=[]; tq=None;
  roots = [n for n,d in G.in_degree() if d==0] 
  logging.debug(f"Roots: {len(roots)}")
  N_inset = GraphCount_InSet(G, nodeSet)
  logging.info(f"nodeSet count: {len(nodeSet)}")
  logging.info(f"Graph nodes in nodeSet: {N_inset}")
  if N_inset==0:
    logging.error("Zero graph nodes in nodeSet.")
    return None
  for n in G.nodes:
    if tq is None: tq = tqdm.tqdm(total=len(G.nodes), unit="nodes")
    tq.update()
    i_node+=1
    level = Level(G, n)
    if level>max_level: continue
    N_sub = SubclassCount(G, n)
    label = nx.get_node_attributes(G, 'label')[n]
    in_set = bool(n in nodeSet)
    #is_root = bool(n in roots)
    N_sub_set = SubclassCount_InSet(G, n, nodeSet)
    logging.debug(f"{i_node}/{G.number_of_nodes()}. {n}: {label}; level={level}; N_sub={N_sub}; N_sub_{setname}={N_sub_set}")
    if N_sub_set >= min_groupsize:
      logging.debug(f"{i_node}/{G.number_of_nodes()}. {n}: {label}; level={level}; N_sub={N_sub}; N_sub_{setname}={N_sub_set} (> min_groupsize={min_groupsize}")
      grouplist.append((n,label,level,in_set,N_sub,N_sub_set))
  groups = pd.DataFrame({
	'Id': [Id for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'label': [label for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'level': [level for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	f'in_{setname}': [in_set for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'N_sub': [N_sub for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	f'N_sub_{setname}': [N_sub_set for Id,label,level,in_set,N_sub,N_sub_set in grouplist]
	}).sort_values(by=[f'N_sub_{setname}', 'N_sub'], ascending=False)
  #print(groups.head(10)) #DEBUG
  print(groups[groups[f'in_{setname}']].head(18)) #DEBUG
  return groups

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description="NetworkX analytics, from edgelist, optional node list, attributes.", epilog="Clustering into groups by common ancestor. Designed for EFO.")
  OPS = ['summary', 'graph2cyjs', 'cluster']
  parser.add_argument("op", choices=OPS, help='OPERATION')
  parser.add_argument("--i_edge", required=True, dest="ifile_edge", help="Input edgelist (TSV)")
  parser.add_argument("--i_node_set", dest="ifile_node_set", help="Input node set (IDs)")
  parser.add_argument("--i_node_attr", dest="ifile_node_attr", help="Input node attributes (TSV)")
  parser.add_argument("--o", dest="ofile", help="output (TSV|CYJS|etc.)")
  parser.add_argument("--setname", default="myset", help="Node set name")
  parser.add_argument("--graphname", help="Assign this name to graph.")
  parser.add_argument("--min_groupsize", type=int, default=100)
  parser.add_argument("--max_level", type=int, default=10)
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG if args.verbose>0 else logging.INFO)

  t0 = time.time()

  logging.info(f"Reading {args.ifile_edge}")
  efo_edges = pd.read_csv(args.ifile_edge, "\t", dtype=str)
  G = from_pandas_edgelist(efo_edges, source="source", target="target", edge_attr="edge_attr", create_using=nx.DiGraph)
  G.graph['name'] = args.graphname

  fout = open(args.ofile, "w") if args.ofile else sys.stdout

  ###
  if args.ifile_node_attr:
    logging.info(f"Reading {args.ifile_node_attr}")
    efo_nodes = pd.read_csv(args.ifile_node_attr, "\t", index_col="id")
    node_attr = efo_nodes.to_dict(orient='index')
    nx.set_node_attributes(G, node_attr)
  #
  ###
  if args.ifile_node_set:
    nodeSetIds = set()
    logging.info(f"Reading {args.ifile_node_set}")
    with open(args.ifile_node_set) as fin:
      for line in fin:
        nodeSetIds.add(line.strip())
    logging.info(f"nodeSetIds in {args.setname}: {len(nodeSetIds)}")
    nx.set_node_attributes(G, {nodeId:{'in_%s'%args.setname:True} for nodeId in nodeSetIds})

  ###
  if args.op == 'summary':
    GraphSummary(G)
  #
  elif args.op == 'graph2cyjs':
    logging.info(f"Writing {args.ofile}")
    Graph2CYJS(G, fout)
  #
  ###
  elif args.op == 'cluster':
    if not args.ifile_node_set:
      parser.error(f"--i_node_set required for: {args.op}")
    GraphSummary(G)
    groups = Cluster(G, nodeSetIds, args.min_groupsize, args.max_level, args.setname)
    Groups2TSV(groups, fout)

  logging.info(('Elapsed time: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)))))


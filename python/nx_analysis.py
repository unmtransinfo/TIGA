#!/usr/bin/env python3
"""
	Load efo.owl as NetworkX (NX) graph, for NX analytics.
	For each class, associate (cluster) with an ancestral group of defined
	size or other criteria.
"""
###
import sys,os,json,logging,argparse,time
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
  groups.to_csv(fout, '\t', index=False)
###
def GraphSummary(G):
  logging.info(nx.info(G))
  logging.info("nodes: {0}; edges: {1}; directed: {2}".format(G.number_of_nodes(), G.number_of_edges(), G.is_directed()))
  logging.info("connected: {0}".format(nx.is_weakly_connected(G)))
  logging.info("connected components: {0}".format(nx.number_weakly_connected_components(G)))
  logging.info("DAG (mono-hierarchy): {0}".format(nx.is_directed_acyclic_graph(G)))
  logging.info("Tree: {0}".format(nx.is_tree(G)))
  logging.info("Forest: {0}".format(nx.is_forest(G)))
  #
  roots = [n for n,d in G.in_degree() if d==0] 
  logging.info("Roots: {0}".format(len(roots)))
  for root in roots:
    logging.info("root: {0}: {1}".format(root, nx.get_node_attributes(G, 'label')[root] if root in nx.get_node_attributes(G, 'label') else ''))
  #
  leafs = [n for n,d in G.out_degree() if d==0] 
  logging.info("Leafs: {0}".format(len(leafs)))
  #
  singles = set(roots) & set(leafs)
  logging.info("Singles: {0}".format(len(singles)))
###
def Cluster(G, nodeSet, min_groupsize, setname):
  roots = [n for n,d in G.in_degree() if d==0] 
  grouplist=[];
  i_node=0;
  for n in G.nodes:
    i_node += 1
    N_sub = SubclassCount(G, n)
    label = nx.get_node_attributes(G, 'label')[n]
    in_set = bool(n in nodeSet)
    #is_root = bool(n in roots)
    level = Level(G, n)
    N_sub_set = SubclassCount_InSet(G, n, nodeSet)
    if N_sub_set >= min_groupsize:
      logging.debug("{0}/{1}. {2}: {3}; level={4}; N_sub={5}; N_sub_{6}={7}".format(i_node, G.number_of_nodes(), n, label, level, N_sub, setname, N_sub_set))
      grouplist.append((n,label,level,in_set,N_sub,N_sub_set))
  groups = pd.DataFrame({
	'Id': [Id for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'label': [label for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'level': [level for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'in_%s'%setname: [in_set for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'N_sub': [N_sub for Id,label,level,in_set,N_sub,N_sub_set in grouplist],
	'N_sub_%s'%setname: [N_sub_set for Id,label,level,in_set,N_sub,N_sub_set in grouplist]
	}).sort_values(by=['N_sub_%s'%setname, 'N_sub'], ascending=False)
  return(groups)

#############################################################################
if __name__=="__main__":
  PROG=os.path.basename(sys.argv[0])
  t0 = time.time()
  logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(
	description="NetworkX analytics, designed for EFO",
	epilog="clustering into groups by common ancestor")
  ops = ['summary', 'graph2cyjs', 'cluster']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--i_edge", dest="ifile_edge", help="input edgelist (TSV)")
  parser.add_argument("--i_node_attr", dest="ifile_node_attr", help="input node attributes (TSV)")
  parser.add_argument("--i_node_set", dest="ifile_node_set", help="input node set (IDs)")
  parser.add_argument("--setname", default="myset", help="node set name")
  parser.add_argument("--o", dest="ofile", help="output (TSV|CYJS|etc.)")
  parser.add_argument("--graphname", help="assign this name to graph")
  parser.add_argument("--min_groupsize", type=int, default=100)
  parser.add_argument("-v", "--verbose", action="count")
  args = parser.parse_args()

  #args.ifile_edge = "data/efo_edgelist.tsv"
  #args.ifile_node_attr = "data/efo_nodelist.tsv"
  #args.ifile_node_set = "data/gwascatalog.efoid"
  #args.ofile = "data/efo_groups.tsv"

  if args.ifile_edge:
    logging.info("Reading {0}".format(args.ifile_edge))
    efo_edges = pd.read_csv(args.ifile_edge, "\t", dtype=str)
    G = from_pandas_edgelist(efo_edges, source="source", target="target", edge_attr="edge_attr", create_using=nx.DiGraph)
    G.graph['name'] = args.graphname
  else:
    parser.error('--i_edge required.')

  if args.ofile:
    fout = open(args.ofile, "w")
  else:
    fout = sys.stdout

  ###
  if args.ifile_node_attr:
    logging.info("Reading {0}".format(args.ifile_node_attr))
    efo_nodes = pd.read_csv(args.ifile_node_attr, "\t", index_col="id")
    node_attr = efo_nodes.to_dict(orient='index')
    nx.set_node_attributes(G, node_attr)
  #
  ###
  if args.ifile_node_set:
    nodeSetIds=set()
    logging.info("Reading {0}".format(args.ifile_node_set))
    with open(args.ifile_node_set) as fin:
      for line in fin:
        nodeSetIds.add(line.strip())
    logging.info("nodeSetIds in {0}: {1}".format(args.setname, len(nodeSetIds)))
    nx.set_node_attributes(G, {nodeId:{'in_%s'%args.setname:True} for nodeId in nodeSetIds})

  ###
  if args.op == 'summary':
    GraphSummary(G)
  #
  elif args.op == 'graph2cyjs':
    #args.ofile = "data/efo_nxgraph.cyjs"
    logging.info("Writing {0}".format(args.ofile))
    Graph2CYJS(G, fout)
  #
  ###
  # Find ancestors with N_subclasses >= MIN_CLUSTER_SIZE.
  elif args.op == 'cluster':
    if not args.ifile_node_set:
      parser.error('--i_node_set required for cluster operation.')
    groups = Cluster(G, nodeSetIds, args.min_groupsize, args.setname)
    Groups2TSV(groups, fout)
    #print(groups.head(10))
    print(groups[groups['in_%s'%setname]].head(18)) #DEBUG

  logging.info(('%s: elapsed time: %s'%(PROG, time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)))))


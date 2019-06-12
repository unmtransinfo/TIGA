#!/usr/bin/env python3
###
import sys,os,argparse,re
import pandas

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(
        description='Pandas utilities for simple datafile transformations.')
  ops = ['csv2tsv', 'summary','showcols','selectcols','uvalcounts','colvalcounts','sortbycols','deduplicate']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--i", dest="ifile", help="input (CSV|TSV)")
  parser.add_argument("--o", dest="ofile", help="output (CSV|TSV)")
  parser.add_argument("--coltags", help="cols specified by tag (comma-separated)")
  parser.add_argument("--cols", help="cols specified by idx (1+) (comma-separated)")
  parser.add_argument("-v", "--verbose", action="count")
  args = parser.parse_args()

  if args.op in ('selectcols', 'uvalcounts', 'colvalcounts', 'sortbycols'):
    if not (args.cols or args.coltags): 
      parser.error('%s requires --cols or --coltags.'%args.op)

  if not args.ifile:
    parser.error('Input file required.')
  fin = open(args.ifile)

  if args.ofile:
    fout = open(args.ofile, "w")
  else:
    fout = sys.stdout

  if args.ifile[-4:].lower()=='.csv': delim=','
  else: delim='\t'

  if args.op == 'csv2tsv':
    df = pandas.read_table(fin, ',')
    df.to_csv(fout, '\t', index=False)
    exit(0)

  cols=None; coltags=None;
  if args.cols:
    cols = [(int(col.strip())-1) for col in re.split(r',', args.cols.strip())]
  elif args.coltags:
    coltags = [coltag.strip() for coltag in re.split(r',', args.coltags.strip())]

  df = pandas.read_table(fin, delim)

  if args.op == 'summary':
    print("rows: %d ; cols: %d"%(df.shape[0], df.shape[1]))
    print("coltags: %s"%(', '.join(['"%s"'%tag for tag in df.columns])))

  elif args.op == 'showcols':
    for j,tag in enumerate(df.columns):
      print('%d. "%s"'%(j+1,tag))

  elif args.op == 'selectcols':
    ds = df[coltags] if coltags else df.iloc[:, cols]
    ds.to_csv(fout, '\t', index=False)

  elif args.op == 'uvalcounts':
    for j,tag in enumerate(df.columns):
      if cols and j not in cols: continue
      if coltags and tag not in coltags: continue
      print('%d. %s: %d'%(j+1,tag,df[tag].nunique()))

  elif args.op == 'colvalcounts':
    for j,tag in enumerate(df.columns):
      if cols and j not in cols: continue
      if coltags and tag not in coltags: continue
      print('%d. %s:'%(j+1, tag))
      for key,val in df[tag].value_counts().iteritems():
        print('\t%6d: %s'%(val, key))

  else:
    parser.error('Unknown operation: %s'%args.op)

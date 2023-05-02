#! /usr/bin/env python3
"""
Set management; set operations on files of strings, one per line
"""
import sys,os,argparse,re,logging
import pandas as pd

#############################################################################
if __name__=='__main__':
  EPILOG='''\
AminusB: A and not-B
AandB: intersection of A and B
'''
  parser = argparse.ArgumentParser(description='Set operations on files of strings, one per line', epilog=EPILOG)
  OPS = ['AminusB', 'AandB']
  parser.add_argument("op", choices=OPS, help='OPERATION')
  parser.add_argument("--iA", dest="ifileA", required=True, help="input file A")
  parser.add_argument("--iB", dest="ifileB", required=True, help="input file B")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("--numerical", action="store_true", help="numerical comparison")
  parser.add_argument("--sort", action="store_true", help="sort output")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  fout = open(args.ofile,"w") if args.ofile else sys.stdout

  logging.debug(f"Input file A: {args.ifileA}")
  dfA = pd.read_table(args.ifileA, sep='\t', header=None)
  logging.debug(f"FileA {args.ifileA} input lines: {dfA.shape[0]}")
  logging.debug(f"Input file B: {args.ifileB}")
  dfB = pd.read_table(args.ifileB, sep='\t', header=None)
  logging.debug(f"FileB {args.ifileB} input lines: {dfB.shape[0]}")

  listA = dfA[0].astype('float' if args.numerical else 'str').values
  setA = set(listA)
  logging.debug(f"FileA unique values: {len(setA)}")

  listB = dfB[0].astype('float' if args.numerical else 'str').values
  setB = set(listB)
  logging.debug(f"FileB unique values: {len(setB)}")

  logging.debug(f"A rows: {dfA.shape[0]}")
  logging.debug(f"A Unique entities: {len(setA)}")
  logging.debug(f"B rows: {dfB.shape[0]}")
  logging.debug(f"B Unique entities: {len(setB)}")
  logging.debug(f"Output file: {args.ofile}")

  if args.op=="AminusB":
    setOut = setA - setB
    df_out = pd.DataFrame({args.op: list(setOut)})
    if args.sort: df_out = df_out.sort_values(by=args.op)
    df_out.to_csv(fout, sep="\t", index=False, header=False)
    logging.debug(f"Output lines: {df_out.shape[0]}")
    logging.info(f"Unique entities in output: {len(setOut)}")

  elif args.op=="AandB":
    setOut = setA & setB
    df_out = pd.DataFrame({args.op: list(setOut)})
    if args.sort: df_out = df_out.sort_values(by=args.op)
    df_out.to_csv(fout, sep="\t", index=False, header=False)
    logging.debug(f"Unique entities in output: {len(setOut)} ({100*len(setOut)/len(setA):.1f} of A, {100*len(setOut)/len(setB):.1f} of B)")
    logging.debug(f"Output lines: {df_out.shape[0]}")
    logging.info(f"Unique entities in output: {len(setOut)}")

  else:
    parser.error(f"Invalid operation: {args.op}")


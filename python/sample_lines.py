#!/usr/bin/env python3
#############################################################################
import sys,os,argparse
import random

PROG=os.path.basename(sys.argv[0])

#############################################################################
def SampleKofN(ifile, k, n, verbose):
  fin=open(ifile)
  if not n:
    n=0
    while True:
      line=fin.readline()
      if not line: break
      n+=1
    fin.close()
    fin=open(ifile)
  print("n = %d"%(n), file=sys.stderr)
  i_sams=set(random.sample(range(n),k))
  n_in=0; n_out=0;
  while True:
    line=fin.readline()
    if not line: break
    if n_in in i_sams:
      sys.stdout.write(line)
      n_out+=1
      if n_out>=k: break
    n_in+=1
  print("n_in: %d"%(n_in), file=sys.stderr)
  print("n_out: %d"%(n_out), file=sys.stderr)

#############################################################################
def SampleP(ifile, p, verbose):
  if ifile=='-':
    fin = sys.stdin
  else:
    fin=open(ifile)
  n_in=0; n_out=0;
  while True:
    line=fin.readline()
    if not line: break
    n_in+=1
    if random.random()<p:
      sys.stdout.write(line)
      n_out+=1
  print("n_in: %d"%(n_in), file=sys.stderr)
  print("n_out: %d"%(n_out), file=sys.stderr)
  print("p: %.2f ; p_sample: %.2f"%(p,float(n_out)/n_in), file=sys.stderr)

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='randomly sample lines from file')
  ops = ['chooseK', 'sampleP']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--i", dest="ifile", help='input file, or "-" for stdin')
  parser.add_argument("--k", type=int, help='choose K lines')
  parser.add_argument("--n", type=int, help='population size [default=all-lines] (disallowed for stdin)')
  parser.add_argument("--p", type=float, help='sample with probability P each line')
  parser.add_argument("-v", "--verbose", default=0, action="count")
  args = parser.parse_args()

  if args.op == 'chooseK':
    if args.ifile=='-':
      parser.error('ERROR: operation %s disallowed with stdin.'%args.op)
    if not args.k:
      parser.error('ERROR: operation %s requires K.'%args.op)
    SampleKofN(args.ifile, args.k, args.n, args.verbose)

  elif args.op == 'sampleP':
    if not args.p:
      parser.error('ERROR: operation %s requires P.'%args.op)
    SampleP(args.ifile, args.p, args.verbose)

  else:
    parser.error('ERROR: unknown operation: %s'%args.op)


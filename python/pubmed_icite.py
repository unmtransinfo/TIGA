#!/usr/bin/env python3
#############################################################################
### PUbMed iCite  REST API client
###
### https://icite.od.nih.gov/api
#############################################################################
### Jeremy Yang
#############################################################################
import sys,os,argparse
import urllib.request
import json,re
#
API_HOST="icite.od.nih.gov"
API_BASE_PATH="/api/pubs"
API_BASE_URL='https://'+API_HOST+API_BASE_PATH
#
#############################################################################
def LOG(msg, file=sys.stdout):
  print(msg, file=file, flush=True)
#
#############################################################################
def GetPmid(base_url,pmids,fout,verbose):
  n_out=0;
  tags=None;
  for pmid in pmids:
    url=base_url+'/%s'%pmid
    f = urllib.request.urlopen(url)
    rval=f.read()
    if not rval:
      if verbose:
        LOG('not found: %s'%(pmid))
      continue
    rval=json.loads(rval.decode('utf-8'), encoding='utf-8')
    pub = rval
    if not tags:
      tags = pub.keys()
      fout.write('\t'.join(tags)+'\n')

    vals = [];
    for tag in tags:
      val=(pub[tag] if pub.has_key(tag) else '')
      vals.append(val)
    fout.write('\t'.join([str(vals) for val in vals])+'\n')
    n_out+=1

  LOG('n_in = %d'%(len(pmids)))
  LOG('n_out = %d'%(n_out))

#############################################################################
def GetPmids(base_url, pmids, fout, verbose):
  n_in=0; n_out=0; tags=None; nchunk=100;

  while True:
    if n_in >= len(pmids):
      break
    pmids_this = pmids[n_in:n_in+nchunk]
    n_in += (nchunk if n_in+nchunk < len(pmids) else len(pmids)-n_in)
    url_this=base_url+'?pmids=%s'%(','.join(pmids_this))
    f = urllib.request.urlopen(url_this)
    rval=f.read()

    if not rval:
      break

    rval=json.loads(rval.decode('utf-8'), encoding='utf-8')
    #LOG('DEBUG: rval="%s"'%rval.decode('utf-8'))

    if verbose:
      LOG('%s'%url_this)

    url_self = rval['links']['self']
    #if verbose: print >>sys.stderr, 'DEBUG: %s'%url_self

    pubs = rval['data']
    for pub in pubs:
      if not tags:
        tags = pub.keys()
        fout.write('\t'.join(tags)+'\n')

      vals=[];
      for tag in tags:
        val=(pub[tag] if tag in pub else '')
        vals.append(val)
      fout.write('\t'.join([str(val) for val in vals])+'\n')
      n_out+=1

  LOG('n_in = %d'%(len(pmids)))
  LOG('n_out = %d'%(n_out))

#############################################################################
if __name__=='__main__':

  parser = argparse.ArgumentParser(
	description='PubMed iCite REST API client utility',
	epilog='Publication metadata.')
  ops = ['get']
  parser.add_argument("op",choices=ops,help='operation')
  parser.add_argument("--pmid",dest="pmid",help="PubMed ID (ex:25533513)")
  parser.add_argument("--pmids",dest="pmids",help="PubMed IDs, comma-separated")
  parser.add_argument("--i",dest="ifile",help="input file, PubMed IDs")
  parser.add_argument("--nmax",help="list: max to return")
  parser.add_argument("--year",help="list: year of publication")
  parser.add_argument("--o",dest="ofile",help="output (TSV)")
  parser.add_argument("-v","--verbose",action="count")

  args = parser.parse_args()

  if args.ofile:
    fout=open(args.ofile,"w",encoding="utf-8")
    if not fout:
      parser.error('ERROR: cannot open outfile: %s'%args.ofile)
      parser.exit()
  else:
    fout=sys.stdout

  pmids=[];
  pmid=None;
  if args.ifile:
    fin=open(args.ifile)
    if not fin:
      parser.error('ERROR: cannot open ifile: %s'%args.ifile)
      parser.exit()
    while True:
      line=fin.readline()
      if not line: break
      pmids.append(line.rstrip())
    if args.verbose:
      LOG('input IDs: %d'%(len(pmids)))
    fin.close()
  elif args.pmids:
    pmids = re.split(r'\s*,\s*',args.pmids.strip())
  elif args.pmid:
    pmids.append(args.pmid)
  pmid=pmids[0]

  if args.op == 'get':
    if not pmids:
      parser.error('get requires PMID[s].')
      parser.exit()
    GetPmids(API_BASE_URL, pmids, fout, args.verbose)

  else:
    parser.print_help()


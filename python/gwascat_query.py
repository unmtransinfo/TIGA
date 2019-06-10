#!/usr/bin/env python3
"""
	GWAS Catalog REST API client.

	https://www.ebi.ac.uk/gwas/docs/api
	https://www.ebi.ac.uk/gwas/rest/docs/api
	https://www.ebi.ac.uk/gwas/rest/docs/sample-scripts
	https://www.ebi.ac.uk/gwas/rest/api

	author: Jeremy Yang
"""
import sys,os,re,argparse,json
import urllib.parse
#
import rest_utils
#
PROG=os.path.basename(sys.argv[0])
#
API_HOST='www.ebi.ac.uk'
API_BASE_PATH='/gwas/rest/api'
#
# https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=28530673
# https://www.ebi.ac.uk/gwas/rest/api/studies/GCST004364
##############################################################################
def SearchStudies(base_url, ids, idtype, fout, verbose):
  url=base_url+'/studies/search'
  if idtype=='pmid':
    url+='/findByPublicationIdPubmedId?pubmedId='
  elif idtype=='gcst':
    url+='/'
  else:
    print('ERROR: Not yet supported: idtype = %s'%(idtype), file=sys.stderr)
    return

  for id_this in ids:
    url_this=url+'%s'%id_this
    rval=rest_utils.GetURL(url_this, parse_json=True, verbose=verbose)
    if not rval:
      continue
    print('DEBUG: %s'%json.dumps(rval,sort_keys=True,indent=2))

##############################################################################
# https://www.ebi.ac.uk/gwas/rest/api/studies/GCST004364/associations?projection=associationByStudy
# https://www.ebi.ac.uk/gwas/rest/api/studies/GCST000227/associations?projection=associationByStudy
### We know that GCST000227 has mapped genes, so why only reported in returned loci data?

def SearchAssociations(base_url, ids, idtype, fout, verbose):
  url=base_url+'/studies'

  if idtype=='gcst':
    url+='/'
  else:
    print('ERROR: Not yet supported: idtype = %s'%(idtype), file=sys.stderr)
    return

  n_assn=0; n_loci=0; tags=[]; tags_loci=[];
  for id_this in ids:
    url_this=url+'%s/associations?projection=associationByStudy'%id_this
    rval=rest_utils.GetURL(url_this, parse_json=True, verbose=verbose)
    if not rval:
      continue
    if '_embedded' in rval and 'associations' in rval['_embedded']:
      assns = rval['_embedded']['associations']
    else:
      print('ERROR: no associations.', file=sys.stderr)
      return

    for assn in assns:
      #print('DEBUG: %s'%json.dumps(assn, sort_keys=True,indent=2), file=sys.stderr)
      n_assn+=1
      if n_assn==1:
        for key,val in assn.items():
          if type(val) not in (list, dict):
            tags.append(key)
        for key,val in assn['loci'][0].items():
          #if type(val) not in (list, dict):
          tags_loci.append(key)

        fout.write('\t'.join(tags+list(map(lambda s: 'loci_'+s, tags_loci)))+'\n')

      vals=[]; vals_loci=[];
      for tag in tags:
        if tag in assn:
          vals.append(str(assn[tag]))
      for loci in assn['loci']:
        n_loci+=1
        for tag in tags_loci:
          if tag in assn['loci']:
            vals_loci.append(str(assn['loci'][tag]))

        fout.write('\t'.join(vals+vals_loci)+'\n')

  print('associations: %d ; loci: %d'%(n_assn, n_loci), file=sys.stderr)

##############################################################################
# https://www.ebi.ac.uk/gwas/rest/api/studies/GCST004364?projection=study

##############################################################################
if __name__=='__main__':

  idtypes=['pmid', 'gcst', 'efo', 'rs']
  parser = argparse.ArgumentParser(
        description='GWAS Catalog query client')
  ops = ['searchStudies', 'searchAssociations']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--ids", dest="ids", help="IDs, comma-separated")
  parser.add_argument("--idtype", choices=idtypes, default="pmid", help="ID type")
  parser.add_argument("--i", dest="ifile", help="input file, IDs")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("--api_host", default=API_HOST)
  parser.add_argument("--api_base_path", default=API_BASE_PATH)
  parser.add_argument("-v", "--verbose", default=0, action="count")
  args = parser.parse_args()

  base_url='https://'+args.api_host+args.api_base_path

  if args.ofile:
    fout=open(args.ofile, 'w')
    if not fout: parser.error('ERROR: cannot open outfile: %s'%args.ofile)
  else:
    fout=sys.stdout

  ids=[]
  if args.ifile:
    fin=open(args.ifile)
    if not fin: parser.error('ERROR: failed to open input file: %s'%args.ifile)
    while True:
      line=fin.readline()
      if not line: break
      ids.append(line.strip())
  elif args.ids:
    ids=re.split(r'\s*,\s*',args.ids.strip())

  if args.op == 'searchStudies':
    SearchStudies(base_url, ids, args.idtype, fout, args.verbose)

  elif args.op == 'searchAssociations':
    SearchAssociations(base_url, ids, args.idtype, fout, args.verbose)

  else:
    parser.error("Unknown operation: %s"%args.op)

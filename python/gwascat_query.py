#!/usr/bin/env python3
"""
	GWAS Catalog REST API client.

	https://www.ebi.ac.uk/gwas/docs/api
	https://www.ebi.ac.uk/gwas/rest/docs/api
		See "Response structure", "Links", etc.
		&page=1&size=50 (default size=20, max 500)
		http://stateless.co/hal_specification.html
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
def StudyAssociations(base_url, ids, fout, verbose):
  """
	/studies/GCST000227/associations?projection=associationByStudy
	Mapped genes via SNP links.
	arg = authorReportedGene
	sra = strongestRiskAllele
  """
  url=base_url+'/studies'
  idtype='gcst'
  url+='/'

  n_id=0;
  n_assn=0; n_loci=0; n_arg=0; n_sra=0; n_snp=0;
  tags=[]; tags_locus=[]; tags_sra=[]; tags_arg=[];
  for id_this in ids:
    n_id+=1
    url_this=url+'%s/associations?projection=associationByStudy'%id_this
    rval=rest_utils.GetURL(url_this, parse_json=True, verbose=verbose)
    if not rval:
      continue
    if '_embedded' in rval and 'associations' in rval['_embedded']:
      assns = rval['_embedded']['associations']
    else:
      print('ERROR: no associations for study: %s'%id_this, file=sys.stderr)
      continue

    for assn in assns:
      n_assn+=1
      if n_assn==1:
        for key,val in assn.items():
          if type(val) not in (list, dict):
            tags.append(key)
        for key in assn['loci'][0].keys():
          if key not in ('strongestRiskAlleles', 'authorReportedGenes'):
            tags_locus.append(key)
        for key in assn['loci'][0]['strongestRiskAlleles'][0].keys():
          if key != '_links':
            tags_sra.append(key)
        for key in assn['loci'][0]['authorReportedGenes'][0].keys():
          tags_arg.append(key)

        fout.write('\t'.join(tags
		+list(map(lambda s: 'locus_'+s, tags_locus))
		+list(map(lambda s: 'reported_'+s, tags_arg))
		+list(map(lambda s: 'allele_'+s, tags_sra))
		+['snp_url'])+'\n')

      vals=[];
      for tag in tags:
        vals.append(str(assn[tag]) if tag in assn and assn[tag] is not None else '')
      for locus in assn['loci']:
        n_loci+=1
        vals_locus=[];
        for tag in tags_locus:
          vals_locus.append(str(locus[tag]) if tag in locus and locus[tag] is not None else '')
        for arg in locus['authorReportedGenes']:
          n_arg+=1
          vals_arg=[];
          for tag in tags_arg:
            if tag=='entrezGeneIds':
              geneids = [g['entrezGeneId'] for g in arg[tag]]
              if None in geneids: geneids.remove(None)
              vals_arg.append(';'.join(geneids) if geneids else '')
            elif tag=='ensemblGeneIds':
              geneids = [g['ensemblGeneId'] for g in arg[tag]]
              if None in geneids: geneids.remove(None)
              vals_arg.append(';'.join(geneids) if geneids else '')
            else:
              vals_arg.append(str(arg[tag]) if tag in arg else '')
          fout.write('\t'.join(vals+vals_locus+vals_arg+['' for tag in tags_sra]+[''])+'\n')

        for sra in locus['strongestRiskAlleles']:
          n_sra+=1
          vals_sra=[]; snp_href='';
          snp_href = sra['_links']['snp']['href'] if '_links' in sra and 'snp' in sra['_links'] and 'href' in sra['_links']['snp'] else ''
          if snp_href: n_snp+=1
          for tag in tags_sra:
            vals_sra.append(str(sra[tag]) if tag in sra and sra[tag] is not None else '')
          fout.write('\t'.join(vals+vals_locus+['' for tag in tags_arg]+vals_sra+[snp_href])+'\n')

  print('%ss: %d ; assns: %d ; loci: %d ; reportedGenes: %d ; alleles: %d ; snps: %d'%(idtype, n_id, n_assn, n_loci, n_arg, n_sra, n_snp), file=sys.stderr)

##############################################################################
def GetSnps(base_url, ids, fout, verbose):
  """
	/singleNucleotidePolymorphisms/rs6085920
	loc = location
	gc = genomicContext
  """
  url=base_url+'/singleNucleotidePolymorphisms'+'/'

  n_snp=0; n_gc=0; n_gene=0; n_gcloc=0; n_loc=0;
  tags=[]; tags_loc=[]; tags_gc=[]; tags_gcloc=[];  tags_gene=[]; 
  for id_this in ids:
    n_snp+=1
    url_this=url+id_this
    snp=rest_utils.GetURL(url_this, parse_json=True, verbose=verbose)
    if not snp:
      continue
    if n_snp==1:
      for key,val in snp.items():
        if type(val) not in (list, dict):
          tags.append(key)
      #for key in snp['locations'][0].keys():
      #  if key != '_links':
      #    tags_loc.append(key)
      for key in snp['genomicContexts'][0].keys():
        if key not in ('gene', 'location', '_links'):
          tags_gc.append(key)
      for key in snp['genomicContexts'][0]['location'].keys():
        if key != '_links':
          tags_gcloc.append(key)
      for key in snp['genomicContexts'][0]['gene'].keys():
        if key != '_links':
          tags_gene.append(key)

      fout.write('\t'.join(tags
		+list(map(lambda s: 'genomicContext_'+s, tags_gc))
		+list(map(lambda s: 'loc_'+s, tags_gcloc))
		+list(map(lambda s: 'gene_'+s, tags_gene)))
		+'\n')

    vals=[];
    for tag in tags:
      vals.append(str(snp[tag]) if tag in snp and snp[tag] is not None else '')
    for gc in snp['genomicContexts']:
      n_gc+=1
      vals_gc=[];
      for tag in tags_gc:
        vals_gc.append(str(gc[tag]) if tag in gc and gc[tag] is not None else '')
      gcloc = gc['location']
      n_gcloc+=1
      vals_gcloc=[];
      for tag in tags_gcloc:
        if tag=='region':
          vals_gcloc.append(gcloc[tag]['name'] if tag in gcloc and gcloc[tag] and 'name' in gcloc[tag] and gcloc[tag]['name'] is not None else '')
        else:
          vals_gcloc.append(str(gcloc[tag]) if tag in gcloc and gcloc[tag] is not None else '')
      gene = gc['gene']
      n_gene+=1
      vals_gene=[];
      for tag in tags_gene:
        if tag=='entrezGeneIds':
          geneids = [g['entrezGeneId'] for g in gene[tag]]
          if None in geneids: geneids.remove(None)
          vals_gene.append(';'.join(geneids) if geneids else '')
        elif tag=='ensemblGeneIds':
          geneids = [g['ensemblGeneId'] for g in gene[tag]]
          if None in geneids: geneids.remove(None)
          vals_gene.append(';'.join(geneids) if geneids else '')
        else:
          vals_gene.append(str(gene[tag]) if tag in gene and gene[tag] is not None else '')

      fout.write('\t'.join(vals+vals_gc+vals_gcloc+vals_gene)+'\n')

  print('SNPs: %d ; genomicContexts: %d ; genes: %d ; locations: %d'%(n_snp, n_gc, n_gene, n_gcloc), file=sys.stderr)

##############################################################################
# https://www.ebi.ac.uk/gwas/rest/api/studies/GCST004364?projection=study

##############################################################################
if __name__=='__main__':

  idtypes=['pmid', 'gcst', 'efo', 'rs']
  parser = argparse.ArgumentParser(
        description='GWAS Catalog query client')
  ops = ['searchStudies', 'studyAssociations', 'getSnps']
  parser.add_argument("op", choices=ops, help='operation')
  parser.add_argument("--ids", dest="ids", help="IDs, comma-separated")
  parser.add_argument("--idtype", choices=idtypes, help="ID type")
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

  elif args.op == 'studyAssociations':
    if args.idtype and args.idtype != 'gcst':
      parser.error('Operation "%s" requires idtype "gcst"'%args.op)
    StudyAssociations(base_url, ids, fout, args.verbose)

  elif args.op == 'getSnps':
    if args.idtype and args.idtype != 'rs':
      parser.error('Operation "%s" requires idtype "rs"'%args.op)
    GetSnps(base_url, ids, fout, args.verbose)

  else:
    parser.error("Unknown operation: %s"%args.op)

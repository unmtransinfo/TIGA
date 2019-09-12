#!/usr/bin/env python3
"""
	From GWAS Catalog associations, generate input files for VEGAS2,
	one file for each study, two columns SNP-id and p-value.
	Input file gwascat_assn.tsv is pre-cleaned by gwascat_assn.R
	from raw GWAS Catalog file, e.g.,
	gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv.
"""
import sys,os,argparse
import pandas

#############################################################################
def WriteVegasInputFiles(df, odir, prefix, verbose):


#############################################################################
if __name__=='__main__':

  parser = argparse.ArgumentParser(
	description='GWAS Catalog to VEGAS2 input files.')
  parser.add_argument("--i", dest="ifile", help="input (TSV)")
  parser.add_argument("--odir", default=".", help="output dir")
  parser.add_argument("--prefix", default="vegas_", help="output prefix")
  parser.add_argument("--study_id_tag", default="STUDY_ACCESSION")
  parser.add_argument("--snp_id_tag", default="SNP_ID_CURRENT")
  parser.add_argument("--pval_tag", default="P-VALUE")
  parser.add_argument("-v", "--verbose", action="count")
  args = parser.parse_args()

  if not args.ifile:
    parser.error('Input file required.')

  df = pandas.read_csv(args.ifile, sep='\t')

  if args.verbose>0:
    print("rows: %d ; cols: %d"%(df.shape[0], df.shape[1]))
    print("coltags: %s"%(', '.join(['"%s"'%tag for tag in df.columns])))
    for j,tag in enumerate(df.columns):
      print('%d. "%s"'%(j+1,tag))

  WriteVegasInputFiles(df, args.odir, args.prefix, args.verbose)


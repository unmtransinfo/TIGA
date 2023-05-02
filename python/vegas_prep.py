#!/usr/bin/env python3
"""
	From GWAS Catalog associations, generate input files for VEGAS2,
	one file for each study, two columns SNP-id and p-value.
	Input file gwascat_assn.tsv is pre-cleaned by gwascat_assn.R
	from raw GWAS Catalog file, e.g.,
	gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv.
"""
import sys,os,os.path,argparse
import pandas

#############################################################################
def WriteVegasInputFiles(df, args):
  accs = df[args.study_id_tag].unique()
  i_acc=0; i_pvs=0;
  snps = set([])
  for acc in accs:
    i_acc+=1
    print("%d/%d. %s"%(i_acc, len(accs), acc), file=sys.stderr)
    df_this = df[df[args.study_id_tag] == acc]
    df_this = df_this[[args.snp_id_tag, args.pval_tag]]
    df_this[args.pval_tag] = df_this[args.pval_tag].astype('float')
    if sum(df_this[args.pval_tag].isna())>0:
      df_this = df_this[df_this[args.pval_tag].notna()]
      print("NAs removed: %d"%(sum(df_this[args.pval_tag].isna())), file=sys.stderr)
    if df_this.shape[0]==0:
      print("Zero p-values; skipping (%s)."%(acc), file=sys.stderr)
      next
    i_pvs+=df_this.shape[0]
    snps |= set(df_this[args.snp_id_tag])
    if args.test:
      fout = sys.stdout
    else:
      fout = open(args.odir+'/'+args.prefix+acc+'.tsv', 'w')
    df_this.to_csv(fout, sep='\t', index=False, header=False)
  print("Total studies (%s): %d"%(args.study_id_tag, len(accs)), file=sys.stderr)
  print("Total files output: %d"%(i_acc), file=sys.stderr)
  print("Total p-values output: %d"%(i_pvs), file=sys.stderr)
  print("Total (unique) SNPs output: %d"%(len(snps)), file=sys.stderr)

#############################################################################
if __name__=='__main__':

  parser = argparse.ArgumentParser(
	description='GWAS Catalog to VEGAS2 input files.')
  parser.add_argument("--i", dest="ifile", help="input (TSV)")
  parser.add_argument("--odir", help="output dir")
  parser.add_argument("--prefix", default="vegas_in_", help="output prefix")
  parser.add_argument("--study_id_tag", default="STUDY_ACCESSION")
  parser.add_argument("--snp_id_tag", default="SNPS")
  parser.add_argument("--pval_tag", default="P-VALUE")
  parser.add_argument("--test", help="test-mode, no output files")
  parser.add_argument("-v", "--verbose", default=0, action="count")
  args = parser.parse_args()

  if not args.ifile:
    parser.error('Input file required.')

  if not args.odir:
    parser.error('Output directory required.')
  else:
    if not os.path.exists(args.odir):
      parser.error('Output directory exists not: %s'%args.odir)

  df = pandas.read_csv(args.ifile, sep='\t', low_memory=False)

  if args.verbose:
    print("rows: %d ; cols: %d"%(df.shape[0], df.shape[1]), file=sys.stderr)
    print("coltags: %s"%(', '.join(['"%s"'%tag for tag in df.columns])), file=sys.stderr)
    print("studies: %d"%(df[args.study_id_tag].nunique()), file=sys.stderr)

  WriteVegasInputFiles(df, args)


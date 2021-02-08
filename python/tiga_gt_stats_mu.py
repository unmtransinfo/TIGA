#!/usr/bin/env python3
"""
Compute MU scores from selected gene-trait variables.
Input: gt_variables.tsv.gz
Output: gt_stats_mu.tsv.gz
Replaces tiga_gt_stats_mu.R, which uses muStat and is a memory hog (50+G for 50k*3 matrix)
This code may be slow but runs with typical workstation memory.
"""
###
import sys,os,re,argparse,time,logging,tqdm
import pandas as pd
import numpy as np

#############################################################################
def ComputeMuScores(df, coltags, fout):
  """Variables larger for higher rank. Pre-convert NaNs to values.""" 
  df = df[coltags].astype('float')
  for tag in df.columns.tolist():
    if df[tag].isna().sum()>0:
      logging.error("NaNs found in column: \"{tag}\"")
      return
  df['nAbove'] = np.nan
  df['nBelow'] = np.nan
  tq = tqdm.tqdm(total=df.shape[0], unit="rows")
  for i in range(df.shape[0]):
    tq.update()
    vals_this = [df[tag].iloc[i] for tag in coltags]
    logging.debug(f"{i}. vals_this: {vals_this}")
    ges_this = pd.Series([True for i in range(df.shape[0])])
    lts_this = pd.Series([True for i in range(df.shape[0])])
    for j in range(len(vals_this)):
      ges_this = (ges_this & (df.iloc[:,j] >= vals_this[j]))
      lts_this = (lts_this & (df.iloc[:,j] <  vals_this[j]))
    logging.debug(f"{i}. sum(ges_this): {sum(ges_this)}")
    logging.debug(f"{i}. sum(lts_this): {sum(lts_this)}")
    df['nBelow'].iloc[i] = sum(ges_this)
    df['nAbove'].iloc[i] = sum(lts_this)
    #if i>1000: break #DEBUG
  tq.close()
  df['muScore'] = df['nBelow'] - df['nAbove'] 
  df.to_csv(fout, "\t", index=False)
  logging.info(f"nBelow range: [{min(df['nBelow'])},{max(df['nBelow'])}]")
  logging.info(f"nAbove range: [{min(df['nAbove'])},{max(df['nAbove'])}]")
  logging.info(f"muScore range: [{min(df['muScore'])},{max(df['muScore'])}]")
  logging.info(f"MU scores computed for {df.shape[0]} rows, using variables: {coltags}")

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Compute MU scores from selected gene-trait variables.", epilog="")
  parser.add_argument("--i", required=True, dest="ifile", help="Input gene-trait variables (TSV)")
  parser.add_argument("--o", dest="ofile", help="output gene-trait MU stats (TSV)")
  parser.add_argument("--coltags", required=True, help="Selected columns for multivariate MU scoring (comma-separated).")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG if args.
verbose>0 else logging.INFO)

  t0 = time.time()

  delim = ',' if re.search('\.csv', args.ifile, re.I) else '\t' if re.search('\.(tsv|tab|txt)', args.ifile, re.I) else '\t'

  fout = open(args.ofile, "w") if args.ofile else sys.stdout

  df = pd.read_csv(args.ifile, sep=delim)

  coltags = [coltag.strip() for coltag in re.split(r',', args.coltags.strip())]

  for tag in coltags:
    if tag not in df.columns:
      logging.error(f"Column not found: \"{tag}\"; columns present: {df.columns}")
      sys.exit()

  ComputeMuScores(df, coltags, fout)

  logging.info(('Elapsed time: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)))))


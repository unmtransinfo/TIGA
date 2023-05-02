#!/usr/bin/env python3
"""
Compute MU scores from selected gene-trait variables.
Input: gt_variables.tsv.gz
Output: gt_stats_mu.tsv.gz
Replaces tiga_gt_stats_mu.R, muStat memory hog (50+G for 50k*3 matrix)
This code may be slow but runs with typical workstation memory.
"""
###
import sys,os,re,argparse,time,logging,tqdm
import pandas as pd
import numpy as np

#############################################################################
def ComputeMuScores(df, mutags, ofile):
  """Variables larger for higher rank. Pre-convert NaNs to values.""" 
  quiet = bool(logging.getLogger().getEffectiveLevel()>15)
  df = df.astype({tag:'float' for tag in mutags})
  for tag in mutags:
    if df[tag].isna().sum()>0:
      logging.error("NaNs found in column: \"{tag}\"")
      return
  df['nAbove'] = np.nan
  df['nBelow'] = np.nan
  if not quiet: tq = tqdm.tqdm(total=df.shape[0], unit="rows")
  for i in range(df.shape[0]):
    if not quiet: tq.update()
    vals_this = {tag:df[tag].iloc[i] for tag in mutags}
    logging.debug(f"{i}. vals_this: {vals_this}")
    ges_this = pd.Series([True for i in range(df.shape[0])])
    lts_this = pd.Series([True for i in range(df.shape[0])])
    for tag in mutags:
      ges_this = (ges_this & (df.loc[:,tag] >= vals_this[tag]))
      lts_this = (lts_this & (df.loc[:,tag] <  vals_this[tag]))
    logging.debug(f"{i}. sum(ges_this): {sum(ges_this)}")
    logging.debug(f"{i}. sum(lts_this): {sum(lts_this)}")
    df.loc[i, 'nBelow'] = sum(ges_this)
    df.loc[i, 'nAbove'] = sum(lts_this)
    #if i==10: break #DEBUG
  if not quiet: tq.close()
  df['muScore'] = df['nBelow'] - df['nAbove'] 
  df.to_csv(ofile, sep="\t", index=False)
  logging.info(f"nBelow range: [{min(df['nBelow'])},{max(df['nBelow'])}]")
  logging.info(f"nAbove range: [{min(df['nAbove'])},{max(df['nAbove'])}]")
  logging.info(f"muScore range: [{min(df['muScore'])},{max(df['muScore'])}]")
  logging.info(f"MU scores computed for {df.shape[0]} rows, using variables: {mutags}")
  return df

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Compute MU scores from selected gene-trait variables.", epilog="")
  parser.add_argument("--i", dest="ifile", required=True, help="Input gene-trait variables (TSV)")
  parser.add_argument("--o", dest="ofile", required=True, help="Output gene-trait MU stats (TSV)")
  parser.add_argument("--mutags", required=True, help="Selected columns for multivariate MU scoring (comma-separated).")
  parser.add_argument("-q", "--quiet", action="store_true", help="Suppress progress notification.")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  # logging.PROGRESS = 15 (custom)
  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>0 else logging.ERROR if args.quiet else 15))

  logging.info(os.path.basename(sys.argv[0]))
  t0 = time.time()

  delim = ',' if re.search('\.csv', args.ifile, re.I) else '\t' if re.search('\.(tsv|tab|txt)', args.ifile, re.I) else '\t'

  df = pd.read_csv(args.ifile, sep=delim)

  mutags = [coltag.strip() for coltag in re.split(r'[,\s]', args.mutags.strip())]

  for tag in mutags:
    if tag not in df.columns:
      logging.error(f"Column not found: \"{tag}\"; columns present: {df.columns}")
      sys.exit()

  ComputeMuScores(df, mutags, args.ofile)

  logging.info(f"""{os.path.basename(sys.argv[0])} elapsed time: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}""")


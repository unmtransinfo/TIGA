#!/usr/bin/env python3
"""
Compare versions of TIGA output.
Inputs: gt_statsA.tsv.gz, gt_statsB.tsv.gz
"""
###
import sys,os,re,argparse,time,logging,tqdm
import pandas as pd
import numpy as np

#############################################################################
def CompareGeneTraitStats(ifileA, ifileB, entity, ofile):
  delimA = ',' if re.search('\.csv', ifileA, re.I) else '\t' if re.search('\.(tsv|tab|txt)', ifileA, re.I) else '\t'
  delimB = ',' if re.search('\.csv', ifileB, re.I) else '\t' if re.search('\.(tsv|tab|txt)', ifileB, re.I) else '\t'
  dfA = pd.read_csv(ifileA, sep=delimA)
  dfB = pd.read_csv(ifileB, sep=delimB)

  logging.debug(f"A: GENE ensemblIds: {dfA['ensemblId'].nunique()}")
  logging.debug(f"A: TRAIT efoIds: {dfA['efoId'].nunique()}")
  logging.debug(f"A: GENE-TRAIT pairs: {dfA[['ensemblId','efoId']].drop_duplicates().shape[0]}")
  logging.debug(f"B: GENE ensemblIds: {dfB['ensemblId'].nunique()}")
  logging.debug(f"B: TRAIT efoIds: {dfB['efoId'].nunique()}")
  logging.debug(f"B: GENE-TRAIT pairs: {dfB[['ensemblId','efoId']].drop_duplicates().shape[0]}")

  AminusB={};
  for tag in ['ensemblId', 'efoId']:
    AminusB[tag] = set(dfA[tag].unique()) - set(dfB[tag].unique())
    print(f"AminusB: {tag}s: {len(AminusB[tag])}")

  fout = open(ofile, "w") if ofile else sys.stdout

  if entity=="gene":
    df_out = dfA[["ensemblId", "geneSymbol", "geneName"]][dfA["ensemblId"].isin(AminusB["ensemblId"])].drop_duplicates()
    df_out.sort_values(by=["ensemblId"], inplace=True)
    df_out.to_csv(fout, "\t", index=False)
  elif entity=="trait":
    df_out = dfA[["efoId", "trait"]][dfA["efoId"].isin(AminusB["efoId"])].drop_duplicates()
    df_out.sort_values(by=["efolId"], inplace=True)
    df_out.to_csv(fout, "\t", index=False)
  else:
    df_out=None

  logging.info(f"N_{entity}: {df_out.shape[0]}")
  return df_out

#############################################################################
def ExplainMissing(df_missing, filterfile, entity):
  """
gene: "ensemblId", "ensemblSymb", "geneName", "geneFamily", "TDL", "reason"
trait: "TRAIT_URI", "TRAIT", "reason"
"""

  df_filter = pd.read_csv(filterfile, sep="\t")
  if entity=="gene":
    df_filter = df_filter[["ensemblId", "reason"]]
  elif entity=="trait":
    df_filter["efoId"] = df_filter["TRAIT_URI"].str.sub(r"^.*/", "")
    df_filter = df_filter[["efoId", "reason"]]
  else:
    logging.error(f"Invalid entity: {entity}")

  df_missing = pd.merge(df_missing, df_filter, how="left", on=("ensemblId" if entity=="gene" else "efoId"))

  counts = df_missing["reason"].value_counts(sort=True, dropna=False)
  for i in range(counts.size):
    print(f'''{counts.values[i]:6d}: "{counts.index[i]}"''')

#############################################################################
if __name__=="__main__":
  ENTITIES = ["gene", "trait"]
  parser = argparse.ArgumentParser(description="Compare GENE or TRAIT coverage for versions of TIGA output (AminusB).", epilog="")
  parser.add_argument("entity", choices=ENTITIES, help='ENTITY')
  parser.add_argument("--iA", dest="ifileA", required=True, help="(Version_A) Input gene-trait stats (TSV)")
  parser.add_argument("--iB", dest="ifileB", required=True, help="(Version_B) Input gene-trait stats (TSV)")
  parser.add_argument("--filterfile", help="(Version_B) GENE|TRAIT filter file, with reasons.")
  parser.add_argument("--explain_missing", action="store_true", help="Explain missing entities using filter file.")
  parser.add_argument("--o", dest="ofile", help="(TSV)")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG if args.verbose>0 else logging.INFO)

  df_missing = CompareGeneTraitStats(args.ifileA, args.ifileB, args.entity, args.ofile)

  if args.explain_missing:
    if not args.filterfile:
      parser.error(f"--explain_missing requires --filterfile.")
    ExplainMissing(df_missing, args.filterfile, args.entity)

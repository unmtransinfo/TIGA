#!/usr/bin/env python3
###
import sys,os,re,logging,argparse
import pandas as pd
import numpy as np

###
# Regarding rs (RefSNP) IDs, see
# https://www.ncbi.nlm.nih.gov/snp/docs/RefSNP_about/
# "dbSNP Reference SNP (rs or RefSNP) number is a locus accession for a
# variant type assigned by dbSNP."
###
# dbSNP API:
# https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/328
###
# Despite its name, RefSNP is assigned to all variation types listed below with precise locations for both common and rare variations, including mutations. Most are typically small variations (<= 50bp).
# 
# * Single nucleotide variation (SNV)
# * Short multi-nucleotide changes (MNV)
# * Small deletions or insertions
# * Small STR repeats
# * retrotransposable element insertions
###

# What about non-RefSNP (non-rs) SNPs?

###
# Input gwascat_assn.tsv, from gwascat_assn.R
###

if __name__=="__main__":
  parser = argparse.ArgumentParser(description='GWAS Catalog SNP2GENE links')
  parser.add_argument("INPUT_GWASCAT_ASSN_FILE", help="Cleaned associations, from gwascat_assn.R")
  parser.add_argument("--o", dest="ofile", required=False, help="OUTPUT_SNP2GENE_REPORTED")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  fout = open(args.ofile,"w") if args.ofile else sys.stdout

  assn = pd.read_csv(args.INPUT_GWASCAT_ASSN_FILE, "\t", low_memory=False)
  logging.info(f"{args.INPUT_GWASCAT_ASSN_FILE} rows: {assn.shape[0]}")
  logging.debug(f"{args.INPUT_GWASCAT_ASSN_FILE} columns: {str(assn.columns)}")
  #assn.info()

  logging.debug(f"Check HbA1c measurement (EFO_0004541): {assn[assn['MAPPED_TRAIT_URI'].str.contains('EFO_0004541', na=False)].shape[0]} rows")
  logging.debug(f"Check HbA1c measurement (EFO_0004541) GCSTs: {assn[assn['MAPPED_TRAIT_URI'].str.contains('EFO_0004541', na=False)]['STUDY_ACCESSION'].drop_duplicates().to_list()}")
  gcsts_test = ("GCST002390", "GCST005145", "GCST006001")
  logging.debug(f"Check HbA1c measurement (previous) GCSTs ({','.join(gcsts_test)}) traits: {assn[assn['STUDY_ACCESSION'].isin(gcsts_test)].shape[0]} rows")
  logging.debug(f"Check HbA1c measurement (previous) GCSTs ({','.join(gcsts_test)}) traits: {assn[assn['STUDY_ACCESSION'].isin(gcsts_test)]['MAPPED_TRAIT'].drop_duplicates().to_list()}")

  # REPORTED
  s2gr = assn[["STUDY_ACCESSION", "SNPS", "REPORTED_GENE(S)"]].drop_duplicates()
  s2gr.columns = ["STUDY_ACCESSION", "SNPS", "GSYMB"]
  logging.debug(f"Initial SNP2GENE (REPORTED) rows: {s2gr.shape[0]}")
  #s2gr.info()
  
  # Split delimited SNPS to multiple rows:
  col_split="SNPS"
  s2gr = s2gr.astype({col_split:str})
  s2gr = s2gr.assign(**{col_split:s2gr[col_split].str.split('[,x ]+')})
  s2gr = pd.DataFrame({
	col:np.repeat(s2gr[col].values, s2gr[col_split].str.len())
	for col in s2gr.columns.difference([col_split])
	}).assign(**{col_split:np.concatenate(s2gr[col_split].values)})[s2gr.columns.tolist()]
  s2gr = s2gr.drop_duplicates()
  logging.info(f"SNP2GENE (REPORTED) after spliting delimited SNPs, rows: {s2gr.shape[0]}")

  # Split delimited GENES to multiple rows:
  col_split="GSYMB"
  s2gr = s2gr.astype({col_split:str})
  s2gr = s2gr.assign(**{col_split:s2gr[col_split].str.split(r'[,/;:| ]+')})
  s2gr = pd.DataFrame({
	col:np.repeat(s2gr[col].values, s2gr[col_split].str.len())
	for col in s2gr.columns.difference([col_split])
	}).assign(**{col_split:np.concatenate(s2gr[col_split].values)})[s2gr.columns.tolist()]
  s2gr = s2gr.drop_duplicates()
  logging.info(f"After spliting delimited GSYMB, rows: {s2gr.shape[0]}")

  for badname in ("intergenic", "NR", "x", "nan"):
    s2gr = s2gr[(s2gr["GSYMB"] != badname)]
    logging.info(f"After removing GSYMB \"{badname}\", rows: {s2gr.shape[0]}")

  # Remove non-alphanumeric and len==1.
  good = s2gr["GSYMB"].str.contains(r'^[A-Za-z][A-Za-z0-9\-\.]+$')
  badset = set(s2gr["GSYMB"][~good].to_list())
  logging.debug(f"Removing GSYMB badset \"{str(badset)}\"")
  s2gr = s2gr[good]
  logging.info(f"Removed GSYMB badset, rows: {s2gr.shape[0]}")

  s2gr = s2gr.sort_values(["GSYMB"])
  s2gr["ENSG"] = ""
  s2gr["MAPPED_OR_REPORTED"] = "r"
  s2gr.columns = ["STUDY_ACCESSION", "SNPS", "GSYMB", "ENSG", "MAPPED_OR_REPORTED"]
  logging.info(f"SNP2GENE (REPORTED) SNPs: {s2gr.SNPS.nunique()}")

  # MAPPED
  s2gm = assn[["STUDY_ACCESSION",
	"SNPS", # (delimited?)
	"MAPPED_GENE", # (symbol[s])
	"UPSTREAM_GENE_ID", # (ENSG)
	"DOWNSTREAM_GENE_ID", # (ENSG)
	"SNP_GENE_IDS" # (comma delimited ENSGs)
	]].drop_duplicates()
  logging.debug(f"Initial SNP2GENE (MAPPED) rows: {s2gm.shape[0]}")

  # Split delimited ENSGs to multiple rows:
  col_split="SNP_GENE_IDS"
  s2gm = s2gm.astype({col_split:str})
  s2gm = s2gm.assign(**{col_split:s2gm[col_split].str.split(r'[, ]+')})
  s2gm = pd.DataFrame({
	col:np.repeat(s2gm[col].values, s2gm[col_split].str.len())
	for col in s2gm.columns.difference([col_split])
	}).assign(**{col_split:np.concatenate(s2gm[col_split].values)})[s2gm.columns.tolist()]
  s2gm = s2gm.drop_duplicates()
  logging.info(f"After spliting delimited ENSGs, rows: {s2gm.shape[0]}")
  s2gm["ENSG"] = s2gm["SNP_GENE_IDS"]

  s2gm["GSYMB"] = ""
  s2gm["MAPPED_OR_REPORTED"] = ""
  s2gm = s2gm.replace(np.nan, "", regex=True)
  s2gm = s2gm.replace("nan", "", regex=True)
  #s2gm[["SNP_GENE_IDS","UPSTREAM_GENE_ID","DOWNSTREAM_GENE_ID"]].to_csv("s2gm_DEBUG.tsv", "\t", index=False)
  n_m = (s2gm["SNP_GENE_IDS"]!="").sum()
  n_mu = (s2gm["UPSTREAM_GENE_ID"]!="").sum()
  n_md = (s2gm["DOWNSTREAM_GENE_ID"]!="").sum()
  logging.info(f"n_m={n_m}; n_mu={n_mu}; n_md={n_md}; n_total={n_m+n_mu+n_md}; n_row={s2gm.shape[0]}")
  #
  s2gm_in = s2gm.loc[s2gm["SNP_GENE_IDS"]!=""]
  s2gm_in.loc[:, "GSYMB"] = s2gm_in.loc[:, "MAPPED_GENE"]
  s2gm_in.loc[:, "MAPPED_OR_REPORTED"] = "m"
  #
  # Intergenic SNPs have BOTH up- and downstream. Such rows transformed into two rows.
  s2gm_up = s2gm.loc[s2gm["UPSTREAM_GENE_ID"]!=""]
  s2gm_up.loc[:, "GSYMB"] = s2gm_up.loc[:, "MAPPED_GENE"]
  s2gm_up.loc[:, "ENSG"] = s2gm_up.loc[:, "UPSTREAM_GENE_ID"]
  s2gm_up.loc[:, "MAPPED_OR_REPORTED"] = "mu"
  #
  s2gm_down = s2gm.loc[s2gm["DOWNSTREAM_GENE_ID"]!=""]
  s2gm_down.loc[:, "GSYMB"] = s2gm_down.loc[:, "MAPPED_GENE"]
  s2gm_down.loc[:, "ENSG"] = s2gm_down.loc[:, "DOWNSTREAM_GENE_ID"]
  s2gm_down.loc[:, "MAPPED_OR_REPORTED"] = "md"
  #
  s2gm = pd.concat([s2gm_in, s2gm_up, s2gm_down])
  #
  tag="MAPPED_OR_REPORTED"
  for key,val in s2gm[tag].value_counts().iteritems(): logging.info(f"Count({tag}==\"{key}\"): {val:6d}")
  #
  s2gm = s2gm[["STUDY_ACCESSION", "SNPS", "GSYMB", "ENSG", "MAPPED_OR_REPORTED"]]
  logging.info(f"SNP2GENE (MAPPED) SNPs: {s2gm.SNPS.nunique()}")
  #
  s2g = pd.concat([s2gr, s2gm])
  s2g = s2g.sort_values(["STUDY_ACCESSION", "SNPS", "GSYMB", "ENSG", "MAPPED_OR_REPORTED"])
  logging.info(f"SNP2GENE (ALL) SNPs: {s2g.SNPS.nunique()}")
  #
  for key,val in s2g[tag].value_counts().iteritems(): logging.info(f"Count({tag}==\"{key}\"): {val:6d}")
  #
  s2g.to_csv(fout, "\t", index=False)
  #

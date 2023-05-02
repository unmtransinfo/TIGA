#!/usr/bin/env python3
###
# Replacement for Go_TIGA_DbCreate.sh which for the workflow
# is needed only to produce gwascat_counts.tsv.
# Also produces trait_counts.tsv, not needed.
###
import sys,os,re,logging,argparse
import pandas as pd
import numpy as np

###
# Input files:
#	gwascat_gwas.tsv	#gwascat_gwas.R
#	gwascat_assn.tsv	#gwascat_assn.R
#	gwascat_snp2gene.tsv	#snp2gene.py
#	gwascat_trait.tsv	#gwascat_trait.R
#	gwascat_icite.tsv	#BioClients.icite.Client (currently not used by this program)
# Output files:
#	gwascat_counts.tsv
#	trait_counts.tsv
###

#############################################################################
def ComputeCounts(gwas, assn, trait, snp2gene, icite, fout_gwas, fout_trait):
  """
gwas_counts table (DataFrame):
	study_accession - study ID
	snp_count - per study
	trait_count - per study
	assn_count - per study
	gene_r_count - reported genes per study
	gene_m_count - mapped genes per study
	pubNstudy - studies per PMID (PMID for this study)
  """
  #
  gwas_counts = gwas[["STUDY_ACCESSION"]]

  snp_counts = snp2gene[["STUDY_ACCESSION", "SNP"]].groupby("STUDY_ACCESSION").agg(snp_count=pd.NamedAgg(column="SNP", aggfunc="nunique"))
  gwas_counts = pd.merge(gwas_counts, snp_counts, how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")

  trait_counts = trait[["STUDY_ACCESSION", "MAPPED_TRAIT_URI"]].groupby("STUDY_ACCESSION").agg(trait_count=pd.NamedAgg(column="MAPPED_TRAIT_URI", aggfunc="nunique"))
  gwas_counts = pd.merge(gwas_counts, trait_counts, how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")

  assn_counts = assn[["STUDY_ACCESSION", "P-VALUE"]].groupby("STUDY_ACCESSION").agg(assn_count=pd.NamedAgg(column="P-VALUE", aggfunc="nunique"))
  gwas_counts = pd.merge(gwas_counts, assn_counts, how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")

  gene_r_counts = snp2gene[snp2gene["MAPPED_OR_REPORTED"]=="r"][["STUDY_ACCESSION", "ENSG"]].groupby("STUDY_ACCESSION").agg(gene_r_count=pd.NamedAgg(column="ENSG", aggfunc="nunique"))
  gwas_counts = pd.merge(gwas_counts, gene_r_counts, how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")

  gene_m_counts = snp2gene[snp2gene["MAPPED_OR_REPORTED"].isin(["m", "mu", "md"])][["STUDY_ACCESSION", "ENSG"]].groupby("STUDY_ACCESSION").agg(gene_m_count=pd.NamedAgg(column="ENSG", aggfunc="nunique"))
  gwas_counts = pd.merge(gwas_counts, gene_m_counts, how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")

  perpub_counts = gwas[["STUDY_ACCESSION", "PUBMEDID"]].groupby("PUBMEDID").agg(pubNstudy=pd.NamedAgg(column="STUDY_ACCESSION", aggfunc="nunique"))
  pub_counts = pd.merge(gwas[["STUDY_ACCESSION", "PUBMEDID"]], perpub_counts, how="left", on="PUBMEDID")
  gwas_counts = pd.merge(gwas_counts, pub_counts[["STUDY_ACCESSION", "pubNstudy"]], how="left", on="STUDY_ACCESSION")
  logging.debug(f"GWAS_COUNTS columns: {str(gwas_counts.columns)}")
  logging.debug(f"GWAS_COUNTS rows: {gwas_counts.shape[0]}; STUDY_ACCESSIONs: {gwas_counts.STUDY_ACCESSION.nunique()}")
  #
  gwas_counts = gwas_counts.fillna(0)
  gwas_counts.info()
  gwas_counts.describe().round(3).to_csv(sys.stderr, sep="\t")
  gwas_counts = gwas_counts.rename(columns={"STUDY_ACCESSION":"study_accession"})
  gwas_counts.to_csv(fout_gwas, sep="\t", index=False)

  # Counts per trait. Not used at this time.
  trait_counts = trait[["STUDY_ACCESSION", "MAPPED_TRAIT_URI", "MAPPED_TRAIT"]].groupby("MAPPED_TRAIT_URI").agg(
	mapped_trait_uri=pd.NamedAgg(column="MAPPED_TRAIT_URI", aggfunc="first"),
	mapped_trait=pd.NamedAgg(column="MAPPED_TRAIT", aggfunc="first"),
	n_study=pd.NamedAgg(column="STUDY_ACCESSION", aggfunc="nunique"))
  trait_counts = trait_counts.fillna(0)
  trait_counts.info()
  trait_counts.describe().round(3).to_csv(sys.stderr, sep="\t")
  trait_counts.to_csv(fout_trait, sep="\t", index=False)

#############################################################################
if __name__=="__main__":
  parser = argparse.ArgumentParser(description='TIGA GWAS Catalog counts')
  parser.add_argument("--ifile_gwas", required=True, help="Cleaned studies, from gwascat_gwas.R")
  parser.add_argument("--ifile_assn", required=True, help="Cleaned associations, from gwascat_assn.R")
  parser.add_argument("--ifile_trait", required=True, help="Cleaned traits, from gwascat_trait.R")
  parser.add_argument("--ifile_snp2gene", required=True, help="Mappings, from snp2gene.py")
  parser.add_argument("--ifile_icite", help="Publication iCite RCRs, from BioClients.icite.Client (currently not used by this program)")
  parser.add_argument("--ofile_gwas", help="Counts per study")
  parser.add_argument("--ofile_trait", help="Counts per trait")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  fout_gwas = open(args.ofile_gwas,"w") if args.ofile_gwas else sys.stdout
  fout_trait = open(args.ofile_trait,"w") if args.ofile_trait else sys.stdout

  gwas = pd.read_csv(args.ifile_gwas, sep="\t", low_memory=False)
  logging.debug(f"{args.ifile_gwas} columns: {str(gwas.columns)}")
  logging.info(f"{args.ifile_gwas} rows: {gwas.shape[0]}; studies: {gwas.STUDY_ACCESSION.nunique()}")

  assn = pd.read_csv(args.ifile_assn, sep="\t", low_memory=False)
  logging.debug(f"{args.ifile_assn} columns: {str(assn.columns)}")
  logging.info(f"{args.ifile_assn} rows: {assn.shape[0]}")

  trait = pd.read_csv(args.ifile_trait, sep="\t", low_memory=False)
  logging.debug(f"{args.ifile_trait} columns: {str(trait.columns)}")
  logging.info(f"{args.ifile_trait} rows: {trait.shape[0]}; traits: {trait.MAPPED_TRAIT_URI.nunique()}")

  snp2gene = pd.read_csv(args.ifile_snp2gene, sep="\t", low_memory=False)
  logging.debug(f"{args.ifile_snp2gene} columns: {str(snp2gene.columns)}")
  logging.info(f"{args.ifile_snp2gene} rows: {snp2gene.shape[0]}; genes: {snp2gene.ENSG.nunique()}; snps: {snp2gene.SNP.nunique()}")

  icite = pd.read_csv(args.ifile_icite, sep="\t", low_memory=False) if args.ifile_icite else None
  if icite is not None:
    logging.debug(f"{args.ifile_icite} columns: {str(icite.columns)}")
    logging.info(f"{args.ifile_icite} rows: {icite.shape[0]}; PMIDs: {icite.pmid.nunique()}")

  ComputeCounts(gwas, assn, trait, snp2gene, icite, fout_gwas, fout_trait)

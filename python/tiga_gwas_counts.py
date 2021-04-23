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
#	gwascat_gwas.tsv #gwascat_gwas.R
#	gwascat_assn.tsv #gwascat_assn.R
#	gwascat_snp2gene.tsv #snp2gene.py
#	gwascat_trait.tsv #gwascat_trait.R
#	gwascat_icite.tsv #BioClients.icite.Client
# Output files:
#	gwascat_counts.tsv
#	trait_counts.tsv
###

if __name__=="__main__":
  parser = argparse.ArgumentParser(description='TIGA GWAS Catalog counts')
  parser.add_argument("--ifile_gwas", required=True, help="Cleaned studies, from gwascat_gwas.R")
  parser.add_argument("--ifile_assn", required=True, help="Cleaned associations, from gwascat_assn.R")
  parser.add_argument("--ifile_trait", required=True, help="Cleaned traits, from gwascat_trait.R")
  parser.add_argument("--ifile_snp2gene", required=True, help="Mappings, from snp2gene.py")
  parser.add_argument("--ifile_icite", required=True, help="Publication iCite RCRs, from BioClients.icite.Client")
  parser.add_argument("--o_gwas", dest="ofile_gwas", help="Counts per study")
  parser.add_argument("--o_trait", dest="ofile_trait", help="Counts per trait")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  fout_gwas = open(args.ofile_gwas,"w") if args.ofile_gwas else sys.stdout
  fout_trait = open(args.ofile_trait,"w") if args.ofile_trait else sys.stdout

  gwas = pd.read_csv(args.ifile_gwas, "\t", low_memory=False)
  logging.debug(f"{args.ifile_gwas} columns: {str(gwas.columns)}")
  logging.info(f"{args.ifile_gwas} rows: {gwas.shape[0]}; studies: {gwas.STUDY_ACCESSION.nunique()}")
  #gwas.info()

  assn = pd.read_csv(args.ifile_assn, "\t", low_memory=False)
  logging.debug(f"{args.ifile_assn} columns: {str(assn.columns)}")
  logging.info(f"{args.ifile_assn} rows: {assn.shape[0]}")

  trait = pd.read_csv(args.ifile_trait, "\t", low_memory=False)
  logging.debug(f"{args.ifile_trait} columns: {str(trait.columns)}")
  logging.info(f"{args.ifile_trait} rows: {trait.shape[0]}; traits: {trait.MAPPED_TRAIT_URI.nunique()}")

  snp2gene = pd.read_csv(args.ifile_snp2gene, "\t", low_memory=False)
  logging.debug(f"{args.ifile_snp2gene} columns: {str(snp2gene.columns)}")
  logging.info(f"{args.ifile_snp2gene} rows: {snp2gene.shape[0]}; genes: {snp2gene.ENSG.nunique()}; snps: {snp2gene.SNP.nunique()}")

  icite = pd.read_csv(args.ifile_icite, "\t", low_memory=False)
  logging.debug(f"{args.ifile_icite} columns: {str(icite.columns)}")
  logging.info(f"{args.ifile_icite} rows: {icite.shape[0]}; PMIDs: {icite.pmid.nunique()}")

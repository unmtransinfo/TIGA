#!/usr/bin/env python3
"""
Check for given trait and/or gene, from Catalog source files;
report studies, publications, SNPs, p-values, 
reported and mapped genes, and genes mapped by Ensembl.
Should this gene|trait|gene-trait be in TIGA dataset?
Inputs: 
 - Catalog studies file
 - Catalog association file (annotated)
 - Ensembl human genes file 
"""
###
import sys,os,re,argparse,time,logging,tqdm
import pandas as pd
import numpy as np

ASSN_COLS = ["MAPPED_TRAIT_URI", "SNPS", "MAPPED_GENE", "SNP_GENE_IDS", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "P-VALUE", "OR or BETA", "STUDY ACCESSION", "PUBMEDID"]

#############################################################################
def CheckGeneSnps(ensemblId, gwas, assn, fout):
  assn_gene_this = assn.loc[(assn.SNP_GENE_IDS.str.find(ensemblId)>=0)]
  assn_gene_this = pd.concat([assn_gene_this, assn.loc[(assn.UPSTREAM_GENE_ID.str.find(ensemblId)>=0)]])
  assn_gene_this = pd.concat([assn_gene_this, assn.loc[(assn.DOWNSTREAM_GENE_ID.str.find(ensemblId)>=0)]])
  if assn_gene_this.empty:
    logging.info(f"ASSN: Gene {ensemblId} not found in associations.")
    return
  logging.info(f"ASSN: Gene {ensemblId} found; associations: {assn_gene_this.shape[0]}")
  assn_gene_this = assn_gene_this.sort_values(["SNPS"])
  assn_gene_this[ASSN_COLS].drop_duplicates().to_csv(fout, "\t", index=False)

#############################################################################
def CheckTraitSnps(efoId, gwas, assn, fout):
  assn_trait_this =  assn.loc[(assn.MAPPED_TRAIT_URI.str.find(efoId)>=0)]
  if assn_trait_this.empty:
    logging.info(f"ASSN: Trait {efoId} not found in associations.")
    return
  logging.info(f"ASSN: Trait {efoId} found; studies: {assn_trait_this['STUDY ACCESSION'].nunique()}; associations: {assn_trait_this.shape[0]}")
  assn_trait_this = assn_trait_this.sort_values(["SNPS"])
  assn_trait_this[ASSN_COLS].drop_duplicates().to_csv(fout, "\t", index=False)

#############################################################################
def CheckGeneTraitSnps(ensemblId, efoId, gwas, assn, fout):
  CheckGeneSnps(ensemblId, gwas, assn, fout)
  CheckTraitSnps(efoId, gwas, assn, fout)
  assn_gt_this =  assn.loc[((assn.MAPPED_TRAIT_URI.str.find(efoId)>=0) & (assn.SNP_GENE_IDS.str.find(ensemblId)>=0))]
  if assn_gt_this.empty:
    logging.info(f"ASSN: Gene-Trait {ensemblId}-{efoId} not found in associations.")
    return
  logging.info(f"ASSN: Gene-Trait {ensemblId}-{efoId} found; associations: {assn_gt_this.shape[0]}")
  assn_gt_this = assn_gt_this.sort_values(["SNPS"])
  assn_gene_this[ASSN_COLS].drop_duplicates().to_csv(fout, "\t", index=False)

#############################################################################
def CheckStudySnps(gwasId, gwas, assn, fout):
  assn_this =  assn.loc[(assn["STUDY ACCESSION"]==gwasId)]
  if assn_this.empty:
    logging.info(f"ASSN: Study {gwasId} not found in associations.")
    return
  assn_this = assn_this.sort_values(["SNPS"])
  ensgs = set()
  for i in range(assn_this.shape[0]):
    logging.debug(f"{i+1}. SNP_GENE_IDS: \"{assn_this.SNP_GENE_IDS.astype('str').values[i]}\"; DOWNSTREAM_GENE_ID: \"{assn_this.DOWNSTREAM_GENE_ID.astype('str').values[i]}\"; UPSTREAM_GENE_ID: \"{assn_this.UPSTREAM_GENE_ID.astype('str').values[i]}\"")
    ensgs_this = set([assn_this.SNP_GENE_IDS.astype('str').values[i]]+re.split(r"[\s,]+", assn_this.UPSTREAM_GENE_ID.astype('str').values[i])+re.split(r"[\s,]+", assn_this.DOWNSTREAM_GENE_ID.astype('str').values[i]))
    ensgs_this -= set(["nan"])
    ensgs |= ensgs_this
    logging.info(f"{i+1}. SNP: {assn_this.SNPS.values[i]}: ensemblIds ({len(ensgs_this)}): {str(ensgs_this)}")
  logging.info(f"ASSN: Study {gwasId} found; associations: {assn_this.shape[0]}; SNPs: {assn_this.SNPS.nunique()}; mapped genes: {len(ensgs)}")
  assn_this[ASSN_COLS].drop_duplicates().to_csv(fout, "\t", index=False)

#############################################################################
if __name__=="__main__":
  SRCDATADIR=os.environ["HOME"]+"/../data/GWASCatalog/releases/2020/12/16/"
  GWAS_FILE=SRCDATADIR+"/gwas-catalog-studies_ontology-annotated.tsv"
  ASSN_FILE=SRCDATADIR+"/gwas-catalog-associations_ontology-annotated.tsv"
  epilog="""\
Example traits: EFO_0004541
Example genes: ENSG00000160785
Example studies: GCST006001, GCST005145, GCST002390
"""
  OPS=["checkStudy", "checkGene", "checkGeneTrait"]
  parser = argparse.ArgumentParser(description="Check GENE and/or TRAIT, for SNPs etc.", epilog=epilog)
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("-g", "--ensemblId", help="Input gene (Ensembl ID)")
  parser.add_argument("-t", "--efoId", help="Input trait (EFO ID)")
  parser.add_argument("-s", "--studyId", help="Input study (GST ID)")
  parser.add_argument("--gwas_file", default=GWAS_FILE, help="GWAS Catalog study file")
  parser.add_argument("--assn_file", default=ASSN_FILE, help="GWAS Catalog association file")
  parser.add_argument("--o", dest="ofile", help="(TSV)")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG if args.verbose>0 else logging.INFO)

  logging.debug(f"gwas_file: {args.gwas_file}")
  logging.debug(f"assn_file: {args.assn_file}")

  fout = open(args.ofile, "w") if args.ofile else sys.stdout
  gwas = pd.read_csv(args.gwas_file, sep="\t", dtype=str, low_memory=False)
  #gwas.info()
  assn = pd.read_csv(args.assn_file, sep="\t", dtype=str, low_memory=False)
  #assn.info()
  logging.info(f"GWAS: TRAIT URIs: {gwas.MAPPED_TRAIT_URI.nunique()}")
  logging.info(f"ASSN: TRAIT URIs: {assn.MAPPED_TRAIT_URI.nunique()}")

  if args.op=="checkStudy":
    CheckStudySnps(args.studyId, gwas, assn, fout)

  else:
    if not (args.efoId or args.ensemblId):
      parser.error("Either -g or -t required.")
    CheckGeneTraitSnps(args.ensemblId, args.efoId, gwas, assn, fout)


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

#############################################################################
def CheckGeneTraitSnps(ensemblId, efoId, gwas_file, assn_file, ofile):
  fout = open(ofile, "w") if ofile else sys.stdout

  gwas = pd.read_csv(gwas_file, sep="\t")
  assn = pd.read_csv(assn_file, sep="\t")
  assn_cols = ["MAPPED_TRAIT_URI", "MAPPED_GENE", "SNP_GENE_IDS", "SNPS", "P-VALUE", "P-VALUE (TEXT)", "OR or BETA", "STUDY ACCESSION", "PUBMEDID"]

  logging.info(f"GWAS: TRAIT URIs: {gwas.MAPPED_TRAIT_URI.nunique()}")
  logging.info(f"ASSN: TRAIT URIs: {assn.MAPPED_TRAIT_URI.nunique()}")
  if efoId:
    assn_trait_this =  assn.loc[(assn.MAPPED_TRAIT_URI.str.find(efoId)>=0)]
    if assn_trait_this.shape[0]>0:
      logging.info(f"ASSN: Trait {efoId} found; associations: {assn_trait_this.shape[0]}")
      assn_trait_this[assn_cols].drop_duplicates().to_csv(fout, "\t", index=False)
    else:
      logging.info(f"ASSN: Trait {efoId} not found in associations.")
  if ensemblId:
    assn_gene_this = assn.loc[(assn.SNP_GENE_IDS.str.find(ensemblId)>=0)]
    assn_gene_this = pd.concat([assn_gene_this, assn.loc[(assn.UPSTREAM_GENE_ID.str.find(ensemblId)>=0)]])
    assn_gene_this = pd.concat([assn_gene_this, assn.loc[(assn.DOWNSTREAM_GENE_ID.str.find(ensemblId)>=0)]])
    if assn_gene_this.shape[0]>0:
      logging.info(f"ASSN: Gene {ensemblId} found; associations: {assn_gene_this.shape[0]}")
      assn_gene_this[assn_cols].drop_duplicates().to_csv(fout, "\t", index=False)
    else:
      logging.info(f"ASSN: Gene {ensemblId} not found in associations.")

  if efoId and ensemblId:
    assn_gt_this =  assn.loc[((assn.MAPPED_TRAIT_URI.str.find(efoId)>=0) & (assn.SNP_GENE_IDS.str.find(ensemblId)>=0))]
    if assn_gt_this.shape[0]>0:
      logging.info(f"ASSN: Gene-Trait {ensemblId}-{efoId} found; associations: {assn_gt_this.shape[0]}")
      assn_gene_this[assn_cols].drop_duplicates().to_csv(fout, "\t", index=False)
    else:
      logging.info(f"ASSN: Gene-Trait {ensemblId}-{efoId} not found in associations.")

#############################################################################
if __name__=="__main__":
  SRCDATADIR=os.environ["HOME"]+"/../data/GWASCatalog/releases/2020/12/16/"
  GWAS_FILE=SRCDATADIR+"/gwas-catalog-studies_ontology-annotated.tsv"
  ASSN_FILE=SRCDATADIR+"/gwas-catalog-associations_ontology-annotated.tsv"
  epilog="""\
Example traits: EFO_0004541
Example genes: ENSG00000160785
"""
  parser = argparse.ArgumentParser(description="Check GENE and/or TRAIT, for SNPs etc.", epilog=epilog)
  parser.add_argument("-g", "--ensemblId", help="Input gene (Ensembl ID)")
  parser.add_argument("-t", "--efoId", help="Input trait (EFO ID)")
  parser.add_argument("--gwas_file", default=GWAS_FILE, help="GWAS Catalog study file")
  parser.add_argument("--assn_file", default=ASSN_FILE, help="GWAS Catalog association file")
  parser.add_argument("--o", dest="ofile", help="(TSV)")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG if args.verbose>0 else logging.INFO)

  if not (args.efoId or args.ensemblId):
    parser.error("Either -g or -t required.")
    sys.exit()

  logging.debug(f"gwas_file: {args.gwas_file}")
  logging.debug(f"assn_file: {args.assn_file}")

  CheckGeneTraitSnps(args.ensemblId, args.efoId, args.gwas_file, args.assn_file, args.ofile)


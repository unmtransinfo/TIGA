#!/bin/bash
#
cwd=$(pwd)
DATADIR=${cwd}/data
#
${cwd}/python/tiga_gt_snpcheck.py checkGeneTrait \
	-g ENSG00000160785 -t EFO_0004541 --o $DATADIR/tiga_gt_snpcheck_GeneTrait_ENSG00000160785-EFO_0004541.tsv
#
${cwd}/python/tiga_gt_snpcheck.py checkStudy \
	-s GCST002390 --o $DATADIR/tiga_gt_snpcheck_Study_GCST002390.tsv 
${cwd}/python/tiga_gt_snpcheck.py checkStudy \
	-s GCST005145 --o $DATADIR/tiga_gt_snpcheck_Study_GCST005145.tsv 
${cwd}/python/tiga_gt_snpcheck.py checkStudy \
	-s GCST006001 --o $DATADIR/tiga_gt_snpcheck_Study_GCST006001.tsv 
#
cat $DATADIR/tiga_gt_snpcheck_Study_GCST002390.tsv |sed '1d' |awk -F '\t' '{print $2}' >$DATADIR/tiga_gt_snpcheck_Study_GCST002390.rs
cat $DATADIR/tiga_gt_snpcheck_Study_GCST005145.tsv |sed '1d' |awk -F '\t' '{print $2}' >$DATADIR/tiga_gt_snpcheck_Study_GCST005145.rs
cat $DATADIR/tiga_gt_snpcheck_Study_GCST006001.tsv |sed '1d' |awk -F '\t' '{print $2}' >$DATADIR/tiga_gt_snpcheck_Study_GCST006001.rs
cat $DATADIR/tiga_gt_snpcheck_Study_GCST00*.rs >$DATADIR/z.rs
python3 -m BioClients.util.pandas.Utils --i data/old/gwascat_Snps-20200813.tsv.gz --coltags rsId,gene_ensemblGeneIds selectcols |grep ENSG00000160785 |sort -u > data/old/gwascat_Snps-20200813_ENSG00000160785.tsv
cat $DATADIR/old/gwascat_Snps-20200813_ENSG00000160785.tsv |awk -F '\t' '{print $1}' >$DATADIR/old/gwascat_Snps-20200813_ENSG00000160785.rs
${cwd}/python/setman.py --iA $DATADIR/z.rs --iB $DATADIR/old/gwascat_Snps-20200813_ENSG00000160785.rs AandB
#
# rs2273833
# rs6684514
# rs144991356
# INFO:Unique entities in output: 3
###
# Notes:
# rs2273833 (GCST006001; pValue=1E-18)
# rs6684514 (GCST002390; pValue=1E-23), (GCST005145; pValue=1E-7)
# rs144991356 (GCST005145; pValue=3E-8)
# 
# gwascat_Snps.tsv from GWAS Catalog REST API via
# BioClients.gwascatalog.Client get_snps, e.g.
# https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs2273833
# https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs6684514
# https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs144991356
#
python3 -m BioClients.gwascatalog.Client get_snps --ids rs2273833,rs6684514,rs144991356 -v -v --o $DATADIR/gwascatalog_Snps_rs2273833-rs6684514-rs144991356.tsv
# ALL THREE SNPs MAPPED TO ENSG00000160785 (SLC25A44), accessed 2021-03-30.

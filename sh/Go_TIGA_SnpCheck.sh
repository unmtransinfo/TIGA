#!/bin/bash
#
cwd=$(pwd)
#
#SRCDATADIR="/home/data/GWASCatalog/releases/2021/02/12"
#DATADIR="${cwd}/data/20210212_GCAPI"
#
SRCDATADIR="/home/data/GWASCatalog/releases/2021/03/29"
DATADIR="${cwd}/data/20210329"
#
gwas_file="${SRCDATADIR}/gwas-catalog-studies_ontology-annotated.tsv"
assn_file="${SRCDATADIR}/gwas-catalog-associations_ontology-annotated.tsv"
#
ensgId="ENSG00000160785" #SLC25A44 (Solute carrier family 25 member 44)
efoId="EFO_0004541" #HbA1c measurement
#
${cwd}/python/tiga_gt_snpcheck.py checkGeneTrait \
	-g ${ensgId} -t ${efoId} --o $DATADIR/tiga_gt_snpcheck_GeneTrait_${ensgId}-${efoId}.tsv
#
gcstIds="GCST002390 GCST005145 GCST006001"
rm -f $DATADIR/z.rs
touch $DATADIR/z.rs
for gcstId in $gcstIds ; do
	${cwd}/python/tiga_gt_snpcheck.py checkStudy \
		--gwas_file $gwas_file --assn_file $assn_file \
		-s ${gcstId} --o $DATADIR/tiga_gt_snpcheck_Study_${gcstId}.tsv 
	cat $DATADIR/tiga_gt_snpcheck_Study_${gcstId}.tsv |sed '1d' |awk -F '\t' '{print $2}' >>$DATADIR/z.rs
done
rsfile="$DATADIR/tiga_gt_snpcheck_Study_$(echo $gcstIds |sed 's/ /_/g').rs"
cat $DATADIR/z.rs |sort -u >$rsfile
printf "SNP count for studies ($(echo ${gcstIds} |sed 's/ /,/g')): $(cat $rsfile |wc -l)\n"
###
#python3 -m BioClients.util.pandas.Utils --i data/old/gwascat_Snps-20200813.tsv.gz --coltags rsId,gene_ensemblGeneIds selectcols |grep ${ensgId} |sort -u > data/old/gwascat_Snps-20200813_${ensgId}.tsv
#cat $DATADIR/old/gwascat_Snps-20200813_${ensgId}.tsv |awk -F '\t' '{print $1}' >$DATADIR/old/gwascat_Snps-20200813_${ensgId}.rs
#${cwd}/python/setman.py --iA $DATADIR/z.rs --iB $DATADIR/old/gwascat_Snps-20200813_${ensgId}.rs AandB
###
python3 -m BioClients.util.pandas.Utils --i $DATADIR/gwascat_snp2gene_API.tsv --coltags rsId,ensemblGeneIds selectcols |grep ${ensgId} |sort -u > $DATADIR/gwascat_snp2gene_API_${ensgId}.tsv
cat $DATADIR/gwascat_snp2gene_API_${ensgId}.tsv |awk -F '\t' '{print $1}' |sort -u >$DATADIR/gwascat_snp2gene_API_${ensgId}.rs
${cwd}/python/setman.py AandB --iA $DATADIR/z.rs --iB $DATADIR/gwascat_snp2gene_API_${ensgId}.rs
###
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
# ALL THREE SNPs MAPPED TO ENSG00000160785 (SLC25A44), accessed 2021-03-30.
python3 -m BioClients.gwascatalog.Client get_snps -q --ids rs2273833,rs6684514,rs144991356 --o $DATADIR/gwascatalog_Snps_rs2273833-rs6684514-rs144991356.tsv
cat $DATADIR/gwascatalog_Snps_rs2273833-rs6684514-rs144991356.tsv |grep ${ensgId}
printf "Associations to ${ensgId}: $(cat $DATADIR/gwascatalog_Snps_rs2273833-rs6684514-rs144991356.tsv |grep ${ensgId}|wc -l)\n"
#

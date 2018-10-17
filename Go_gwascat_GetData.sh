#!/bin/sh
#############################################################################
### Go_gwascat_GetData.sh - Create CSV files for import to RDB.
### gt_stats.tsv is used by GWAX web app.
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
#############################################################################
### GWAS Catalog studies each have a "study_accession".
### Also are associated with a publication (PubMedID), but not uniquely.
### See https://www.ebi.ac.uk/gwas/docs/fileheaders
#############################################################################
### MAPPED GENE(S): Gene(s) mapped to the strongest SNP. If the SNP is located
### within a gene, that gene is listed. If the SNP is intergenic, the upstream
### and downstream genes are listed, separated by a hyphen.
#############################################################################
### "OR_or_BETA: Reported odds ratio or beta-coefficient associated with
### strongest SNP risk allele. Note that if an OR <1 is reported this is
### inverted, along with the reported allele, so that all ORs included in
### the Catalog are >1. Appropriate unit and increase/decrease are included
### for beta coefficients."
#############################################################################
### Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H,
### Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H. The NHGRI
### GWAS Catalog, a curated resource of SNP-trait associations. Nucleic
### Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006.
#############################################################################
### iCite annotations from iCite API, with all PMIDs from GWASCatalog.
#############################################################################
### Jeremy Yang
#############################################################################
#
cwd=$(pwd)
#
SRCDATADIR="/home/data/gwascatalog/data"
DATADIR="${cwd}/data"
#
#Source files:
#gwasfile="${SRCDATADIR}/gwas_catalog_v1.0.1-studies_r2017-10-10.tsv"
gwasfile="${SRCDATADIR}/gwas_catalog_v1.0.2-studies_r2018-09-30.tsv"
#
#assnfile="${SRCDATADIR}/gwas_catalog_v1.0.1-associations_e90_r2017-10-10.tsv"
assnfile="${SRCDATADIR}/gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv"
###
#Output files:
tsvfile_gwas="${DATADIR}/gwascat_gwas.tsv"
tsvfile_assn="${DATADIR}/gwascat_assn.tsv"
###
#Clean studies:
${cwd}/R/gwascat_gwas.R $gwasfile $tsvfile_gwas
#
###
#Clean, separate OR_or_beta into oddsratio, beta columns:
${cwd}/R/gwascat_assn.R $assnfile $tsvfile_assn
#
#############################################################################
#trait links:
traitfile="${DATADIR}/gwascat_trait.tsv"
#
#SNP to gene links:
snpgenefile="${DATADIR}/gwascat_snp2gene.tsv"
#
#############################################################################
### TRAITS:
####
printf "STUDY_ACCESSION\tMAPPED_TRAIT\tMAPPED_TRAIT_URI\n" >${traitfile}
#
cat $tsvfile_gwas \
	|sed -e '1d' \
	|perl -n perl/gwas2trait.pl \
	>>${traitfile}
#
#############################################################################
### REPORTED GENES:
#
printf "STUDY_ACCESSION\tGSYMB\tSNP\tREPORTED_OR_MAPPED\n" >${snpgenefile}
#
# "REPORTED_GENE(S),SNPS,STUDY_ACCESSION" (14, 22, 37)
###
cat $tsvfile_assn \
	|sed -e '1d' \
	|perl -n perl/assn2snpgene_reported.pl \
	>>${snpgenefile}
#
#############################################################################
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
# "MAPPED_GENE(S),SNPS,STUDY_ACCESSION" (15, 22, 37)
###
cat $tsvfile_assn \
	|sed -e '1d' \
	|perl -n perl/assn2snpgene_mapped.pl \
	>>${snpgenefile}
#
#############################################################################
### PMIDs:
cat $tsvfile_gwas \
	|sed -e '1d' \
	|awk -F '\t' '{print $2}' \
	|sort -nu >$DATADIR/gwascat.pmid
#
${cwd}/python/pubmed_icite.py get \
	--i $DATADIR/gwascat.pmid \
	--o $DATADIR/gwascat_icite.tsv
#
#############################################################################
### Entrez gene IDs: UPSTREAM_GENE_ID, DOWNSTREAM_GENE_ID, SNP_GENE_IDS
cat $tsvfile_assn \
	|sed -e '1d' \
	|awk -F '\t' '{print $16}' \
	|egrep -v '(^$|^NA$)' \
	|sort -nu \
	>$DATADIR/gwascat_upstream.geneid
cat $tsvfile_assn \
	|sed -e '1d' \
	|awk -F '\t' '{print $17}' \
	|egrep -v '(^$|^NA$)' \
	|sort -nu \
	>$DATADIR/gwascat_downstream.geneid
cat $tsvfile_assn \
	|sed -e '1d' \
	|awk -F '\t' '{print $18}' \
	|egrep -v '(^$|^NA$)' \
	|perl -ne 'print join("\n",split(/, */))' \
	|sort -nu \
	>$DATADIR/gwascat_snp.geneid
cat \
	$DATADIR/gwascat_upstream.geneid \
	$DATADIR/gwascat_downstream.geneid \
	$DATADIR/gwascat_snp.geneid \
	|sort -nu \
	>$DATADIR/gwascat.geneid
#
#############################################################################
# gt_stats.tsv generated by this R script:
#
${cwd}/R/gwascat_gt_stats.R \
	$DATADIR/gwascat_assn.tsv \
	$DATADIR/gwascat_snp2gene.tsv \
	$DATADIR/gwascat_trait.tsv \
	$DATADIR/tcrd_targets.csv \
	$DATADIR/gt_stats.tsv
#

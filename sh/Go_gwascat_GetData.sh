#!/bin/sh
#############################################################################
### Go_gwascat_GetData.sh - Create CSV files for import to RDB.
### gt_stats.tsv is used by GWAX web app.
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
#############################################################################
### According to EBI, additional information is available via the API,
### see https://www.ebi.ac.uk/gwas/docs/api,https://www.ebi.ac.uk/gwas/docs/faq#E,
### particularly on "Genomic Mappings" (SNP-gene) .
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
### Additional information via the API:
### But exactly what should be used?
${cwd}/python/sample_lines.py \
	--i $DATADIR/gwascat_gwas.tsv \
	--k 20 \
	chooseK \
	|awk -F '\t' '{print $15}' \
	|sed -e '1d' \
	>$DATADIR/z.gcst
${cwd}/python/gwascat_query.py \
	--i $DATADIR/z.gcst \
	--idtype "gcst" \
	--o $DATADIR/gwascat_query_studyAssociations_out.tsv \
	getStudyAssociations
#
${cwd}/python/sample_lines.py \
	--i $DATADIR/gwascat_snp2gene.tsv \
	--k 100 \
	chooseK \
	|awk -F '\t' '{print $3}' \
	|sed -e '1d' \
	>$DATADIR/z.rs
${cwd}/python/gwascat_query.py \
	--i $DATADIR/z.rs \
	--idtype "rs" \
	--o $DATADIR/gwascat_query_getSnps_out.tsv \
	getSnps
#
${cwd}/python/pandas_utils.py --i $DATADIR/gwascat_query_getSnps_out.tsv showcols
${cwd}/python/pandas_utils.py --i $DATADIR/gwascat_query_getSnps_out.tsv --coltag functionalClass colvalcounts
${cwd}/python/pandas_utils.py --i $DATADIR/gwascat_query_getSnps_out.tsv --coltag genomicContext_mappingMethod colvalcounts
${cwd}/python/pandas_utils.py --i $DATADIR/gwascat_query_getSnps_out.tsv --coltag genomicContext_isClosestGene colvalcounts
#
#############################################################################
# Gene-trait statistics:
# gt_stats.tsv generated by this R script:
#
${cwd}/R/gwascat_gt_stats.R \
	$DATADIR/gwascat_gwas.tsv \
	$DATADIR/gwascat_assn.tsv \
	$DATADIR/gwascat_snp2gene.tsv \
	$DATADIR/gwascat_trait.tsv \
	$DATADIR/gwascat_icite.tsv \
	$DATADIR/tcrd_targets.csv \
	$DATADIR/gt_stats.tsv
#

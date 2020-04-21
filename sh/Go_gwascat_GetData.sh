#!/bin/bash
#############################################################################
### Go_gwascat_GetData.sh - Create CSV files for import to RDB.
### gt_stats.tsv is used by TIGA web app.
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
#############################################################################
### According to EBI, additional information is available via the API,
### see https://www.ebi.ac.uk/gwas/docs/api,https://www.ebi.ac.uk/gwas/docs/faq#E,
### particularly on "Genomic Mappings" (SNP-gene) .
### Also it appears beta and OR values are separated better.
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
snp2genefile="${DATADIR}/gwascat_snp2gene.tsv"
#
#############################################################################
### TRAITS:
# Now generated by gwascat_trait.R
####
# Old way:
#printf "STUDY_ACCESSION\tMAPPED_TRAIT\tMAPPED_TRAIT_URI\n" >${traitfile}
#cat $tsvfile_gwas |sed -e '1d' |perl -n perl/gwas2trait.pl >>${traitfile}
####
# New way:
${cwd}/R/gwascat_trait.R $gwasfile $traitfile
#
#############################################################################
### REPORTED GENES:
#
printf "STUDY_ACCESSION\tGSYMB\tSNP\tREPORTED_OR_MAPPED\n" >${snp2genefile}
#
# "REPORTED_GENE(S),SNPS,STUDY_ACCESSION" (14, 22, 37)
###
cat $tsvfile_assn \
	|sed -e '1d' \
	|perl -n perl/snp2gene_reported.pl \
	>>${snp2genefile}
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
	|perl -n perl/snp2gene_mapped.pl \
	>>${snp2genefile}
#
#############################################################################
### PMIDs:
cat $tsvfile_gwas \
	|sed -e '1d' \
	|awk -F '\t' '{print $2}' \
	|sort -nu >$DATADIR/gwascat.pmid
###
# https://github.com/jeremyjyang/BioClients
python3 -m BioClients.icite.Client get_stats \
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
#
cat $DATADIR/gwascat_gwas.tsv \
	|sed -e '1d' |awk -F '\t' '{print $15}' \
	>$DATADIR/gwascat_gwas.gcst
###
# INFO:INPUT RCSTs: 5774; OUTPUT RCSTs: 4879 ; assns: 72366 ; loci: 72529 ; alleles: 73269 ; snps: 73269
# ~3hr
#${cwd}/python/gwascat_query.py \
python3 -m BioClients.gwascatalog.Client get_studyAssociations \
	--i $DATADIR/gwascat_gwas.gcst \
	--o $DATADIR/gwascat_StudyAssociations.tsv
gzip -f $DATADIR/gwascat_StudyAssociations.tsv
#
# re.match(r'(rs|SNP|snp|chr)', val)
cat $DATADIR/gwascat_snp2gene.tsv \
	|sed -e '1d' |awk -F '\t' '{print $3}' \
	|sort -u \
	>$DATADIR/gwascat_snp2gene.snpId
cat $DATADIR/gwascat_snp2gene.snpId |grep '^rs' \
	>$DATADIR/gwascat_snp2gene.rs
#
###
# ~15hr
#${cwd}/python/gwascat_query.py \
python3 -m BioClients.gwascatalog.Client get_snps \
	--i $DATADIR/gwascat_snp2gene.rs \
	--o $DATADIR/gwascat_Snps.tsv
gzip -f $DATADIR/gwascat_Snps.tsv
#
##
# (Split semicolon-separated multi-ENSG's.)
${cwd}/python/pandas_utils.py selectcols \
	--i $DATADIR/gwascat_Snps.tsv.gz \
	--coltags gene_ensemblGeneIds \
	|perl -pe 's/;/\n/g' \
	|sed -e '1d' |grep '^ENSG' |sort -u \
	>$DATADIR/gwascat_Snps.ensg
printf "Ensembl ID count: %d\n" "$(cat $DATADIR/gwascat_Snps.ensg |wc -l)"
#
###
# ~13hr
python3 -m BioClients.ensembl.Client get_info \
	--i $DATADIR/gwascat_Snps.ensg \
	--o $DATADIR/gwascat_Snps_EnsemblInfo.tsv
gzip -f $DATADIR/gwascat_Snps_EnsemblInfo.tsv
#
###
# tcrd_targets.tsv from:
# python3 -m BioClients.idg.tcrd.Client listTargets --dbname "tcrd610" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw=""
#############################################################################
# Gene-trait statistics:
# tiga_gt_stats.R: INPUT: 9 files; OUTPUT: gt_stats.tsv
#
${cwd}/R/tiga_gt_stats.R \
	$DATADIR/gwascat_gwas.tsv \
	$DATADIR/gwascat_counts.tsv \
	$DATADIR/gwascat_assn.tsv \
	$DATADIR/gwascat_snp2gene.tsv \
	$DATADIR/gwascat_trait.tsv \
	$DATADIR/gwascat_icite.tsv \
	$DATADIR/gwascat_Snps.tsv.gz \
	$DATADIR/gwascat_Snps_EnsemblInfo.tsv.gz \
	$DATADIR/tcrd_targets.tsv \
	$DATADIR/gt_stats.tsv
#

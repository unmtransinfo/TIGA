#!/bin/bash
#############################################################################
### Go_gwascat_GetData.sh - Create TSV files.
### gt_stats.tsv is used by TIGA web app.
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
#############################################################################
### Download releases from:
### ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/{YYYY}/{MM}/{DD}/
### Note that "v1.0.1", "v1.0.2", "v1.0.3" refer to formats, not releases.
#############################################################################
### According to EBI, additional information is available via the API,
### see https://www.ebi.ac.uk/gwas/docs/api,https://www.ebi.ac.uk/gwas/docs/faq,
### especially "Genomic Mappings" (SNP-gene). The API additional data includes
### Ensembl IDs for mapped genes, from which we query the Ensembl API
### for additional annotations including gene biotype, thereby filtering for
### protein_coding.
### Also it appears beta and OR values separated better.
#############################################################################
#
set -e
#
cwd=$(pwd)
#
SRCDATADIR="$HOME/../data/GWASCatalog/releases/2020/07/15"
DATADIR="${cwd}/data"
#
#Source files:
gwasfile="${SRCDATADIR}/gwas-catalog-studies_ontology-annotated.tsv"
#
#assnfile="${SRCDATADIR}/gwas_catalog_v1.0.2-associations_e94_r2018-09-30.tsv"
assnfile="${SRCDATADIR}/gwas-catalog-associations_ontology-annotated.tsv"
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
### TRAITS:
#
tsvfile_trait="${DATADIR}/gwascat_trait.tsv"
###
# EFO:
EFO_DIR="$HOME/../data/EFO/data"
OWLFILE="$EFO_DIR/efo.owl"
#EFO_URL="https://github.com/EBISPOT/efo/releases/download/v3.20.0/efo-base.owl"
EFO_URL="https://github.com/EBISPOT/efo/releases/download/v3.20.0/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
LIBDIR="$HOME/../app/lib"
efofile="${DATADIR}/efo.tsv"
###
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${efofile}
#
tsvfile_trait_sub="${DATADIR}/efo_sub_gwas.tsv"
###
${cwd}/R/gwascat_trait.R $gwasfile $efofile $tsvfile_trait $tsvfile_trait_sub
#
#############################################################################
### GENES:
#SNP to gene links:
snp2genefile="${DATADIR}/gwascat_snp2gene.tsv"
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
### Esp. "Genomic Mappings" (SNP-gene)
#
cat $DATADIR/gwascat_gwas.tsv \
	|sed -e '1d' |awk -F '\t' '{print $15}' \
	>$DATADIR/gwascat_gwas.gcst
###
# ~3hr
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
# TCRD:
python3 -m BioClients.idg.tcrd.Client listTargets \
	--dbname "tcrd660" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $DATADIR/tcrd_targets.tsv
#
###
# gwascat_counts.tsv from Go_gwascat_DbCreate.sh
###
# Pre-process and filter. Studies, genes and traits may be removed
# due to insufficient evidence.
${cwd}/R/tiga_gt_prepfilter.R \
	$DATADIR/gwascat_gwas.tsv \
	$DATADIR/gwascat_counts.tsv \
	$DATADIR/gwascat_assn.tsv \
	$DATADIR/gwascat_snp2gene.tsv \
	$DATADIR/gwascat_trait.tsv \
	$DATADIR/gwascat_icite.tsv \
	$DATADIR/gwascat_Snps.tsv.gz \
	$DATADIR/gwascat_Snps_EnsemblInfo.tsv.gz \
	$DATADIR/tcrd_targets.tsv \
	$DATADIR/gt_prepfilter.Rdata
###
# Provenance for gene-trait pairs (STUDY_ACCESSION, PUBMEDID).
${cwd}/R/tiga_gt_provenance.R \
	$DATADIR/gt_prepfilter.Rdata \
	$DATADIR/gt_provenance.tsv.gz
###
# Generates variables, statistics, evidence features for gene-trait pairs.
${cwd}/R/tiga_gt_variables.R \
	$DATADIR/gt_prepfilter.Rdata \
	$DATADIR/gt_variables.tsv.gz
###
# Scores and ranks gene-trait pairs based on selected variables.
${cwd}/R/tiga_gt_stats.R \
	$DATADIR/gt_variables.tsv.gz \
	$DATADIR/gt_stats.tsv.gz
#

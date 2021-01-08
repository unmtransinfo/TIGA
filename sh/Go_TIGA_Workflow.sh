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
DATADIR="${cwd}/data"
#
###
# GWASCatalog release:
GC_REL_Y="2020"
GC_REL_M="12"
GC_REL_D="16"
#
#SRCDATADIR="$HOME/../data/GWASCatalog/releases/2020/07/15"
SRCDATADIR="$HOME/../data/GWASCatalog/releases/${GC_REL_Y}/${GC_REL_M}/${GC_REL_D}"
#
printf "${GC_REL_Y}-${GC_REL_M}-${GC_REL_D}\n" >${DATADIR}/gwascat_release.txt
#
#Source files:
gwasfile="${SRCDATADIR}/gwas-catalog-studies_ontology-annotated.tsv"
if [ ! -f "${gwasfile}" ]; then
	echo "ERROR: FILE NOT FOUND: ${gwasfile}"
	exit
fi
#
assnfile="${SRCDATADIR}/gwas-catalog-associations_ontology-annotated.tsv"
if [ ! -f "${assnfile}" ]; then
	echo "ERROR: FILE NOT FOUND: ${assnfile}"
	exit
fi
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
EFO_RELEASE="3.25.0"
printf "${EFO_RELEASE}\n" >${DATADIR}/efo_release.txt
#
#EFO_URL="https://github.com/EBISPOT/efo/releases/download/v3.20.0/efo.owl"
EFO_URL="https://github.com/EBISPOT/efo/releases/download/v${EFO_RELEASE}/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
LIBDIR="$HOME/../app/lib"
efofile="${DATADIR}/efo.tsv"
###
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${efofile}
#
###
# From efo.tsv create GraphML file:
graphmlfile="${DATADIR}/efo_graph.graphml"
${cwd}/R/efo_graph.R ${efofile} ${graphmlfile}
gzip ${graphmlfile}
#
###
tsvfile_trait_sub="${DATADIR}/efo_sub_gwas.tsv"
#
${cwd}/R/gwascat_trait.R $gwasfile $efofile $tsvfile_trait $tsvfile_trait_sub
#
#############################################################################
### GENES:
#SNP to gene links:
snp2genefile="${DATADIR}/gwascat_snp2gene.tsv"
#
#############################################################################
### REPORTED GENES (ignored by TIGA):
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
printf "PMIDS: %d\n" $(cat $DATADIR/gwascat.pmid |wc -l)
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
	|sort -u \
	>$DATADIR/gwascat_upstream.ensg
cat $tsvfile_assn \
	|sed -e '1d' \
	|awk -F '\t' '{print $17}' \
	|egrep -v '(^$|^NA$)' \
	|sort -u \
	>$DATADIR/gwascat_downstream.ensg
cat $tsvfile_assn \
	|sed -e '1d' \
	|awk -F '\t' '{print $18}' \
	|egrep -v '(^$|^NA$)' \
	|perl -ne 'print join("\n",split(/, */))' \
	|sort -u \
	>$DATADIR/gwascat_snp.ensg
cat $DATADIR/gwascat_upstream.ensg $DATADIR/gwascat_downstream.ensg $DATADIR/gwascat_snp.ensg \
	|sort -u \
	>$DATADIR/gwascat.ensg
#
#############################################################################
### Additional information via the API:
### Esp. "Genomic Mappings" (SNP-gene)
### However, in 2020 there are ENSGs in downloads, so API not needed for ENSGs.
#
#cat $DATADIR/gwascat_gwas.tsv \
#	|sed -e '1d' |awk -F '\t' '{print $15}' \
#	>$DATADIR/gwascat_gwas.gcst
###
# ~3hr
#python3 -m BioClients.gwascatalog.Client get_studyAssociations \
#	--i $DATADIR/gwascat_gwas.gcst \
#	--o $DATADIR/gwascat_StudyAssociations.tsv
#gzip -f $DATADIR/gwascat_StudyAssociations.tsv
#
# re.match(r'(rs|SNP|snp|chr)', val)
cat $DATADIR/gwascat_snp2gene.tsv \
	|sed -e '1d' |awk -F '\t' '{print $3}' \
	|sort -u \
	>$DATADIR/gwascat_snp2gene.snpId
printf "SNPs: %d\n" $(cat $DATADIR/gwascat_snp2gene.snpId |wc -l)
cat $DATADIR/gwascat_snp2gene.snpId |grep '^rs' \
	>$DATADIR/gwascat_snp2gene.rs
printf "SNPs (rs*): %d\n" $(cat $DATADIR/gwascat_snp2gene.rs |wc -l)
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
# (No longer needed since ENSGs in download files.)
#python3 -m BioClients.util.pandas.Utils selectcols \ 
#	--i $DATADIR/gwascat_Snps.tsv.gz \
#	--coltags gene_ensemblGeneIds \
#	|perl -pe 's/;/\n/g' \
#	|sed -e '1d' |grep '^ENSG' |sort -u \
#	>$DATADIR/gwascat_Snps.ensg
#printf "Ensembl ID count: %d\n" "$(cat $DATADIR/gwascat_Snps.ensg |wc -l)"
#
###
# ~13hr
#python3 -m BioClients.ensembl.Client get_info \
#	--i $DATADIR/gwascat_Snps.ensg \
#	--o $DATADIR/gwascat_Snps_EnsemblInfo.tsv
#gzip -f $DATADIR/gwascat_Snps_EnsemblInfo.tsv
#
python3 -m BioClients.ensembl.Client get_info -v \
	--i $DATADIR/gwascat.ensg \
	--o $DATADIR/gwascat_EnsemblInfo.tsv
gzip -f $DATADIR/gwascat_EnsemblInfo.tsv
#
###
# TCRD:
TCRD_DBNAME="tcrd660"
python3 -m BioClients.idg.tcrd.Client listTargets \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $DATADIR/tcrd_targets.tsv
python3 -m BioClients.idg.tcrd.Client info \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $DATADIR/tcrd_info.tsv
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
	$DATADIR/gwascat_EnsemblInfo.tsv.gz \
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
cp \
	${DATADIR}/gwascat_gwas.tsv \
	${DATADIR}/filtered_studies.tsv \
	${DATADIR}/filtered_genes.tsv \
	${DATADIR}/filtered_traits.tsv \
	${DATADIR}/gt_provenance.tsv.gz \
	${DATADIR}/gt_stats.tsv.gz \
	${DATADIR}/efo_graph.graphml.gz \
	${DATADIR}/gwascat_release.txt \
	${DATADIR}/efo_release.txt \
	${DATADIR}/tcrd_info.tsv \
	${cwd}/R/tiga/data/

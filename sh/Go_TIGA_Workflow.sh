#!/bin/bash
#############################################################################
### Go_TIGA_Workflow.sh - TSVs for TIGA web app, DISEASES, TCRD.
#############################################################################
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
### ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/{YYYY}/{MM}/{DD}/
### Note that "v1.0.1", "v1.0.2", "v1.0.3" refer to formats, not releases.
#############################################################################
### Previously (2018), additional information is available via the API,
### "Genomic Mappings" with EnsemblIDs for mapped genes, now available via
### download assn file.
### From EnsemblIDs, we query Ensembl API for annotations including gene biotype,
### thereby filtering for protein_coding.
#############################################################################
# Install BioClients from https://github.com/jeremyjyang/BioClients
# or with "pip3 install BioClients".
#############################################################################
#
set -e
#
T0=$(date +%s)
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
###
# GWASCatalog release:
GC_REL_Y="2020"
#GC_REL_M="07"
GC_REL_M="12"
#GC_REL_D="15"
GC_REL_D="16"
#
ODIR="${DATADIR}/${GC_REL_Y}${GC_REL_M}${GC_REL_D}"
if [ ! -d $ODIR ]; then
	mkdir -p $ODIR
fi
#
SRCDIR="$HOME/../data/GWASCatalog/releases/${GC_REL_Y}/${GC_REL_M}/${GC_REL_D}"
#
printf "${GC_REL_Y}-${GC_REL_M}-${GC_REL_D}\n" >${ODIR}/gwascat_release.txt
#
#Source files:
gwasfile="${SRCDIR}/gwas-catalog-studies_ontology-annotated.tsv"
if [ ! -f "${gwasfile}" ]; then
	echo "ERROR: FILE NOT FOUND: ${gwasfile}"
	exit
fi
#
assnfile="${SRCDIR}/gwas-catalog-associations_ontology-annotated.tsv"
if [ ! -f "${assnfile}" ]; then
	echo "ERROR: FILE NOT FOUND: ${assnfile}"
	exit
fi
###
#Output files:
tsvfile_gwas="${ODIR}/gwascat_gwas.tsv"
tsvfile_assn="${ODIR}/gwascat_assn.tsv"
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
tsvfile_trait="${ODIR}/gwascat_trait.tsv"
###
# EFO:
EFO_DIR="$HOME/../data/EFO/data"
OWLFILE="$EFO_DIR/efo.owl"
#EFO_RELEASE="3.20.0"
EFO_RELEASE="3.25.0"
printf "${EFO_RELEASE}\n" >${ODIR}/efo_release.txt
#
EFO_URL="https://github.com/EBISPOT/efo/releases/download/v${EFO_RELEASE}/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
LIBDIR="$HOME/../app/lib"
efofile="${ODIR}/efo.tsv"
###
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${efofile}
#
###
tsvfile_trait_sub="${ODIR}/efo_sub_gwas.tsv"
#
${cwd}/R/gwascat_trait.R $gwasfile $efofile $tsvfile_trait $tsvfile_trait_sub
#
###
# From efo.tsv create GraphML file:
graphmlfile="${ODIR}/efo_graph.graphml"
${cwd}/R/efo_graph.R ${efofile} ${tsvfile_trait_sub} ${graphmlfile}
gzip -f ${graphmlfile}
#
#
#############################################################################
### GENES:
#SNP to gene links:
snp2genefile="${ODIR}/gwascat_snp2gene.tsv"
#
#############################################################################
### REPORTED GENES (ignored by TIGA):
#
printf "STUDY_ACCESSION\tSNP\tGSYMB\tENSG\tREPORTED_OR_MAPPED\n" >${snp2genefile}
#
# "REPORTED_GENE(S),SNPS,STUDY_ACCESSION" (14, 22, 37)
###
cat $tsvfile_assn |sed -e '1d' \
	|perl -n perl/snp2gene_reported.pl \
	>>${snp2genefile}
#
#############################################################################
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
# UPSTREAM_GENE_ID,DOWNSTREAM_GENE_ID,SNP_GENE_IDS,SNPS,STUDY_ACCESSION (16,17,18,22,37)
###
cat $tsvfile_assn |sed -e '1d' \
	|perl -n perl/snp2gene_mapped.pl \
	>>${snp2genefile}
#
#############################################################################
### PMIDs:
cat $tsvfile_gwas \
	|sed -e '1d' \
	|awk -F '\t' '{print $2}' \
	|sort -nu >$ODIR/gwascat.pmid
printf "PMIDS: %d\n" $(cat $ODIR/gwascat.pmid |wc -l)
###
if [ ! -f "$ODIR/gwascat_icite.tsv" ]; then
	python3 -m BioClients.icite.Client get_stats \
		--i $ODIR/gwascat.pmid \
		--o $ODIR/gwascat_icite.tsv
fi
#
#############################################################################
### Entrez gene IDs: UPSTREAM_GENE_ID, DOWNSTREAM_GENE_ID, SNP_GENE_IDS
if [ ! -e $ODIR/gwascat_EnsemblInfo.tsv.gz ]; then
	cat $tsvfile_assn |sed -e '1d' \
		|awk -F '\t' '{print $16}' \
		|egrep -v '(^$|^NA$)' \
		|sort -u \
		>$ODIR/gwascat_upstream.ensg
	cat $tsvfile_assn |sed -e '1d' \
		|awk -F '\t' '{print $17}' \
		|egrep -v '(^$|^NA$)' \
		|sort -u \
		>$ODIR/gwascat_downstream.ensg
	cat $tsvfile_assn |sed -e '1d' \
		|awk -F '\t' '{print $18}' \
		|egrep -v '(^$|^NA$)' \
		|perl -ne 'print join("\n",split(/, */))' \
		|sort -u \
		>$ODIR/gwascat_snp.ensg
	cat $ODIR/gwascat_upstream.ensg $ODIR/gwascat_downstream.ensg $ODIR/gwascat_snp.ensg \
		|sort -u \
		>$ODIR/gwascat.ensg
	#
	# ~13hr
	python3 -m BioClients.ensembl.Client get_info -v \
		--i $ODIR/gwascat.ensg |gzip -c \
		>$ODIR/gwascat_EnsemblInfo.tsv.gz
fi
#
###
# TCRD:
TCRD_DBNAME="tcrd684"
python3 -m BioClients.idg.tcrd.Client listTargets \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $ODIR/tcrd_targets.tsv
python3 -m BioClients.idg.tcrd.Client info \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $ODIR/tcrd_info.tsv
#
###
# Go_TIGA_DbCreate.sh
# input:
#	gwascat_gwas.tsv
#	gwascat_assn.tsv
#	gwascat_snp2gene.tsv
#	gwascat_trait.tsv
#	gwascat_icite.tsv
# output:
#	gwascat_counts.tsv
${cwd}/sh/Go_TIGA_DbCreate.sh ${ODIR}
#
###
# Pre-process and filter. Studies, genes and traits may be removed
# due to insufficient evidence.
# $ODIR/gwascat_Snps.tsv.gz (NO LONGER NEEDED)
${cwd}/R/tiga_gt_prepfilter.R \
	$ODIR/gwascat_gwas.tsv \
	$ODIR/gwascat_counts.tsv \
	$ODIR/gwascat_assn.tsv \
	$ODIR/gwascat_snp2gene.tsv \
	$ODIR/gwascat_trait.tsv \
	$ODIR/gwascat_icite.tsv \
	$ODIR/gwascat_EnsemblInfo.tsv.gz \
	$ODIR/tcrd_targets.tsv \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/filtered_studies.tsv \
	$ODIR/filtered_traits.tsv \
	$ODIR/filtered_genes.tsv
###
# Provenance for gene-trait pairs (STUDY_ACCESSION, PUBMEDID).
${cwd}/R/tiga_gt_provenance.R \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/gt_provenance.tsv.gz
###
# Generates variables, statistics, evidence features for gene-trait pairs.
${cwd}/R/tiga_gt_variables.R \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/gt_variables.tsv.gz
###
# Scores and ranks gene-trait pairs based on selected variables.
${cwd}/R/tiga_gt_stats.R \
	$ODIR/gt_variables.tsv.gz \
	$ODIR/gt_stats.tsv.gz
#
cp \
	${ODIR}/gwascat_gwas.tsv \
	${ODIR}/filtered_studies.tsv \
	${ODIR}/filtered_genes.tsv \
	${ODIR}/filtered_traits.tsv \
	${ODIR}/gt_provenance.tsv.gz \
	${ODIR}/gt_stats.tsv.gz \
	${ODIR}/efo_graph.graphml.gz \
	${ODIR}/tcrd_info.tsv \
	${ODIR}/gwascat_release.txt \
	${ODIR}/efo_release.txt \
	${cwd}/R/tiga/data/
#
printf "Elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
#

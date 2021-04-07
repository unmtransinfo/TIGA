#!/bin/bash
#############################################################################
### Go_TIGA_Workflow_GCAPI.sh - Variant using catalog API.
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
### Issue: API provides many more "Ensembl_pipeline" gene mappings than download
### file, so we merge both.
#############################################################################
### From EnsemblIDs, we query Ensembl API for annotations including gene biotype,
### including "protein_coding", but prefer TCRD mappings to define protein
### coding.
#############################################################################
# Install BioClients from https://github.com/jeremyjyang/BioClients
# or with "pip3 install BioClients".
#############################################################################
# Issue: ENSEMBL API USE (SLOW): Streamline via FTP comprehensive geneslist.
# http://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.103.entrez.tsv.gz
#############################################################################
#
set -e
#
T0=$(date +%s)
#
cwd=$(pwd)
#
GWASCATALOGDIR="$(cd $HOME/../data/GWASCatalog; pwd)"
DATADIR="${cwd}/data"
###
# GWASCatalog release:
#GC_REL="2020-07-15"
#GC_REL="2020-12-16"
#GC_REL="2021-02-12"
#GC_REL="2021-03-29"
if [ $# -eq 1 ]; then
	GC_REL=$1
else
	printf "ERROR: syntax $(basename $0) \"YYYY-MM-DD\"\n"
	printf "$(date +'%Y') releases:\n"
	ls -l ${GWASCATALOGDIR}/releases/$(date +'%Y')/*
	exit
fi
#
GC_REL_Y=$(echo $GC_REL |sed 's/-.*$//')
GC_REL_M=$(echo $GC_REL |sed 's/^.*-\(.*\)-.*$/\1/')
GC_REL_D=$(echo $GC_REL |sed 's/^.*-//')
#
printf "GWASCatalog release: \"%s\" (Y=%s,M=%s,D=%s)\n" "$GC_REL" "$GC_REL_Y" "$GC_REL_M" "$GC_REL_D"
#
if [ ! "$GC_REL_Y" -o  ! "$GC_REL_M" -o  ! "$GC_REL_D" ]; then
	printf "ERROR: Badly formed GWASCatalog release (YYYY-MM-DD): \"%s\"\n" "$GC_REL"
	exit
fi 
#
ODIR="${DATADIR}/${GC_REL_Y}${GC_REL_M}${GC_REL_D}"
if [ ! -d $ODIR ]; then
	mkdir -p $ODIR
fi
#exit #DEBUG
#ODIR="$DATADIR/20210212_GCAPI" #DEBUG
#
SRCDIR="$GWASCATALOGDIR/releases/${GC_REL_Y}/${GC_REL_M}/${GC_REL_D}"
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
#SNP to gene links, from download association file:
snp2genefile_file="${ODIR}/gwascat_snp2gene_FILE.tsv"
#
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
### (REPORTED GENES not used for TIGA scoring.)
#
${cwd}/python/snp2gene.py $tsvfile_assn --o ${snp2genefile_file}
#
###
# (Alternative to download file; may be incomplete/different.)
# SNPs, SNP2GENE, via API:
cat $tsvfile_assn |sed -e '1d' \
	|awk -F '\t' '{print $22}' \
	|perl -pe 's/[; ]+/\n/g' \
	|perl -pe 's/ x /\n/g' \
	|grep '^rs' \
	|sort -u \
	>${ODIR}/gwascat_snp.rs
printf "SNPs: %d\n" $(cat $ODIR/gwascat_snp.rs |wc -l)
python3 -m BioClients.gwascatalog.Client get_snps -q \
	--i ${ODIR}/gwascat_snp.rs \
	--o ${ODIR}/gwascat_snp_API.tsv
#
#SNP to gene links, from API:
snp2genefile_api="${ODIR}/gwascat_snp2gene_API.tsv"
python3 -m BioClients.util.pandas.Utils selectcols \
	--i $ODIR/gwascat_snp_API.tsv \
	--coltags "rsId,isIntergenic,isUpstream,isDownstream,distance,source,mappingMethod,isClosestGene,chromosomeName,chromosomePosition,geneName,ensemblGeneIds" \
	--o ${snp2genefile_api}
#
#Merge FILE and API snp2gene files:
snp2genefile_merged="${ODIR}/gwascat_snp2gene_MERGED.tsv"
${cwd}/R/snp2gene_merge.R \
	${snp2genefile_file} \
	${snp2genefile_api} \
	${snp2genefile_merged}
#
#############################################################################
# http://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.103.entrez.tsv.gz
ENTREZGENEFILE="Homo_sapiens.GRCh38.103.entrez.tsv.gz"
if [ ! -e $ODIR/gwascat_EnsemblInfo.tsv.gz ]; then
	wget -O - "http://ftp.ensembl.org/pub/current_tsv/homo_sapiens/$ENTREZGENEFILE" >$ODIR/$ENTREZGENEFILE
	gunzip -c $ODIR/$ENTREZGENEFILE |sed '1d' |awk -F '\t' '{print $1}' |sort -u \
		>$ODIR/ensembl_human_genes.ensg
	python3 -m BioClients.ensembl.Client get_info -q \
		--i $ODIR/ensembl_human_genes.ensg |gzip -c \
		>$ODIR/gwascat_EnsemblInfo.tsv.gz
fi
#############################################################################
### PMIDs:
cat $tsvfile_gwas \
	|sed -e '1d' |awk -F '\t' '{print $2}' |sort -nu \
	>$ODIR/gwascat.pmid
printf "PMIDS: %d\n" $(cat $ODIR/gwascat.pmid |wc -l)
###
if [ ! -f "$ODIR/gwascat_icite.tsv" ]; then
	python3 -m BioClients.icite.Client get_stats -q \
		--i $ODIR/gwascat.pmid \
		--o $ODIR/gwascat_icite.tsv
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
# Generate counts via MySql:
${cwd}/sh/Go_TIGA_DbCreate.sh \
	"${ODIR}" \
	"tiga_${GC_REL_Y}${GC_REL_M}${GC_REL_D}" \
	${ODIR}/gwascat_gwas.tsv \
	${ODIR}/gwascat_assn.tsv \
	${snp2genefile_merged} \
	${ODIR}/gwascat_trait.tsv \
	${ODIR}/gwascat_icite.tsv \
	${ODIR}/gwascat_counts.tsv \
	${ODIR}/trait_counts.tsv
#
###
# Pre-process and filter. Studies, genes and traits may be removed
# due to insufficient evidence.
${cwd}/R/tiga_gt_prepfilter.R \
	$ODIR/gwascat_gwas.tsv \
	$ODIR/gwascat_counts.tsv \
	$ODIR/gwascat_assn.tsv \
	${snp2genefile_merged} \
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
# Mu scores for benchmark comparision.
${cwd}/python/tiga_gt_stats_mu.py --mutags "pvalue_mlog_max,rcras,n_snpw" \
	-q \
	--i $ODIR/gt_variables.tsv.gz \
	--o $ODIR/gt_stats_mu.tsv.gz
###
printf "Copy files for TIGA web app with command:\n"
printf "cp ${ODIR}/gwascat_gwas.tsv ${ODIR}/filtered_studies.tsv ${ODIR}/filtered_genes.tsv ${ODIR}/filtered_traits.tsv ${ODIR}/gt_provenance.tsv.gz ${ODIR}/gt_stats.tsv.gz ${ODIR}/efo_graph.graphml.gz ${ODIR}/tcrd_info.tsv ${ODIR}/gwascat_release.txt ${ODIR}/efo_release.txt ${cwd}/R/tiga/data/\n"
printf "Remove TIGA web app Rmd with command:\n"
printf "rm -f ${cwd}/R/tiga/tiga.Rmd\n"
#
printf "Elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
#

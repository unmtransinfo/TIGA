#!/bin/bash
#############################################################################
### Go_TIGA_Workflow.sh - Produces TSVs for TIGA web app, DISEASES, TCRD.
#############################################################################
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
### ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/{YYYY}/{MM}/{DD}/
### Note that "v1.0.1", "v1.0.2", "v1.0.3" refer to formats, not releases.
#############################################################################
### Using catalog API, for additional SNP-gene mappings.
### Genomic Mappings with EnsemblIDs for mapped genes, available via
### download assn file, but API provides many more "Ensembl_pipeline"
### gene mappings than download file, so we currently merge both.
### Issue: Catalog API slow (~50hrs for 154656 SNPs)
#############################################################################
### From EnsemblIDs, we query Ensembl API for annotations including gene biotype,
### including "protein_coding", but prefer TCRD mappings to define protein
### coding.
#############################################################################
# Dependency: https://github.com/jeremyjyang/BioClients
# ("pip3 install BioClients")
#############################################################################
#
set -e
#
function MessageBreak {
  printf "============================================\n"
  printf "=== [%s] %s\n" "$(date +'%Y-%m-%d:%H:%M:%S')" "$1"
}
#
T0=$(date +%s)
#
cwd=$(pwd)
#
GWASCATALOGDIR="$(cd $HOME/../data/GWASCatalog; pwd)"
DATADIR="${cwd}/data"
###
MessageBreak "Starting $(basename $0)"
###
# GWASCatalog release:
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
### OUTPUT FILES:
tsvfile_gwas="${ODIR}/gwascat_gwas.tsv"
tsvfile_assn="${ODIR}/gwascat_assn.tsv"
efofile="${ODIR}/efo.tsv"
tsvfile_trait="${ODIR}/gwascat_trait.tsv"
tsvfile_trait_sub="${ODIR}/efo_sub_gwas.tsv"
tsvfile_icite="${ODIR}/gwascat_icite.tsv"
snp2genefile_file="${ODIR}/gwascat_snp2gene_FILE.tsv"
snp2genefile_api="${ODIR}/gwascat_snp2gene_API.tsv"
snp2genefile_merged="${ODIR}/gwascat_snp2gene_MERGED.tsv"
###
MessageBreak "Clean studies:"
#Clean studies:
${cwd}/R/gwascat_gwas.R $gwasfile $tsvfile_gwas
#
###
MessageBreak "Clean associations:"
#Clean, separate OR_or_beta into oddsratio, beta columns:
${cwd}/R/gwascat_assn.R $assnfile $tsvfile_assn
#
#############################################################################
### TRAITS:
#
MessageBreak "TRAITS:"
###
# EFO:
EFO_DIR="$(cd $HOME/../data/EFO/data; pwd)"
OWLFILE="$EFO_DIR/efo.owl"
#EFO_RELEASE="3.20.0"
EFO_RELEASE="3.25.0"
printf "${EFO_RELEASE}\n" >${ODIR}/efo_release.txt
#
EFO_URL="https://github.com/EBISPOT/efo/releases/download/v${EFO_RELEASE}/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
LIBDIR="$(cd $HOME/../app/lib; pwd)"
###
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${efofile}
#
###
#
MessageBreak "Clean traits:"
${cwd}/R/gwascat_trait.R $gwasfile $efofile $tsvfile_trait $tsvfile_trait_sub
#
###
MessageBreak "Create EFO GraphML file:"
# From efo.tsv create GraphML file:
graphmlfile="${ODIR}/efo_graph.graphml"
${cwd}/R/efo_graph.R ${efofile} ${tsvfile_trait_sub} ${graphmlfile}
gzip -f ${graphmlfile}
#
#
#############################################################################
### GENES:
MessageBreak "GENES:"
#
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
### (REPORTED GENES not used for TIGA scoring.)
#
MessageBreak "SNP2GENE (from association file):"
#SNP to gene links, from download association file:
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
MessageBreak "GWASCATALOG API REQUESTS (get_snps):"
python3 -m BioClients.gwascatalog.Client get_snps -q \
	--i ${ODIR}/gwascat_snp.rs \
	--o ${ODIR}/gwascat_snp_API.tsv
#
#SNP2GENE, from API:
python3 -m BioClients.util.pandas.Utils selectcols \
	--i $ODIR/gwascat_snp_API.tsv \
	--coltags "rsId,isIntergenic,isUpstream,isDownstream,distance,source,mappingMethod,isClosestGene,chromosomeName,chromosomePosition,geneName,ensemblGeneIds" \
	--o ${snp2genefile_api}
#
#Merge FILE and API snp2gene files:
${cwd}/R/snp2gene_merge.R \
	${snp2genefile_file} \
	${snp2genefile_api} \
	${snp2genefile_merged}
#
#############################################################################
# Download latest Ensembl human gene file.
# ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/
#ENTREZGENEFILE="Homo_sapiens.GRCh38.104.entrez.tsv.gz"
ENTREZGENEFILE=$(lftp ftp://anonymous:@ftp.ensembl.org -e "cd pub/current_tsv/homo_sapiens/ ; ls *.entrez.tsv.gz; quit" |sed 's/^.* //')
printf "ENTREZGENEFILE: %s\n" "${ENTREZGENEFILE}"
ensemblinfofile="$ODIR/gwascat_EnsemblInfo.tsv"
if [ ! -e ${ensemblinfofile} ]; then
	wget -O - "ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/$ENTREZGENEFILE" >$ODIR/$ENTREZGENEFILE
	gunzip -c $ODIR/$ENTREZGENEFILE |sed '1d' |awk -F '\t' '{print $1}' |sort -u \
		>$ODIR/ensembl_human_genes.ensg
	MessageBreak "ENSEMBL API REQUESTS (get_info):"
	python3 -m BioClients.ensembl.Client get_info -q \
		--i $ODIR/ensembl_human_genes.ensg \
		--o ${ensemblinfofile}
fi
#############################################################################
### PMIDs:
MessageBreak "PUBLICATIONS (iCite):"
cat $tsvfile_gwas \
	|sed -e '1d' |awk -F '\t' '{print $2}' |sort -nu \
	>$ODIR/gwascat.pmid
printf "PMIDS: %d\n" $(cat $ODIR/gwascat.pmid |wc -l)
###
if [ ! -f "${tsvfile_icite}" ]; then
	python3 -m BioClients.icite.Client get_stats -q \
		--i $ODIR/gwascat.pmid \
		--o ${tsvfile_icite}
fi
#
###
# TCRD:
MessageBreak "IDG (TCRD):"
TCRD_DBNAME="tcrd6110"
python3 -m BioClients.idg.tcrd.Client listTargets \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $ODIR/tcrd_targets.tsv
python3 -m BioClients.idg.tcrd.Client info \
	--dbname "${TCRD_DBNAME}" --dbhost="tcrd.kmc.io" --dbusr="tcrd" --dbpw="" \
	--o $ODIR/tcrd_info.tsv
#
###
MessageBreak "Generate counts:"
# Generate counts via Python:
${cwd}/python/tiga_gwas_counts.py \
	--ifile_gwas ${tsvfile_gwas} \
	--ifile_assn ${tsvfile_assn} \
	--ifile_trait ${tsvfile_trait} \
	--ifile_snp2gene ${snp2genefile_merged} \
	--ifile_icite ${tsvfile_icite} \
	--ofile_gwas $ODIR/gwascat_gwas_counts.tsv \
	--ofile_trait $ODIR/gwascat_trait_counts.tsv
#
###
MessageBreak "PREPFILTER:"
# Pre-process and filter. Studies, genes and traits may be removed
# due to insufficient evidence.
${cwd}/R/tiga_gt_prepfilter.R \
	${tsvfile_gwas} \
	$ODIR/gwascat_gwas_counts.tsv \
	${tsvfile_assn} \
	${snp2genefile_merged} \
	${tsvfile_trait} \
	${tsvfile_icite} \
	${ensemblinfofile} \
	$ODIR/tcrd_targets.tsv \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/filtered_studies.tsv \
	$ODIR/filtered_traits.tsv \
	$ODIR/filtered_genes.tsv
###
MessageBreak "PROVENANCE:"
# Provenance for gene-trait pairs (STUDY_ACCESSION, PUBMEDID).
${cwd}/R/tiga_gt_provenance.R \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/gt_provenance.tsv.gz
###
MessageBreak "VARIABLES:"
# Generates variables, statistics, evidence features for gene-trait pairs.
${cwd}/R/tiga_gt_variables.R \
	$ODIR/gt_prepfilter.Rdata \
	$ODIR/gt_variables.tsv.gz
###
MessageBreak "STATS:"
# Scores and ranks gene-trait pairs based on selected variables.
${cwd}/R/tiga_gt_stats.R \
	$ODIR/gt_variables.tsv.gz \
	$ODIR/gt_stats.tsv.gz
###
MessageBreak "STATS_MU:"
# Mu scores for benchmark comparision.
${cwd}/python/tiga_gt_stats_mu.py --mutags "pvalue_mlog_max,rcras,n_snpw" \
	-q \
	--i $ODIR/gt_variables.tsv.gz \
	--o $ODIR/gt_stats_mu.tsv.gz
###
printf "Copy files for TIGA web app with command:\n"
printf "cp ${ODIR}/gwascat_gwas.tsv ${ODIR}/filtered_studies.tsv ${ODIR}/filtered_genes.tsv ${ODIR}/filtered_traits.tsv ${ODIR}/gt_provenance.tsv.gz ${ODIR}/gt_stats.tsv.gz ${ODIR}/efo_graph.graphml.gz ${ODIR}/tcrd_info.tsv ${ODIR}/gwascat_release.txt ${ODIR}/efo_release.txt ${cwd}/R/tiga/data/\n"
printf "Remove TIGA web app Rdata with command:\n"
printf "rm -f ${cwd}/R/tiga/tiga.Rdata\n"
#
###
printf "Also remember to copy TIGA download files to: https://unmtid-shinyapps.net/download/TIGA/${GC_REL_Y}${GC_REL_M}${GC_REL_D} with symlink https://unmtid-shinyapps.net/download/TIGA/latest.\n"
#
printf "Elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
MessageBreak "Done $(basename $0)"
#

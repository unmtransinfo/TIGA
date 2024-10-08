#!/bin/bash
#############################################################################
### Go_TIGA_Workflow.sh - Produces TSVs for TIGA web app, DISEASES, TCRD.
#############################################################################
### NHGRI-EBI GWAS Catalog: http://www.ebi.ac.uk/gwas/
### ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/{YYYY}/{MM}/{DD}/
### Note that "v1.0.1", "v1.0.2", "v1.0.3" refer to formats, not releases.
#############################################################################
### Need to manually check https://github.com/EBISPOT/efo/releases
### for latest EFO release prior to GWC release.
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
if [ -f "${cwd}/LATEST_RELEASE_GWC.txt" ]; then
	GC_REL=$(cat ${cwd}/LATEST_RELEASE_GWC.txt)
else
	printf "ERROR: not found: ${cwd}/LATEST_RELEASE_GWC.txt\n"
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
###
# Should be latest EFO release prior to GWC release.
if [ -f "${cwd}/LATEST_RELEASE_EFO.txt" ]; then
	EFO_RELEASE=$(cat ${cwd}/LATEST_RELEASE_EFO.txt)
else
	printf "ERROR: not found: ${cwd}/LATEST_RELEASE_EFO.txt\n"
	exit
fi
printf "EFO release: ${EFO_RELEASE}\n"
printf "${EFO_RELEASE}\n" >${ODIR}/efo_release.txt
#
###
#Source files:
#FILESET CHANGE: 2024-01-22 (see email from Elliot Sollis)
#gwasfile="${SRCDIR}/gwas-catalog-studies_ontology-annotated.tsv"
gwasfile="${SRCDIR}/gwas-catalog-studies-download-alternative-v1.0.2.1.txt"
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
# Activate Virtual Environment
source ${cwd}/venv/bin/activate
###
# TCRD:
# Version specified here:
TCRD_DBNAME="tcrd6134pharos2"
TCRD_DBHOST="tcrd.ncats.io"
TCRD_DBUSR="tcrd"
TCRD_DBPW=""
#
MessageBreak "IDG (TCRD):"
if [ ! -s $ODIR/tcrd_targets.tsv ]; then
python3 -m BioClients.idg.tcrd.Client listTargets \
	--dbname "${TCRD_DBNAME}" --dbhost="${TCRD_DBHOST}" --dbusr="${TCRD_DBUSR}" --dbpw="${TCRD_DBPW}" \
	--o $ODIR/tcrd_targets.tsv
fi
if [ ! -s $ODIR/tcrd_info.tsv ]; then
python3 -m BioClients.idg.tcrd.Client info \
	--dbname "${TCRD_DBNAME}" --dbhost="${TCRD_DBHOST}" --dbusr="${TCRD_DBUSR}" --dbpw="${TCRD_DBPW}" \
	--o $ODIR/tcrd_info.tsv
fi
#
###
### OUTPUT FILES:
tsvfile_gwas="${ODIR}/gwascat_gwas.tsv" #gwascat_gwas.R
tsvfile_assn="${ODIR}/gwascat_assn.tsv" #gwascat_assn.R
efofile="${ODIR}/efo.tsv"
tsvfile_trait="${ODIR}/gwascat_trait.tsv" #gwascat_trait.R
tsvfile_trait_sub="${ODIR}/efo_sub_gwas.tsv" #gwascat_trait.R
tsvfile_icite="${ODIR}/gwascat_icite.tsv"
snp2genefile_file="${ODIR}/gwascat_snp2gene_FILE.tsv"
snpfile_api="${ODIR}/gwascat_snp_API.tsv"
snp2genefile_api="${ODIR}/gwascat_snp2gene_API.tsv"
snp2genefile_merged="${ODIR}/gwascat_snp2gene_MERGED.tsv"
###
#Clean studies:
MessageBreak "Clean studies:"
${cwd}/R/gwascat_gwas.R
#
###
#Clean, separate OR_or_beta into oddsratio, beta columns:
MessageBreak "Clean associations:"
${cwd}/R/gwascat_assn.R
#
#############################################################################
### TRAITS:
#
MessageBreak "TRAITS:"
###
# EFO:
EFO_DIR="$(cd $HOME/../data/EFO/data; pwd)"
OWLFILE="$EFO_DIR/efo.owl"
###
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
${cwd}/R/gwascat_trait.R
#
###
# From efo.tsv create GraphML file:
MessageBreak "Create EFO GraphML file:"
graphmlfile="${ODIR}/efo_graph.graphml" #efo_graph.R
${cwd}/R/efo_graph.R $GC_REL_Y $GC_REL_M $GC_REL_D
gzip -f ${graphmlfile}
#
#
#############################################################################
### GENES:
#
### MAPPED GENES:
### Separate mapped into up-/down-stream.
# "m" - mapped within gene
# "mu" - mapped to upstream gene
# "md" - mapped to downstream gene
### (REPORTED GENES not used for TIGA scoring.)
MessageBreak "GENES:"
#
#SNP to gene links, from download association file:
MessageBreak "SNP2GENE (from association file):"
${cwd}/python/snp2gene.py --o ${snp2genefile_file} $tsvfile_assn
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
if [ -f "${snpfile_api}" ]; then
	printf "File exists, not regenerated: %s (May have required manual effort due to API issues.)\n" ${snpfile_api}
else
	#python3 -m BioClients.gwascatalog.Client get_snps \
	python3 -m BioClients.gwascatalog.Client get_snps -q \
		--i ${ODIR}/gwascat_snp.rs \
		--o ${snpfile_api}
fi
#
#SNP2GENE, from API:
python3 -m BioClients.util.pandas.App selectcols \
	--i ${snpfile_api} \
	--coltags "rsId,isIntergenic,isUpstream,isDownstream,distance,source,mappingMethod,isClosestGene,chromosomeName,chromosomePosition,geneName,ensemblGeneIds" \
	--o ${snp2genefile_api}
#
#Merge FILE and API snp2gene files:
${cwd}/R/snp2gene_merge.R
#
#############################################################################
# Download latest Ensembl human gene file.
# ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/
#ENTREZGENEFILE="Homo_sapiens.GRCh38.104.entrez.tsv.gz"
ENTREZGENEFILE=$(lftp ftp://anonymous:@ftp.ensembl.org -e "cd pub/current_tsv/homo_sapiens/ ; ls *.entrez.tsv.gz; quit" |sed 's/^.* //')
printf "ENTREZGENEFILE: %s\n" "${ENTREZGENEFILE}"
ensemblinfofile="$ODIR/gwascat_EnsemblInfo.tsv"
if [ ! -s ${ensemblinfofile} ]; then
	wget -O - "ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/$ENTREZGENEFILE" >$ODIR/$ENTREZGENEFILE
	gunzip -c $ODIR/$ENTREZGENEFILE |sed '1d' |awk -F '\t' '{print $1}' |sort -u \
		>$ODIR/ensembl_human_genes.ensg
	MessageBreak "ENSEMBL API REQUESTS (get_info):"
	#python3 -m BioClients.ensembl.Client get_info -q \
	python3 -m BioClients.ensembl.Client get_info \
		--i $ODIR/ensembl_human_genes.ensg \
		--o ${ensemblinfofile}
else
	printf "File exists, not regenerated: %s (May have required manual effort due to API issues.)\n" ${ensemblinfofile}
fi
#############################################################################
### PMIDs:
MessageBreak "PUBLICATIONS (iCite):"
###
if [ ! -s "${tsvfile_icite}" ]; then
	cat $tsvfile_gwas \
		|sed -e '1d' |awk -F '\t' '{print $2}' |sort -nu \
		>$ODIR/gwascat.pmid
	printf "PMIDS: %d\n" $(cat $ODIR/gwascat.pmid |wc -l)
	#python3 -m BioClients.icite.Client get_stats -q \
	python3 -m BioClients.icite.Client get_stats \
		--i $ODIR/gwascat.pmid \
		--o ${tsvfile_icite}
else
	printf "File exists, not regenerated: %s (May have required manual effort due to API issues.)\n" ${tsvfile_icite}
fi
#
###
# Generate counts via Python:
MessageBreak "Generate counts:"
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
# Pre-process and filter. Studies, genes and traits may be removed
# due to insufficient evidence.
MessageBreak "PREPFILTER:"
${cwd}/R/tiga_gt_prepfilter.R
###
# Provenance for gene-trait pairs (STUDY_ACCESSION, PUBMEDID).
MessageBreak "PROVENANCE:"
${cwd}/R/tiga_gt_provenance.R
###
# Generates variables, statistics, evidence features for gene-trait pairs.
MessageBreak "VARIABLES:"
${cwd}/R/tiga_gt_variables.R
###
# Scores and ranks gene-trait pairs based on selected variables.
MessageBreak "STATS:"
${cwd}/R/tiga_gt_stats.R
###
# Mu scores for benchmark comparision. (Not needed for routine updates.)
#MessageBreak "STATS_MU:"
#${cwd}/python/tiga_gt_stats_mu.py --mutags "pvalue_mlog_max,rcras,n_snpw" \
#	-q \
#	--i $ODIR/gt_variables.tsv.gz \
#	--o $ODIR/gt_stats_mu.tsv.gz
###
# Generate final output files.
# INPUTS: gwascat_gwas.tsv filtered_studies.tsv filtered_genes.tsv filtered_traits.tsv gt_provenance.tsv.gz gt_stats.tsv.gz efo_graph.graphml.gz tcrd_info.tsv gwascat_release.txt efo_release.txt
# OUTPUTS: tiga.Rdata tiga_gene-trait_stats.tsv tiga_gene-trait_provenance.tsv tiga_genes.tsv tiga_traits.tsv
#
MessageBreak "FINAL OUTPUT FILES:"
${cwd}/R/tiga_final_files.R
###
# Show commands for installing updated datafiles.
printf "Copy Rdata for TIGA web app with command:\n"
printf "cp ${ODIR}/tiga.Rdata ${cwd}/R/tiga/\n"
printf "Copy TIGA download files to: unmtid-dbs.net/download/TIGA/${GC_REL_Y}${GC_REL_M}${GC_REL_D} with revised symlink \"latest\":\n"
printf "ssh unmtid-dbs.net mkdir /var/www/html/download/TIGA/${GC_REL_Y}${GC_REL_M}${GC_REL_D}\n"
printf "scp $ODIR/tiga_gene-trait_stats.tsv $ODIR/tiga_gene-trait_provenance.tsv $ODIR/tiga_genes.tsv $ODIR/tiga_traits.tsv unmtid-dbs.net:/var/www/html/download/TIGA/${GC_REL_Y}${GC_REL_M}${GC_REL_D}\n"
printf "ssh unmtid-dbs.net rm /var/www/html/download/TIGA/latest\n"
printf "ssh unmtid-dbs.net ln -s /var/www/html/download/TIGA/${GC_REL_Y}${GC_REL_M}${GC_REL_D} /var/www/html/download/TIGA/latest\n"
#
s=$[$(date +%s) - ${T0}]
printf "Elapsed time: %ds (%s)\n" "$s" $(${cwd}/python/nicetime.py $s)
MessageBreak "DONE $(basename $0)"
#

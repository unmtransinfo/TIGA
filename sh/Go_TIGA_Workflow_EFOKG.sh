#!/bin/bash
#############################################################################
### Go_TIGA_Workflow_EFOKG.sh - Steps for EFO relationships.
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
#
###
efofile="${ODIR}/efo.tsv"
tsvfile_trait_sub="${ODIR}/efo_sub_gwas.tsv.gz" # Large, thus gzip.
###
#############################################################################
###
${cwd}/R/gwascat_trait_efokg.R
#
s=$[$(date +%s) - ${T0}]
printf "Elapsed time: %ds (%s)\n" "$s" $(${cwd}/python/nicetime.py $s)
MessageBreak "DONE $(basename $0)"
#

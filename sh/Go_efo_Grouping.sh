#!/bin/bash
###
# See https://github.com/IUIDSL/(iu_idsl_util,iu_idsl_jena)
###
cwd=$(pwd)
#
DATADIR=${cwd}/data
#
EFO_DIR="$(cd $HOME/../data; pwd)/EFO/data"
OWLFILE="$EFO_DIR/efo.owl"
#
LIBDIR="$(cd $HOME/../app; pwd)/lib"
#
# GWASCatalog release:
if [ $# -eq 1 ]; then
	GC_REL=$1
else
	printf "ERROR: syntax $(basename $0) \"YYYY-MM-DD\"\n"
	exit
fi
#
GC_REL_Y=$(echo $GC_REL |sed 's/-.*$//')
GC_REL_M=$(echo $GC_REL |sed 's/^.*-\(.*\)-.*$/\1/')
GC_REL_D=$(echo $GC_REL |sed 's/^.*-//')
#
ODIR="${DATADIR}/${GC_REL_Y}${GC_REL_M}${GC_REL_D}"
printf "Output dir: ${ODIR}\n"
#
###
EFO_URL="https://github.com/EBISPOT/efo/releases/download/v${EFO_RELEASE}/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
#############################################################################
EFO_RELEASE="3.40.0"
printf "${EFO_RELEASE}\n" >${ODIR}/efo_release.txt
#
#############################################################################
#
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${ODIR}/efo.tsv
#
# From efo.tsv create GraphML file:
graphmlfile="${ODIR}/efo_graph.graphml"
${cwd}/R/efo_graph.R $GC_REL_Y $GC_REL_M $GC_REL_D
gzip -f ${graphmlfile}
#
#############################################################################
# Exports:
###
# CYJS: Not needed now but maybe later?
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2cyjs -o ${ODIR}/efo.cyjs
###
# Edgelist and nodelist
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2edgelist -o ${ODIR}/efo_edgelist.tsv
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2nodelist -o ${ODIR}/efo_nodelist.tsv
#
#############################################################################
###
# Grouping:
cat $ODIR/gwascat_trait.tsv \
	|sed -e '1d' |awk -F '\t' '{print $5}' |sort -u \
	>$ODIR/gwascatalog.efoid
printf "EFO IDs: %d\n" $(cat $ODIR/gwascatalog.efoid |grep '^EFO' |wc -l)
###
${cwd}/python/nx_analysis.py cluster --min_groupsize 2 --max_level 10 \
	--i_edge $ODIR/efo_edgelist.tsv \
	--i_node_attr $ODIR/efo_nodelist.tsv \
	--i_node_set $ODIR/gwascatalog.efoid \
	--setname gwc --graphname "EFO" \
	--o $ODIR/efo_groups.tsv
#
#
${cwd}/R/efo_groups.R \
	"mood disorder" \
	$ODIR/efo.tsv \
	$ODIR/efo_groups.tsv \
	$ODIR/efo_subgraph.graphml \
###
# To do:
#   * Groups file usable by TIGA app for marker color and selection.
#   * Associate child nodes with selected ancestor nodes, defining groups.
###

#!/bin/bash
#
cwd=$(pwd)
#
DATADIR=${cwd}/data
#
if [ "$(uname -s)" = "Darwin" ]; then
	EFO_DIR="/Users/data/EFO/data"
else
	EFO_DIR="/home/data/EFO/data"
fi
#
OWLFILE="$EFO_DIR/efo.owl"
EFO_URL="https://github.com/EBISPOT/efo/releases/download/current/efo.owl"
wget -q -O $OWLFILE $EFO_URL
#
if [ "$(uname -s)" = "Darwin" ]; then
	LIBDIR="/Users/app/lib"
else
	LIBDIR="/home/app/lib"
fi
#
###
# See https://github.com/IUIDSL/iu_idsl_util and
# https://github.com/IUIDSL/iu_idsl_jena
java -classpath $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar edu.indiana.sice.idsl.jena.jena_utils \
	-ontfile ${OWLFILE} -vv -ont2tsv -o ${DATADIR}/efo.tsv
#
###
# CYJS: Not needed now but maybe later?
java -classpath $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar edu.indiana.sice.idsl.jena.jena_utils \
	-ontfile ${OWLFILE} -vv -ont2cyjs -o ${DATADIR}/efo.cyjs
java -classpath $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar edu.indiana.sice.idsl.jena.jena_utils \
	-ontfile ${OWLFILE} -vv -ont2edgelist -o ${DATADIR}/efo_edgelist.tsv
java -classpath $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar edu.indiana.sice.idsl.jena.jena_utils \
	-ontfile ${OWLFILE} -vv -ont2nodelist -o ${DATADIR}/efo_nodelist.tsv
#
###
#
cat $DATADIR/gwascat_trait.tsv \
	|awk -F '\t' '{print $4}' \
	|sed -e '1d' \
	|sort -u \
	>$DATADIR/gwascatalog.efoid
###
${cwd}/python/nx_analysis.py cluster --min_groupsize 2 --max_level 10 \
	--i_edge $DATADIR/efo_edgelist.tsv \
	--i_node_attr $DATADIR/efo_nodelist.tsv \
	--i_node_set $DATADIR/gwascatalog.efoid \
	--setname gwc --graphname "EFO" \
	--o $DATADIR/efo_groups.tsv
#
###
# To do:
#   * Groups file usable by TIGA app for marker color and selection.
#   * Associate child nodes with selected ancestor nodes, defining groups.
###

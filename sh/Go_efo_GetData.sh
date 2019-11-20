#!/bin/bash
#
EFODIR="/home/data/EFO/data"
#
OWLFILE="$EFODIR/efo.owl"
#
cwd=$(pwd)
#
DATADIR=${cwd}/data
#
cd $HOME/src/iu_idsl_jena
#
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2tsv -o ${DATADIR}/efo.tsv"
#
###
# CYJS: Not needed now but maybe later?
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2cyjs -o ${DATADIR}/efo.cyjs"
#
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2edgelist -o ${DATADIR}/efo_edgelist.tsv"
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2nodelist -o ${DATADIR}/efo_nodelist.tsv"
#
cd ${cwd}
${cwd}/python/nx_analysis.py cluster \
	--i_edge $DATADIR/efo_edgelist.tsv \
	--i_node_attr $DATADIR/efo_nodelist.tsv \
	--i_node_set $DATADIR/gwascatalog.efoid \
	--setname gwc --graphname "EFO" \
	--o $DATADIR/efo_groups.tsv
#
###
# To do: groups file usable by GWAX app for marker color and selection.

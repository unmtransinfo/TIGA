#!/bin/bash
#
EFODIR="/home/data/EFO/data"
#
OWLFILE="$EFODIR/efo.owl"
#
cwd=$(pwd)
#
cd $HOME/src/iu_idsl_jena
#
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2tsv -o ${cwd}/data/efo.tsv"
#
#mvn exec:java \
#	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
#	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2cyjs -o ${cwd}/data/efo.cyjs"
#
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2edgelist -o ${cwd}/data/efo_edgelist.tsv"
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${OWLFILE} -vv -ont2nodelist -o ${cwd}/data/efo_nodelist.tsv"
#

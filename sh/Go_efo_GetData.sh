#!/bin/bash
#
EFODIR="/home/data/EFO/data"
cwd=$(pwd)
#
cd $HOME/src/iu_idsl_jena
#
mvn exec:java \
	-Dexec.mainClass="edu.indiana.sice.idsl.jena.jena_utils" \
	-Dexec.args="-ontfile ${EFODIR}/efo.owl -vv -ont2tsv -o ${cwd}/data/efo.tsv"
#	-Dexec.args="-ontfile ${EFODIR}/efo.owl -vv -ont2cyjs -o ${cwd}/data/efo.cyjs"
#

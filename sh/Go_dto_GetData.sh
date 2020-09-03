#!/bin/bash
###
# See https://github.com/IUIDSL/(iu_idsl_util,iu_idsl_jena)
###
cwd=$(pwd)
#
DATADIR=${cwd}/data
#
DTO_DIR="$HOME/../data/DTO/data"
OWLFILE="$DTO_DIR/dto_complete_merged.owl"
DTO_URL="https://github.com/DrugTargetOntology/DTO/releases/download/V1.1/dto_complete_merged.owl"
#
#
wget -O $OWLFILE $DTO_URL
#
LIBDIR="$HOME/../app/lib"
#
dtofile="${DATADIR}/dto.tsv"
###
# TSV
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar \
	-ifile_ont ${OWLFILE} -vv -ont2tsv -o ${dtofile}
###

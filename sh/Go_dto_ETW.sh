#!/bin/bash
###
# ETW = Extract, Transform & Wrangle
# See https://github.com/IUIDSL/(iu_idsl_util,iu_idsl_jena)
###
cwd=$(pwd)
#
DATADIR=${cwd}/data
#
# DTO v1.1 latest since 2017
DTO_URL="https://github.com/DrugTargetOntology/DTO/releases/download/V1.1/dto_complete_merged.owl"
#
DTO_DIR=$(cd $HOME/../data/DTO/data; pwd)
DTO_OWLFILE="$DTO_DIR/dto_complete_merged.owl"
#
wget -O $DTO_OWLFILE $DTO_URL
#
LIBDIR=$(cd $HOME/../app/lib; pwd)
#
dto_tsvfile="${DATADIR}/dto.tsv"
dto_tlscm_file="${DATADIR}/dto_complete_merged_toplevelsuperclassmembership.tsv"
#
set -x
###
# TSV
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar -v \
	-ifile_ont ${DTO_OWLFILE} -ont2tsv -o ${dto_tsvfile}
###
# TLSCM
java -jar $LIBDIR/iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar -v \
	-ifile_ont ${DTO_OWLFILE} \
	-list_toplevelsuperclassmembership -o ${dto_tlscm_file}
###

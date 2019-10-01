#!/bin/bash
#
DBHOST="tcrd.kmc.io"
DBNAME="tcrd600"
DBUSR="tcrd"
DBPW=""
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
#
${cwd}/sh/runsql_my.sh -c \
	-h $DBHOST -n $DBNAME -u $DBUSR -p "$DBPW" \
	-f ${cwd}/sql/tcrd_targets.sql \
	>${DATADIR}/tcrd_targets.tsv
#

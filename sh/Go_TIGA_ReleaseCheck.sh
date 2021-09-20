#!/bin/bash
#############################################################################
### Go_TIGA_ReleaseCheck.sh - Check for new and recent releases
### auto-downloaded by LFTP_GWASCatalog.sh.
#############################################################################
#
#
MSPAN=3 #How many months back?
#
GWASCATALOGDIR="$(cd $HOME/../data/GWASCatalog; pwd)"
###
Y="$(date +'%Y')"
M="$(date +'%m')"
###
function YearMonthCheck {
	y="$1"
	m="$2"
	m=$(printf "%02d" $m)
	if [ -e  "${GWASCATALOGDIR}/releases/${y}/${m}" ]; then
		for d in $(ls ${GWASCATALOGDIR}/releases/${y}/${m}) ; do
			dir="${GWASCATALOGDIR}/releases/${y}/${m}/${d}"
			fcount="$(ls -1 $dir |wc -l)"
			printf "\t${y}-${m}-${d}\t${dir} (${fcount} files)\n"
		done
	fi
}
###
if [ $M -lt $MSPAN ]; then
	y="$(($Y - 1))"
	for m in $(seq $((12 + 10#$M - $MSPAN)) 12) ; do
		YearMonthCheck $y $m
	done
	y="$Y"
	for m in $(seq 1 $((10#$M))) ; do
		YearMonthCheck $y $m
	done
else
	y="$Y"
	for m in $(seq $((10#$M - $MSPAN)) $M) ; do
		YearMonthCheck $y $m
	done
fi
#

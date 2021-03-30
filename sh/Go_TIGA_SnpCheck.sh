#!/bin/bash
#
cwd=$(pwd)
DATADIR=${cwd}/data
#
${pwd}/python/tiga_gt_snpcheck.py -v -g ENSG00000160785 -t EFO_0004541 \
	--o $DATADIR/tiga_gt_snpcheck_ENSG00000160785-EFO_0004541.tsv
#

#!/bin/bash

#
cwd=$(pwd)
DATADIR="${cwd}/data"

RELEASEA="20200715"
RELEASEB="20210212"
#
${cwd}/python/tiga_gt_compare.py gene \
	--iA data/$RELEASEA/gt_stats.tsv.gz \
	--iB data/$RELEASEB/gt_stats.tsv.gz \
	--explain_missing --filterfile data/$RELEASEB/filtered_genes.tsv
#
${cwd}/python/tiga_gt_compare.py trait \
	--iA data/$RELEASEA/gt_stats.tsv.gz \
	--iB data/$RELEASEB/gt_stats.tsv.gz \
	--explain_missing --filterfile data/$RELEASEB/filtered_traits.tsv
#

#!/bin/bash
###
# https://vegas2.qimrberghofer.edu.au/
# https://www.cog-genomics.org/plink2
###
# Input: TSV|SSV, no-header, 2-cols: SNP rs-number and p-value. NA values not allowed.
###
# Output (SSV):
# Columns: Chromosome, Gene, Number of SNPs, Number of simulations, Start position, Stop position, Gene-based test statistic, P-value, Top-SNP, Top-SNP pvalue.
###
# VEGAS2 gene based scores can be computed for each GWAS study, from the
# SNP+pval files provided by the GWAS Catalog.
###
# ~250GB RAM required.
# Parallelize?
###
#
cwd=$(pwd)
#
DATADIR="${cwd}/data/vegas"
#
VEGAS="${cwd}/perl/vegas2v2.pl"
#
cd $DATADIR
#
###
# Example file has 26 genes; full file ~26k.
#GLISTFILE="/home/data/VEGAS/data/VEGAS2v2example/example.glist"
GLISTFILE="/home/data/VEGAS/data/glist-hg19.txt"
#
###
# CUSTOM_DIR should contain: ${CUSTOM_PREFIX}.(bed,bim,fam)
###
CUSTOM_DIR="/home/data/VEGAS/data/g1000p3"
CUSTOM_PREFIX="g1000p3_EUR"
#
OUT_PREFIX="vegas_example"
#
date
T0=$(date +%s)
#
n_file=$(ls -1 vegas_in_*.tsv |wc -l)
i=0
#
for ifile in $(ls vegas_in_*.tsv) ; do
	i=$(($i+1))
	acc=$(echo $ifile |perl -pe 's/vegas_in_(.*)\.tsv/$1/')
	printf "%d/%d. %s\n" "${i}" "${n_file}" "${acc}"
	out_prefix="vegas_out_${acc}"
	T0_this=$(date +%s)
	$VEGAS -G \
		-glist $GLISTFILE \
		-custom ${CUSTOM_DIR}/${CUSTOM_PREFIX} \
		-snpandp $ifile \
		-out $out_prefix
	printf "%s: elapsed time: %ds\n" "${acc}" "$[$(date +%s) - ${T0_this}]"
done
#
printf "TOTAL elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
date
#


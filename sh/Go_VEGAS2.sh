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
###
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
#
VEGAS="${cwd}/perl/vegas2v2.pl"
#
cd $DATADIR
#
###
INFILE="/home/data/VEGAS/data/VEGAS2v2example/example.txt"
GLISTFILE="/home/data/VEGAS/data/VEGAS2v2example/example.glist"
#
# CUSTOM_DIR should contain: ${CUSTOM_PREFIX}.(bed,bim,fam)
#CUSTOM_DIR="/home/data/VEGAS/data/VEGAS2v2example"
#CUSTOM_PREFIX="example"
CUSTOM_DIR="/home/data/VEGAS/data/g1000p3"
CUSTOM_PREFIX="g1000p3_EUR"
#
OUT_PREFIX="vegas_example"
#
$VEGAS -G \
	-snpandp $INFILE \
	-glist $GLISTFILE \
	-custom ${CUSTOM_DIR}/${CUSTOM_PREFIX} \
	-out $OUT_PREFIX
#

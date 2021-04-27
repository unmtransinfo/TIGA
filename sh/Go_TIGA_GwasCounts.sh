#!/bin/bash
###

cwd=$(pwd)
DATADIR="${cwd}/data"
#
ODIR="${DATADIR}/20210329"
#
#	-v -v \
${cwd}/python/tiga_gwas_counts.py \
	--ifile_gwas $ODIR/gwascat_gwas.tsv \
	--ifile_assn $ODIR/gwascat_assn.tsv \
	--ifile_trait $ODIR/gwascat_trait.tsv \
	--ifile_snp2gene $ODIR/gwascat_snp2gene_MERGED.tsv \
	--ifile_icite $ODIR/gwascat_icite.tsv \
	--ofile_gwas $ODIR/gwascat_gwas_counts.tsv \
	--ofile_trait $ODIR/gwascat_trait_counts.tsv
#

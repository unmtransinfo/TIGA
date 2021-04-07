#!/bin/bash
#
cwd=$(pwd)
#
DATADIR="${cwd}/data/20210212"
#
# Generate counts table as in paper.
# Input: TCRD file from BioClients.idg.tcrd.Client listTargets.
# Input: TCRD2DTO file from BioClients.idg.tcrd.Client listTargetsByDTO.
# Input: Gene-traits stats file from tiga_gt_stats.R.
#
#${cwd}/R/tiga_idg.R \
#	${DATADIR}/tcrd_targets.tsv \
#	${DATADIR}/tcrd2dto.tsv \
#	${DATADIR}/gt_stats.tsv.gz \
#	${DATADIR}/tdl_fam_counts_TIGA.tsv \
#	${DATADIR}/tdl_fam_counts_MERGED.tsv \
#	${DATADIR}/tdl_dto_counts_TIGA.tsv \
#	${DATADIR}/tdl_dto_counts_MERGED.tsv
##
# TCRDv6.4_DTO.tsv from TCRDv6.4_DTO.xlsx, from TIO. 
${cwd}/R/tiga_idg_dtoL2.R \
	${DATADIR}/tcrd_targets.tsv \
	${DATADIR}/TCRDv6.4_DTO.tsv \
	${DATADIR}/gt_stats.tsv.gz \
	${DATADIR}/tdl_dto_counts_TIGA.tsv \
	${DATADIR}/tdl_dto_counts_MERGED.tsv
#

#!/bin/bash
#############################################################################
### Go_gwascat_NetWAS.sh
### Generate file for input to VEGAS+NetWAS, for a specified disease/trait.
### https://vegas2.qimrberghofer.edu.au/
### https://hb.flatironinstitute.org/netwas
### Format: <SNP>\t<PVAL>
### EFO_0002508: "Parkinson's disease"
### http://www.ebi.ac.uk/efo/EFO_0002508
#############################################################################
#
cwd=$(pwd)
#
DBNAME="gwascatalog"
#
DATADIR="${cwd}/data"
ofile="${DATADIR}/pd_snps.tsv"
#
sql="\
SELECT snps, p_value FROM assn \
WHERE mapped_trait_uri = 'http://www.ebi.ac.uk/efo/EFO_0002508' \
"
###
mysql -ABr --execute="${sql}" ${DBNAME} \
	|sed -e '1d' \
	>${ofile}
###
printf "Output file: %s\n" "${ofile}"
printf "SNPs: %d\n" "$(cat ${ofile} |wc -l)"
#

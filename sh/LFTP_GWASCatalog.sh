#!/bin/bash
#
YEAR="$(date +'%Y')"
#
lftp ftp://anonymous:@ftp.ebi.ac.uk -e "mirror /pub/databases/gwas/releases/${YEAR} /home/data/GWASCatalog/releases/${YEAR} ; quit"
#

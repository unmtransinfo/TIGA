#!/bin/bash
#############################################################################
### Go_TIGA_DbCreate.sh (MySql) 
#############################################################################
### OBSOLETE -- NOW IMPLEMENTED ALL IN Python/Pandas.
#############################################################################
### EFO-subclass study-study pairs efo_sub_gwas.tsv now from gwascat_trait.R.
#############################################################################
### gwas_counts table:
###   study_accession - study ID
###   snp_count - per study
###   trait_count - per study
###   assn_count - per study
###   gene_r_count - reported genes per study
###   gene_m_count - mapped genes per study
###   study_perpmid_count - per PMID (the PMID for this study)
#############################################################################
#
cwd=$(pwd)
#
printf "Started: %s\n" "$(date)"
#
if [ $# -gt 0 ]; then
	DATADIR=$1
else
	DATADIR="${cwd}/data"
fi
printf "Input DATADIR: %s\n" "${DATADIR}"
#
if [ $# -gt 1 ]; then
	DBNAME=$2
else
	DBNAME="tiga"
fi
printf "DBNAME: %s\n" "${DBNAME}"
#
# INPUT FILES:
if [ $# -gt 2 ]; then
	gwasfile=$3
else
	gwasfile="${DATADIR}/gwascat_gwas.tsv" #gwascat_gwas.R
fi
if [ $# -gt 3 ]; then
	assnfile=$4
else
	assnfile="${DATADIR}/gwascat_assn.tsv" #gwascat_assn.R
fi
if [ $# -gt 4 ]; then
	snp2genefile=$5
else
	snp2genefile="${DATADIR}/gwascat_snp2gene.tsv" #Go_TIGA_Workflow.sh,snp2gene_reported.pl,snp2gene_mapped.pl
fi
if [ $# -gt 5 ]; then
	traitfile=$6
else
	traitfile="${DATADIR}/gwascat_trait.tsv" #gwascat_trait.R
fi
if [ $# -gt 6 ]; then
	icitefile=$7
else
	icitefile="${DATADIR}/gwascat_icite.tsv" #BioClients.icite.Client
fi
# OUTPUT FILES:
if [ $# -gt 7 ]; then
	ofile_gwas=$8
else
	ofile_gwas="${DATADIR}/gwascat_counts.tsv"
fi
if [ $# -gt 8 ]; then
	ofile_trait=$9
else
	ofile_trait="${DATADIR}/trait_counts.tsv"
fi
#
printf "Input gwasfile: %s\n" "$gwasfile"
printf "Input assnfile: %s\n" "$assnfile"
printf "Input snp2genefile: %s\n" "$snp2genefile"
printf "Input traitfile: %s\n" "$traitfile"
printf "Input icitefile: %s\n" "$icitefile"
printf "Output gwas counts: %s\n" "$ofile_gwas"
printf "Output trait counts: %s\n" "$ofile_trait"
#
error=""
for f in $gwasfile $assnfile $snp2genefile $traitfile $icitefile ; do
	if [ ! -f "${f}" ]; then
		printf "ERROR: Input file not found: %s\n" "${f}"
		error="TRUE"
	fi
done
if [ "$error" ]; then
	printf "Quitting due to missing file[s].\n"
	exit
fi
#
mysql -v -e "DROP DATABASE $DBNAME"
mysql -v -e "CREATE DATABASE $DBNAME"
#
mysql -D $DBNAME <${cwd}/sql/tiga_create_tables.sql
###
mysql -D $DBNAME -e "SET GLOBAL local_infile = true"
#
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${gwasfile}' INTO TABLE gwas FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${assnfile}' INTO TABLE assn FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${snp2genefile}' INTO TABLE snp2gene FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${traitfile}' INTO TABLE trait2study FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${icitefile}' INTO TABLE icite FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
#
###
mysql -D $DBNAME -e "UPDATE trait2study SET mapped_trait = NULL WHERE mapped_trait IN ('', 'NA')"
mysql -D $DBNAME -e "UPDATE trait2study SET mapped_trait_uri = NULL WHERE mapped_trait_uri IN ('', 'NA')"
###
mysql -D $DBNAME -e "ALTER TABLE gwas COMMENT = 'GWAS Catalog studies (from raw file, all cols)'"
mysql -D $DBNAME -e "ALTER TABLE assn COMMENT = 'GWAS Catalog associations (OR and beta separated)'"
mysql -D $DBNAME -e "ALTER TABLE snp2gene COMMENT = 'GWAS Catalog  gene associations'"
mysql -D $DBNAME -e "ALTER TABLE icite COMMENT = 'GWAS Catalog iCite publication annotations'"
mysql -D $DBNAME -e "ALTER TABLE trait2study COMMENT = 'GWAS Catalog study traits (EFO, GO, HP) mapped to study accession IDs'"
###
mysql -D $DBNAME -e "UPDATE assn SET upstream_gene_id = NULL WHERE upstream_gene_id = ''"
mysql -D $DBNAME -e "UPDATE assn SET downstream_gene_id = NULL WHERE downstream_gene_id = ''"
mysql -D $DBNAME -e "UPDATE assn SET snp_gene_ids = NULL WHERE snp_gene_ids = ''"
###
mysql -D $DBNAME -e "CREATE INDEX g_study_accession_idx ON gwas (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX ga_study_accession_idx ON assn (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_study_accession_idx ON snp2gene (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_snp_idx ON snp2gene (snp)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_gsymb_idx ON snp2gene (gsymb)"
#
###
printf "Creating gwas_counts table.\n"
# Tables required: gwas, snp2gene, trait2study, assn.
mysql -D $DBNAME <${cwd}/sql/create_gwas_counts_table.sql
#
printf "Writing gwas_counts table to: %s\n" "${ofile_gwas}"
mysql -D $DBNAME -ABr --execute="SELECT * FROM gwas_counts" >${ofile_gwas}
###
printf "Writing trait_counts table to: %s\n" "${ofile_trait}"
(mysql -D $DBNAME <<__EOF__
SELECT
	mapped_trait_uri,
	mapped_trait,
	COUNT(DISTINCT study_accession) AS "n_study"
FROM
	trait2study
GROUP BY
	mapped_trait_uri,
	mapped_trait
ORDER BY
	n_study DESC
	;
__EOF__
) >${ofile_trait}
#
#
printf "Done: %s\n" "$(date)"
#

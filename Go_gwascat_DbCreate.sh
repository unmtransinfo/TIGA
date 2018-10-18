#!/bin/bash
#############################################################################
### Go_gwascat_DbCreate.sh - MySql version
#############################################################################
### To drop db:
### mysql -v -e "DROP DATABASE $DBNAME"
#############################################################################
### gwas_counts table:
###   study_accession - study ID
###   snp_count - per study
###   trait_count - per study
###   assn_count - per study
###   gene_r_count - reported genes per study
###   gene_m_count - mapped genes per study
###   study_perpmid_count - per PMID (the PMID for this study)
### Question: Can we compute a fractional citation count to a gene-trait association?
#############################################################################
### Jeremy Yang
#############################################################################
#
set -x
#
cwd=$(pwd)
#
printf "Started: %s\n" "$(date)"
#
DBNAME="gwascatalog"
#
DATADIR="${cwd}/data"
#
mysql -v -e "CREATE DATABASE $DBNAME"
#
mysql -D $DBNAME <${cwd}/sql/gwascatalog_create_tables.sql
###
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gwascat_gwas.tsv' INTO TABLE gwas FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gwascat_assn.tsv' INTO TABLE assn FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gwascat_snp2gene.tsv' INTO TABLE snp2gene FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gwascat_trait.tsv' INTO TABLE trait FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gwascat_icite.tsv' INTO TABLE icite FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/gt_stats.tsv' INTO TABLE gt_stats FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
###
mysql -D $DBNAME -e "ALTER TABLE gwas COMMENT = 'GWAS Catalog studies (from raw file, all cols)'"
mysql -D $DBNAME -e "ALTER TABLE assn COMMENT = 'GWAS Catalog associations (OR and beta separated)'"
mysql -D $DBNAME -e "ALTER TABLE snp2gene COMMENT = 'GWAS Catalog  gene associations'"
mysql -D $DBNAME -e "ALTER TABLE icite COMMENT = 'GWAS Catalog iCite pub annotations'"
mysql -D $DBNAME -e "ALTER TABLE trait COMMENT = 'GWAS Catalog study traits (EFO, GO, HP)'"
###
mysql -D $DBNAME -e "UPDATE assn SET upstream_gene_id = NULL WHERE upstream_gene_id = ''"
mysql -D $DBNAME -e "UPDATE assn SET downstream_gene_id = NULL WHERE downstream_gene_id = ''"
mysql -D $DBNAME -e "UPDATE assn SET snp_gene_ids = NULL WHERE snp_gene_ids = ''"
###
mysql -D $DBNAME -e "CREATE INDEX g_study_accession_idx ON gwas (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX ga_study_accession_idx ON assn (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_study_accession_idx ON snp2gene (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_snp_idx ON snp2gene (snp (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_gsymb_idx ON snp2gene (gsymb)"
#
#
mysql -D $DBNAME -e "ALTER TABLE gt_stats COMMENT = 'GWAS gene-trait stats, used by GWAX web app'"
mysql -D $DBNAME -e "UPDATE gt_stats SET name = NULL WHERE name = 'NA'"
mysql -D $DBNAME -e "UPDATE gt_stats SET fam = NULL WHERE fam = 'NA'"
mysql -D $DBNAME -e "UPDATE gt_stats SET tdl = NULL WHERE tdl = 'NA'"
#
###
printf "Creating gwas_counts table:\n"
#
mysql -D $DBNAME <${cwd}/sql/create_gwas_counts_table.sql
#
printf "Done: %s\n" "$(date)"
#

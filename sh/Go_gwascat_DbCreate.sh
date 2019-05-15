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
mysql -D $DBNAME -e "UPDATE trait SET mapped_trait = NULL WHERE mapped_trait IN ('', 'NA')"
mysql -D $DBNAME -e "UPDATE trait SET mapped_trait_uri = NULL WHERE mapped_trait_uri IN ('', 'NA')"
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
mysql -D $DBNAME -e "CREATE INDEX g_study_accession_idx ON gwas (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX ga_study_accession_idx ON assn (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_study_accession_idx ON snp2gene (study_accession)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_snp_idx ON snp2gene (snp)"
mysql -D $DBNAME -e "CREATE INDEX gs2g_gsymb_idx ON snp2gene (gsymb)"
#
###
#EFO
mysql -D $DBNAME <<__EOF__
CREATE TABLE efo (
	node_or_edge VARCHAR(8),
	id VARCHAR(16),
	label VARCHAR(16),
	comment VARCHAR(128),
	source VARCHAR(64),
	target VARCHAR(64),
	uri VARCHAR(64)
);
__EOF__
###
mysql -D $DBNAME -e "LOAD DATA LOCAL INFILE '${DATADIR}/efo.tsv' INTO TABLE efo FIELDS TERMINATED BY '\t' IGNORE 1 LINES;"
#
mysql -D $DBNAME -e "CREATE TABLE efo_classes AS SELECT id,label,comment,uri FROM efo WHERE node_or_edge='node'"
mysql -D $DBNAME -e "CREATE TABLE efo_sub AS SELECT source AS trait_uri, target AS subclass_uri FROM efo WHERE node_or_edge='edge'"
mysql -D $DBNAME -e "DROP TABLE efo"
mysql -D $DBNAME -e "CREATE INDEX efosub_t_idx ON efo_sub (trait_uri)"
mysql -D $DBNAME -e "CREATE INDEX efosub_s_idx ON efo_sub (subclass_uri)"
# Link studies by EFO subclass relationships.
mysql -D $DBNAME <<__EOF__
CREATE TABLE efo_sub_gwas AS
SELECT
	t1.study_accession,
	t1.mapped_trait AS trait,
	efo_sub.trait_uri,
	t2.study_accession AS study_accession_subclass,
	t2.mapped_trait AS subclass_trait,
	efo_sub.subclass_uri
FROM
	trait t1,
	trait t2,
	efo_sub
WHERE
	t1.mapped_trait_uri = efo_sub.trait_uri
	AND t2.mapped_trait_uri = efo_sub.subclass_uri
__EOF__
###
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
# EFO-subclass-based GWAS study-study associations:
runsql_my.sh -q 'SELECT * FROM efo_sub_gwas' -c >${DATADIR}/efo_sub_gwas.tsv
#
printf "Done: %s\n" "$(date)"
#

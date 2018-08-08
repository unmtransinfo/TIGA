#!/bin/sh
#############################################################################
### Go_gwascat_DbCreate.sh - MySql version
### 
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
### 13 Nov 2017
#############################################################################
#
# set -e
set -x
#
DBNAME="gwascatalog"
#
DATADIR="data"
#
mysql -v -e "CREATE DATABASE $DBNAME"
#
#
csvfiles="\
${DATADIR}/gwascat_gwas.tsv \
${DATADIR}/gwascat_assn.csv \
${DATADIR}/gwascat_snp2gene.tsv \
${DATADIR}/gwascat_trait.tsv \
${DATADIR}/gwascat_icite.csv \
${DATADIR}/gt_stats.csv"
#
for csvfile in $csvfiles ; do
	#
	if [ $(echo $csvfile |sed -e 's/^.*\.//') = "tsv" ]; then
		opts="--tsv"
	else
		opts=""
	fi
	#
	opts="$opts --maxchar 4096"
	#
	if [ $(echo "$csvfile" | grep '^.*gwascat_[a-z0-9]*\.[ct]sv') ]; then
		tname=$(echo "$csvfile" |sed -e 's/^.*gwascat_\(.*\)\.[ct]sv/\1/')
		opts="$opts --tablename ${tname}"
	fi
	#
	csv2sql.py $opts \
		--dbsystem "mysql" \
		--i $csvfile \
		--create \
		--fixtags \
		|mysql -D $DBNAME
	#
	csv2sql.py $opts \
		--dbsystem "mysql" \
		--i $csvfile \
		--insert \
		--fixtags \
		|mysql -D $DBNAME
	#
done
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
mysql -D $DBNAME -e "UPDATE assn SET upstream_gene_distance = NULL WHERE upstream_gene_distance = ''"
mysql -D $DBNAME -e "UPDATE assn SET downstream_gene_distance = NULL WHERE downstream_gene_distance = ''"
mysql -D $DBNAME -e "ALTER TABLE assn MODIFY COLUMN oddsratio FLOAT"
mysql -D $DBNAME -e "ALTER TABLE assn MODIFY COLUMN beta FLOAT"
###
mysql -D $DBNAME -e "CREATE INDEX g_study_accession_idx ON gwas (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX ga_study_accession_idx ON assn (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_study_accession_idx ON snp2gene (study_accession (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_snp_idx ON snp2gene (snp (32))"
mysql -D $DBNAME -e "CREATE INDEX gs2g_gsymb_idx ON snp2gene (gsymb (32))"
#
mysql -D $DBNAME -e "ALTER TABLE gwas MODIFY COLUMN date DATE"
mysql -D $DBNAME -e "ALTER TABLE gwas MODIFY COLUMN date_added_to_catalog DATE"
mysql -D $DBNAME -e "ALTER TABLE gwas MODIFY COLUMN association_count INTEGER"
#
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN n_traits_g INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN n_genes_t INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN n_snp INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN n_study INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN pvalue_mlog_median FLOAT"
mysql -D $DBNAME -e "ALTER TABLE gt_stats MODIFY COLUMN or_median FLOAT"
mysql -D $DBNAME -e "UPDATE gt_stats SET pvalue_mlog_median = NULL WHERE pvalue_mlog_median = 'Inf'"
mysql -D $DBNAME -e "UPDATE gt_stats SET name = NULL WHERE name = 'NA'"
mysql -D $DBNAME -e "UPDATE gt_stats SET fam = NULL WHERE fam = 'NA'"
mysql -D $DBNAME -e "UPDATE gt_stats SET tdl = NULL WHERE tdl = 'NA'"
#
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN citation_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN citations_per_year FLOAT"
mysql -D $DBNAME -e "UPDATE icite SET nih_percentile = NULL WHERE nih_percentile = ''"
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN nih_percentile FLOAT"
mysql -D $DBNAME -e "UPDATE icite SET field_citation_rate = NULL WHERE field_citation_rate = ''"
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN field_citation_rate FLOAT"
mysql -D $DBNAME -e "UPDATE icite SET expected_citations_per_year = NULL WHERE expected_citations_per_year = ''"
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN expected_citations_per_year FLOAT"
mysql -D $DBNAME -e "UPDATE icite SET relative_citation_ratio = NULL WHERE relative_citation_ratio = ''"
mysql -D $DBNAME -e "ALTER TABLE icite MODIFY COLUMN relative_citation_ratio FLOAT"
#
###
#
mysql -D $DBNAME <<__EOF__
CREATE TABLE
	gwas_counts
AS SELECT
	g.study_accession,
	t1.snp_count,
	t2.trait_count,
	t3.assn_count,
	t4.gene_count AS "gene_r_count",
	t5.gene_count AS "gene_m_count",
	t6.study_perpmid_count
FROM
	gwas g
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.snp) AS "snp_count",
		gs2g.study_accession
	FROM
		snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped = 'r'
	GROUP BY
		gs2g.study_accession
	) t1 ON t1.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gt.trait_uri) AS "trait_count",
		gt.study_accession
	FROM
		trait gt
	GROUP BY
		gt.study_accession
	) t2 ON t2.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(ga.p_value) AS "assn_count",
		ga.study_accession
	FROM
		assn ga
	GROUP BY
		ga.study_accession
	) t3 ON t3.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.gsymb) AS "gene_count",
		gs2g.study_accession
	FROM
		snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped = 'r'
	GROUP BY
		gs2g.study_accession
	) t4 ON t4.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.gsymb) AS "gene_count",
		gs2g.study_accession
	FROM
		snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped IN ('m', 'mu', 'md')
	GROUP BY
		gs2g.study_accession
	) t5 ON t5.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		pubmedid,
		COUNT(DISTINCT study_accession) AS "study_perpmid_count"
	FROM    
		gwas g
	GROUP BY
		pubmedid
	) t6 ON t6.pubmedid = g.pubmedid
ORDER BY
	t3.assn_count DESC
	;
__EOF__
#
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN snp_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN trait_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN assn_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN gene_r_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN gene_m_count INTEGER"
mysql -D $DBNAME -e "ALTER TABLE gwas_counts MODIFY COLUMN study_perpmid_count INTEGER"
#
mysql -D $DBNAME -e "UPDATE gwas_counts SET assn_count = 0 WHERE assn_count IS NULL"
mysql -D $DBNAME -e "UPDATE gwas_counts SET snp_count = 0 WHERE snp_count IS NULL"
mysql -D $DBNAME -e "UPDATE gwas_counts SET gene_r_count = 0 WHERE gene_r_count IS NULL"
mysql -D $DBNAME -e "UPDATE gwas_counts SET gene_m_count = 0 WHERE gene_m_count IS NULL"
#

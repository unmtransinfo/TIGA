--
CREATE TABLE
	gwas_counts
AS SELECT
	gwas.study_accession,
	t1.snp_count,
	t2.trait_count,
	t3.assn_count,
	t4.gene_count AS "gene_r_count",
	t5.gene_count AS "gene_m_count",
	t6.study_perpmid_count
FROM
	gwas
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT snp) AS "snp_count",
		study_accession
	FROM
		snp2gene
	WHERE
		reported_or_mapped = 'r'
	GROUP BY
		study_accession
	) t1 ON t1.study_accession = gwas.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT mapped_trait_uri) AS "trait_count",
		study_accession
	FROM
		trait
	GROUP BY
		study_accession
	) t2 ON t2.study_accession = gwas.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(p_value) AS "assn_count",
		study_accession
	FROM
		assn
	GROUP BY
		study_accession
	) t3 ON t3.study_accession = gwas.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gsymb) AS "gene_count",
		study_accession
	FROM
		snp2gene
	WHERE
		reported_or_mapped = 'r'
	GROUP BY
		study_accession
	) t4 ON t4.study_accession = gwas.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gsymb) AS "gene_count",
		study_accession
	FROM
		snp2gene
	WHERE
		reported_or_mapped IN ('m', 'mu', 'md')
	GROUP BY
		study_accession
	) t5 ON t5.study_accession = gwas.study_accession
LEFT OUTER JOIN
	( SELECT
		pubmedid,
		COUNT(DISTINCT study_accession) AS "study_perpmid_count"
	FROM    
		gwas
	GROUP BY
		pubmedid
	) t6 ON t6.pubmedid = gwas.pubmedid
ORDER BY
	t3.assn_count DESC
	;
--
UPDATE gwas_counts SET assn_count = 0 WHERE assn_count IS NULL ;
UPDATE gwas_counts SET snp_count = 0 WHERE snp_count IS NULL ;
UPDATE gwas_counts SET gene_r_count = 0 WHERE gene_r_count IS NULL ;
UPDATE gwas_counts SET gene_m_count = 0 WHERE gene_m_count IS NULL ;
--

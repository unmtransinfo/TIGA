--
SELECT
	COUNT(DISTINCT study_accession) AS "study_accession_count",
	COUNT(DISTINCT pubmedid) AS "pubmedid_count"
FROM
	gwas
	;
--
SELECT
	COUNT(*) AS "assn_count",
	COUNT(DISTINCT study_accession) AS "study_accession_count",
	COUNT(DISTINCT pubmedid) AS "pubmedid_count"
FROM
	assn
	;
--
SELECT
	COUNT(*) AS "snp2gene_count",
	CASE
		WHEN reported_or_mapped = 'r' THEN 'Reported'
		WHEN reported_or_mapped = 'm' THEN 'MappedWithin'
		WHEN reported_or_mapped = 'mu' THEN 'MappedUpstream'
		WHEN reported_or_mapped = 'md' THEN 'MappedDownstream'
		ELSE 'Unknown'
	END AS "reported_or_mapped"
FROM
	snp2gene
GROUP BY
	reported_or_mapped
	;
--
-- mysql: '\\1'
-- psql: '\1'
--
SELECT
	COUNT(DISTINCT mapped_trait_uri) AS "trait_count",
	REGEXP_REPLACE(mapped_trait_uri, '^.*/(.*)_.*$','\\1') AS "ontology"
FROM
	trait
GROUP BY
	ontology
	;
--
SELECT
	t.study_accession_count AS "gwas_per_paper",
	COUNT(pubmedid) AS "paper_count"
FROM
	( SELECT
		pubmedid,
		COUNT(DISTINCT study_accession) AS "study_accession_count"
	FROM
		gwas
	GROUP BY
		pubmedid
	) t
GROUP BY
	t.study_accession_count
ORDER BY
	t.study_accession_count DESC
	;
--
SELECT
	t.pubmedid_count AS "paper_per_gwas",
	COUNT(study_accession) AS "study_count"
FROM
	( SELECT
		COUNT(DISTINCT pubmedid) AS "pubmedid_count",
		study_accession
	FROM
		gwas
	GROUP BY
		study_accession
	) t
GROUP BY
	t.pubmedid_count
ORDER BY
	t.pubmedid_count DESC
	;
--
SELECT
	COUNT(DISTINCT study_accession) AS "study_count",
	'TOTAL' AS "platform"
FROM
	gwas
GROUP BY
	platform
UNION
SELECT
	COUNT(DISTINCT study_accession) AS "study_count",
	REGEXP_REPLACE(platform_snps_passing_qc, ' \\[.*$', '') AS "platform"
FROM
	gwas
GROUP BY
	platform
ORDER BY
	platform
	;
--
SELECT
	COUNT(DISTINCT gwas.study_accession) AS "study_count",
	icite.journal
FROM
	gwas,
	icite
WHERE
	gwas.pubmedid = icite.pmid
GROUP BY
	icite.journal
ORDER BY
	study_count DESC
LIMIT 15
	;
--

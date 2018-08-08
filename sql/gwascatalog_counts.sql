-- 
--
SELECT
	COUNT(DISTINCT study_accession) AS "study_accession_count",
	COUNT(DISTINCT pubmedid) AS "pubmedid_count"
FROM
	gwascat g
	;
--
SELECT
	COUNT(*) AS "assn_count",
	COUNT(DISTINCT study_accession) AS "study_accession_count",
	COUNT(DISTINCT pubmedid) AS "pubmedid_count"
FROM
	gwascat_assn ga
	;
--
SELECT
	COUNT(*) AS "snp2gene_count",
	CASE
		WHEN reported_or_mapped = 'r' THEN 'Reported'
		WHEN reported_or_mapped = 'm' THEN 'MappedIn'
		WHEN reported_or_mapped = 'mu' THEN 'MappedUpstream'
		WHEN reported_or_mapped = 'md' THEN 'MappedDownstream'
		ELSE 'Unknown'
	END AS "reported_or_mapped"

FROM
	gwascat_snp2gene gs2g
GROUP BY
	reported_or_mapped
	;
--
-- mysql: '\\1'
-- psql: '\1'
--
SELECT
	COUNT(DISTINCT trait_uri) AS "trait_count",
	REGEXP_REPLACE(trait_uri, '^.*/(.*)_.*$','\1') AS "ontology"
FROM
	gwascat_trait gt
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
		gwascat g
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
		gwascat g
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
	gwascat g
GROUP BY
	platform
UNION
SELECT
	COUNT(DISTINCT study_accession) AS "study_count",
	REGEXP_REPLACE(platform_snps_passing_qc, ' \[.*$', '') AS "platform"
FROM
	gwascat g
GROUP BY
	platform
ORDER BY
	platform
	;
--
SELECT
	COUNT(DISTINCT g.study_accession) AS "study_count",
	gi.journal
FROM
	gwascat g,
	gwascat_icite gi
WHERE
	g.pubmedid = gi.pmid
GROUP BY
	gi.journal
ORDER BY
	study_count DESC
LIMIT 15
	;
--

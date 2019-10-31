--
-- SELECT
-- 	REGEXP_REPLACE(mapped_trait_uri,'^.*/([^_]*)_.*$','\\1') AS "trait_ontology",
-- 	COUNT(DISTINCT mapped_trait),
-- 	COUNT(DISTINCT study_accession)
-- FROM
-- 	trait2study
-- GROUP BY
-- 	trait_ontology
-- 	;
--
SELECT
	mapped_trait_uri,
	mapped_trait,
	COUNT(DISTINCT study_accession) AS "n_study"
FROM
	trait2study
GROUP BY
	mapped_trait
ORDER BY
	n_study DESC
	;
--

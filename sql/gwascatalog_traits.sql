--
SELECT
	SUBSTR(REGEXP_REPLACE(mapped_trait_uri,'^.*/','')||': '||mapped_trait,1,54) AS "mapped_trait",
	COUNT(DISTINCT study_accession) AS "n_gwas"
FROM
	trait
GROUP BY
	mapped_trait
ORDER BY
	n_gwas DESC
	;
--
SELECT
	REGEXP_REPLACE(mapped_trait_uri,'^.*/([^_]*)_.*$','\\1') AS "trait_ontology",
	COUNT(DISTINCT mapped_trait),
	COUNT(DISTINCT study_accession)
FROM
	trait
GROUP BY
	trait_ontology
	;
--

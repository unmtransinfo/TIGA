--
SELECT
	SUBSTR(REGEXP_REPLACE(trait_uri,'^.*/','')||': '||trait,1,54) AS "mapped_trait",
	COUNT(DISTINCT study_accession) AS "n_gwas"
FROM
	trait
GROUP BY
	mapped_trait
ORDER BY
	n_gwas DESC
	;
--

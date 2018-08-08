-- 
-- Find statistically unusual GWAS.
--
--
SELECT
	g.study_accession,
	g.date,
	g.pubmedid,
	g.study,
	t2.trait_count AS "traits_per_study",
	t3.assn_count AS "assns_per_study",
	t1.snp_count AS "snps_per_study",
	t5.gene_count AS "genesr_per_study",
	t6.gene_count AS "genesm_per_study",
	t4.study_count AS "study_per_pub"
FROM
	gwascat g
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gt.trait_uri) AS "trait_count",
		gt.study_accession
	FROM
		gwascat_trait gt
	GROUP BY
		gt.study_accession
	) t2 ON t2.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		pubmedid,
		COUNT(DISTINCT study_accession) AS "study_count"
	FROM    
		gwascat g
	GROUP BY
		pubmedid
	) t4 ON t4.pubmedid = g.pubmedid
LEFT OUTER JOIN
	( SELECT
		COUNT(ga.p_value) AS "assn_count",
		ga.study_accession
	FROM
		gwascat_assn ga
	GROUP BY
		ga.study_accession
	) t3 ON t3.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.snp) AS "snp_count",
		gs2g.study_accession
	FROM
		gwascat_snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped = 'r'
	GROUP BY
		gs2g.study_accession
	) t1 ON t1.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.gsymb) AS "gene_count",
		gs2g.study_accession
	FROM
		gwascat_snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped = 'r'
	GROUP BY
		gs2g.study_accession
	) t5 ON t5.study_accession = g.study_accession
LEFT OUTER JOIN
	( SELECT
		COUNT(DISTINCT gs2g.gsymb) AS "gene_count",
		gs2g.study_accession
	FROM
		gwascat_snp2gene gs2g
	WHERE
		gs2g.reported_or_mapped IN ('m', 'mu', 'md')
	GROUP BY
		gs2g.study_accession
	) t6 ON t6.study_accession = g.study_accession
ORDER BY
	t3.assn_count DESC
	;
--

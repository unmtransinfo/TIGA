-- 
-- 
SELECT COUNT(DISTINCT g.study) AS "study_count" FROM grasp.grasp g ;
SELECT COUNT(DISTINCT g.phenotype_description) AS "phenotype_desc_count" FROM grasp.grasp g ;
SELECT COUNT(DISTINCT g.phenotype_categories_assigned) AS "phenotype_category_count" FROM grasp.grasp g ;
SELECT COUNT(DISTINCT g.journal) AS "journal_count" FROM grasp.grasp g ;
SELECT COUNT(DISTINCT g.pubmedid) AS "pmid_count" FROM grasp.grasp g ;
--
SELECT
	COUNT(DISTINCT g.study) AS "study_count",
	g.in_nhgri_gwas_catalog_82614 AS "in_gwascatalog"
FROM
	grasp.grasp g
GROUP BY
	in_gwascatalog
ORDER BY
	in_gwascatalog DESC
	;
--
SELECT
	COUNT(DISTINCT g.study) AS "study_count",
	g.gwas_description AS "gwas_description"
FROM
	grasp.grasp g
GROUP BY
	gwas_description
ORDER BY
	study_count DESC
	;
--
SELECT
	COUNT(DISTINCT g.pubmedid) AS "pmid_count",
	g.journal AS "journal"
FROM
	grasp.grasp g
GROUP BY
	journal
ORDER BY
	pmid_count DESC
LIMIT 25
	;
--
SELECT COUNT(DISTINCT t.id) AS "tcrd_tid_count" FROM tcrd.tcrd_all t ;
SELECT COUNT(DISTINCT t.id) AS "tcrd2grasp_tid_count" FROM tcrd.tcrd_all t
WHERE t.protein_sym IN (SELECT DISTINCT gsymb FROM grasp.grasp_gene2pmid) ;
--
SELECT COUNT(DISTINCT t.protein_id) AS "tcrd2grasp_pid_count" FROM tcrd.tcrd_all t
WHERE t.protein_sym IN (SELECT DISTINCT gsymb FROM grasp.grasp_gene2pmid) ;
--
SELECT COUNT(DISTINCT t.protein_sym) AS "tcrd2grasp_psymb_count" FROM tcrd.tcrd_all t
WHERE t.protein_sym IN (SELECT DISTINCT gsymb FROM grasp.grasp_gene2pmid) ;
--
SELECT
	COUNT(DISTINCT t.id) AS "tid_count",
        t.tdl
FROM
        tcrd.tcrd_all t
WHERE
	t.protein_sym IN (SELECT DISTINCT gsymb FROM grasp.grasp_gene2pmid)
GROUP BY t.tdl
ORDER BY t.tdl
	;
--

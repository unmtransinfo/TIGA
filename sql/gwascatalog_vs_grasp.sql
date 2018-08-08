-- 
SELECT COUNT(DISTINCT gwascatalog_gene2pmid.gene_symb) AS "gwascatalog_gene_count" FROM gwascatalog.gwascatalog_gene2pmid ;
-- 
SELECT COUNT(DISTINCT gwascatalog_gene2pmid.gene_symb) AS "gwascatalog_only_gene_count" FROM gwascatalog.gwascatalog_gene2pmid
WHERE gwascatalog_gene2pmid.gene_symb NOT IN (SELECT DISTINCT grasp_gene2pmid.gsymb FROM grasp.grasp_gene2pmid ) ;
-- 
SELECT COUNT(DISTINCT grasp_gene2pmid.gsymb) AS "grasp_gene_count" FROM grasp.grasp_gene2pmid ;
-- 
SELECT COUNT(DISTINCT grasp_gene2pmid.gsymb) AS "grasp_only_gene_count" FROM grasp.grasp_gene2pmid
WHERE grasp.grasp_gene2pmid.gsymb NOT IN (SELECT DISTINCT gwascatalog_gene2pmid.gene_symb FROM gwascatalog.gwascatalog_gene2pmid ) ;
-- 
SELECT COUNT(DISTINCT gwascatalog_gene2pmid.gene_symb) AS "common_gene_count" FROM gwascatalog.gwascatalog_gene2pmid
WHERE gwascatalog_gene2pmid.gene_symb IN (SELECT DISTINCT grasp_gene2pmid.gsymb FROM grasp.grasp_gene2pmid ) ;
----
SELECT COUNT(DISTINCT gwascatalog_gene2pmid.pmid) AS "gwascatalog_pmid_count" FROM gwascatalog.gwascatalog_gene2pmid ;
-- 
SELECT COUNT(DISTINCT grasp_gene2pmid.pmid) AS "grasp_pmid_count" FROM grasp.grasp_gene2pmid ;
-- 
SELECT COUNT(DISTINCT gwascatalog_gene2pmid.pmid) AS "common_pmid_count" FROM gwascatalog.gwascatalog_gene2pmid
WHERE gwascatalog_gene2pmid.pmid IN (SELECT DISTINCT grasp_gene2pmid.pmid FROM grasp.grasp_gene2pmid ) ;
--
SELECT COUNT(DISTINCT g.study) AS "grasp_study_count", g.in_nhgri_gwas_catalog_82614 AS "in_gwascatalog" FROM grasp.grasp g
GROUP BY in_gwascatalog ORDER BY in_gwascatalog DESC ;
--

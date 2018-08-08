-- 
-- Export CSV file for GWAS Catalog data.
-- 
SELECT
	t.id::INTEGER AS "tid",
	t.protein_sym AS "psymb",
	t.tdl AS "tdl",
	p.pmid,
	COUNT(g.gwas_cat_diseasetrait) AS "disease_count",
	COUNT(g.gwas_cat_reported_genes) AS "gene_count",
	COUNT(g.gwas_cat_snp_gene_ids) AS "snp_gene_count",
	COUNT(g.gwas_cat_snps) AS "snp_count"
FROM
	tcrd.tcrd_all t,
	gwascatalog.gwascatalog_gene2pmid p,
	gwascatalog.gwascatalog g
WHERE
	t.protein_sym = p.gene_symb
	AND p.pmid = g.gwas_cat_pubmedid
GROUP BY t.id::INTEGER,t.protein_sym,t.tdl,p.pmid
ORDER BY t.id::INTEGER,t.protein_sym,t.tdl,p.pmid
	;
--
--

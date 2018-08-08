-- 
-- Export CSV file for GWAS Catalog data.
-- 
SELECT DISTINCT
	t.id AS "tid",
	t.protein_id AS "pid",
	t.protein_sym AS "psymb",
	t.tdl AS "tdl",
	p.pmid,
	'"'||g.gwas_cat_journal||'"',
	g.gwas_cat_date,
	'"'||g.gwas_cat_diseasetrait||'"',
	'"'||g.gwas_cat_reported_genes||'"',
	'"'||g.gwas_cat_snp_gene_ids||'"',
	g.gwas_cat_snps,
	g.gwas_cat_p_value
FROM
	tcrd.tcrd_all t,
	gwascatalog.gwascatalog_gene2pmid p,
	gwascatalog.gwascatalog g
WHERE
	t.protein_sym = p.gene_symb
	AND p.pmid = g.gwas_cat_pubmedid
ORDER BY t.id,p.pmid
	;
--
--

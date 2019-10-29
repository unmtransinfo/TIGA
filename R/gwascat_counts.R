#!/usr/bin/env Rscript
#############################################################################
### gwascat has one row per study, identified by study_accession.
### gwascat_assn has one row per snp-association, linked to gwascat by study_accession.
### http://www.ebi.ac.uk/gwas/docs/fileheaders
#############################################################################
#iCite annotations via iCite API
#############################################################################
library(readr)
library(RMySQL, quietly = T)
library(dplyr, quietly = T)
library(plotly, quietly = T)
library(webshot, quietly=T)

t0 <- proc.time()

###
#gwascat <- read_delim("~/projects/idg/gwas/data/gwascat_gwas.tsv", "\t", col_types = cols(DATE = col_date(format = "%Y-%m-%d"), DATE_ADDED_TO_CATALOG = col_date(format = "%Y-%m-%d")))
#gwascat_assn <- read_delim("~/projects/idg/gwas/data/gwascat_assn.tsv", "\t", col_types = cols(DATE = col_date(format = "%Y-%m-%d"), DATE_ADDED_TO_CATALOG = col_date(format = "%Y-%m-%d")))
#gwascat_snp2gene <- read_delim("~/projects/idg/gwas/data/gwascat_snp2gene.tsv", "\t", col_types = cols(reported_or_mapped = col_factor(c("r","m","md","mu"))))
#gwascat_trait <- read_delim("~/projects/idg/gwas/data/gwascat_trait.tsv", "\t")
#gwascat_icite <- read_csv("~/projects/idg/gwas/data/gwascat_icite.csv", col_types = cols(is_research_article = col_character()))
###
gwascat <- read.delim("~/projects/idg/gwas/data/gwascat_gwas.tsv", sep="\t")
gwascat$DATE <- as.Date(gwascat$DATE, format="%Y-%m-%d")
gwascat$DATE_ADDED_TO_CATALOG <- as.Date(gwascat$DATE_ADDED_TO_CATALOG, format="%Y-%m-%d")
gwascat_assn <- read.delim("~/projects/idg/gwas/data/gwascat_assn.tsv", sep="\t")
gwascat_assn$DATE <- as.Date(gwascat_assn$DATE, format="%Y-%m-%d")
gwascat_assn$DATE_ADDED_TO_CATALOG <- as.Date(gwascat_assn$DATE_ADDED_TO_CATALOG, format="%Y-%m-%d")
gwascat_snp2gene <- read.delim("~/projects/idg/gwas/data/gwascat_snp2gene.tsv", sep="\t")
gwascat_snp2gene$reported_or_mapped <- as.factor(gwascat_snp2gene$reported_or_mapped)
gwascat_trait <- read.delim("~/projects/idg/gwas/data/gwascat_trait.tsv", sep="\t")
gwascat_icite <- read.csv("~/projects/idg/gwas/data/gwascat_icite.csv")
gwascat_icite$is_research_article <- as.character(gwascat_icite$is_research_article)
###
colnames(gwascat) <- tolower(colnames(gwascat))
colnames(gwascat_assn) <- tolower(colnames(gwascat_assn))
colnames(gwascat_trait) <- tolower(colnames(gwascat_trait))
#
#
###
n_total <- nrow(gwascat)
n_gwas <- length(unique(gwascat$study_accession))
print(sprintf("Total GWAS Catalog rows: %d\n", n_total))
print(sprintf("Total GWAS Catalog unique study accessions: %d\n", n_gwas))
print(sprintf("Total GWAS Catalog unique PubMedIDs: %d\n", length(unique(gwascat$pubmedid))))
#
###
print(sprintf("Highest RCR GWAS publications:\n"))
for (i in 1:30)
{
  print(sprintf("%d. RCR = %.1f ; C/yr = %.1f ; C_count = %d ; nih_pctl = %.1f, (%d) %s: %s\n",
	i, gwascat_icite$relative_citation_ratio[i], gwascat_icite$citations_per_year[i],
	gwascat_icite$citation_count[i], gwascat_icite$nih_percentile[i],
	gwascat_icite$year[i], gwascat_icite$journal[i], gwascat_icite$title[i]))
}

###
#Clean & transform:
gwascat$snps_passing_qc <- as.integer(sub("^.*\\[[^[0-9]*([0-9]*)\\].*$", "\\1", gwascat$`platform_.snps_passing_qc.`))
gwascat$snps_passing_qc_rel <- sub("^.*\\[([^[0-9]*)[0-9]*\\].*$", "\\1", gwascat$`platform_.snps_passing_qc.`)
gwascat$snps_passing_qc_rel[gwascat$snps_passing_qc_rel == ''] <- NA
gwascat_trait <- gwascat_trait[!is.na(gwascat_trait$trait_uri),]

gwascat$sample_size <- gsub(' +', '+', gsub('^ *(.*[^ ]) *$', '\\1', gsub('[^0-9 ]', '', gwascat$initial_sample_size)))
for (i in 1:nrow(gwascat))
{
  gwascat$sample_size[i] <- eval(parse(text = gwascat$sample_size[i]))
}

###
#gwascat_counts table via SQL in Go_gwascat_DbCreate.sh (~30min):
# study_accession
# trait_count: traits per study
# assn_count: associations per study
# snp_count: SNPs per study
# gene_r_count: genes-reported per study
# gene_m_count: genes-mapped per study
# study_count: studies per publication
###
dbcon <- dbConnect(MySQL(), host="localhost", dbname="gwascatalog")
gwascat_counts <- dbGetQuery(dbcon, "SELECT * FROM gwas_counts")
ok <- dbDisconnect(dbcon)
###
write.csv(gwascat_counts, file="data/gwascat_counts.csv", row.names=F)
#
gwascat_counts <- merge(gwascat_counts, gwascat[,c("study_accession", "pubmedid", "sample_size")], all.x = T, all.y = F,  by = "study_accession")

t <- table(gwascat$journal)
v <- as.vector(t)
names(v) <- names(t)
v <- sort(v, decreasing = T)

p1 <- subplot(nrows=2, margin=0.1,
  plot_ly(x=as.numeric(format(gwascat$date, "%Y")), type="histogram"),
  plot_ly(type="bar", orientation="h", y=names(v[1:10]), x=v[1:10]) %>%
    layout(yaxis=list(tickfont=list(size=11)))
  ) %>%
  layout(title = paste0("GWAS Catalog counts (N_gwas = ", n_gwas, ")"),
         font=list(family="Arial", size=14),
         legend=list(x=-0.3, y=1.1),
         annotations = list(x=c(0.4,0.4), y=c(1.0,0.5),
                            xanchor="left", xref="paper", yref="paper", showarrow=F,
                            text=c("gwas_per_year", "gwas_per_journal (top 10)"),
                            font=list(size=12)),
         margin=list(t=100, l=160), showlegend=F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d %H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
p1
export(p=p1, file="data/gwascat_counts_p1.png")
###
#Platform and year.  Fractional counts for multi-platform studies.
plat_year <- data.frame(year = as.numeric(format(gwascat$date, "%Y")), 
                platform = as.character(sub(' \\[.*$', '', gwascat$`platform_.snps_passing_qc.`)))

plat_year_multi <- plat_year[grepl(',', plat_year$platform), ]
plat_year <- plat_year[!grepl(',', plat_year$platform), ]
plat_year$fcount <- as.numeric(1.0)

for (i in 1:nrow(plat_year_multi))
{
  plats <- strsplit(as.character(plat_year_multi$platform[i]), '[, ]+', perl = T)
  n_plats <- length(plats[[1]])
  for (plat in plats)
  {
    plat_year <- rbind(plat_year, data.frame(year = plat_year_multi$year[i], platform = plat, fcount = 1.0/n_plats))
    #print(sprintf("DEBUG: %d\t%s\t%f\n", plat_year_multi$year[i], plat, 1.0/n_plats))
  }
}

plat_year_fcount <- data.frame(year = as.integer(sort(unique(plat_year$year))),
                               Affymetrix = as.numeric(NA),
                               Illumina = as.numeric(NA),
                               Perlegen = as.numeric(NA),
                               NR = as.numeric(NA))

plats <- unique(as.character(plat_year$platform))
for (i in 1:nrow(plat_year_fcount))
{
  year <- plat_year_fcount$year[i]
  for (plat in plats)
  {
    plat_year_fcount[[plat]][i] <- sum(plat_year$fcount[plat_year$year == year & plat_year$platform == plat])
  }
}

p2 <- plot_ly(plat_year_fcount, x = ~year, y = ~NR, name = "NR",
                type = "scatter", mode = "markers", fill = "tonexty") %>%
  add_trace(y = ~Perlegen, name = "Perlegen") %>%
  add_trace(y = ~Affymetrix, name = "Affymetrix") %>%
  add_trace(y = ~Illumina, name = "Illumina") %>%
  layout(title = 'GWAS counts by year|vendor',
         xaxis = list(title = "", showgrid = FALSE),
         yaxis = list(title = "N_study", showgrid = FALSE),
         margin = list(t = 100, l = 100, r = 60),
         font = list(family = "Arial", size = 20),
         legend = list(x = 0.1, y = 0.9),
         showlegend = T) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d\n%H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
p2
export(p=p2, file="data/gwascat_counts_p2.png")
###

max_snp_count <- max(gwascat_counts$snp_count)
max_gene_count <- max(gwascat_counts$gene_r_count)
max_assn_count <- max(gwascat_counts$assn_count)
max_trait_count <- max(gwascat_counts$trait_count)
max_study_perpmid_count <- max(gwascat_counts$study_perpmid_count)

#Merge icite metadata
gwascat_counts <- merge(gwascat_counts, gwascat_icite, all.x=T, all.y=F, by.x="pubmedid", by.y="pmid")

p3 <- subplot(nrows = 2,
  plot_ly(type = "histogram", x = gwascat_counts$snp_count),
  plot_ly(type = "histogram", x = gwascat_counts$gene_r_count),
  plot_ly(type = "histogram", x = gwascat_counts$assn_count),
  plot_ly(type = "histogram", x = gwascat_counts$trait_count)) %>%
  layout(title = paste0("GWAS Catalog: per-study counts<BR> (N_gwas = ", n_gwas, ")"),
	annotations = list(x = c(0.1, 0.6, 0.1, 0.6), y = c(0.9, 0.9, 0.4, 0.4),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("snp_count", "gene_count", "assn_count", "trait_count")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d\n%H:%M:%S"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p3
export(p=p3, file="data/gwascat_counts_p3.png")
#
p4 <- subplot(nrows = 2, margin = 0.08,
  plot_ly(type = "histogram", x = gwascat_icite$citation_count),
  plot_ly(type = "histogram", x = gwascat_icite$nih_percentile),
  plot_ly(type = "histogram", x = gwascat_icite$citations_per_year),
  plot_ly(type = "histogram", x = gwascat_icite$relative_citation_ratio)) %>%
  layout(title = paste0("GWAS Catalog: publication iCite stats<BR> (N_pubs = ", nrow(gwascat_icite), " ; N_gwas = ", n_gwas, ")"),
	annotations = list(x = c(0.1, 0.7, 0.1, 0.7), y = c(0.9, 0.9, 0.3, 0.3),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("citation_count", "nih_percentile", "citations_per_year", "relative_citation_ratio")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d %H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
#
p4
export(p=p4, file="data/gwascat_counts_p4.png")
#
### ROC showing cumulative associations 
### Add hline for top 1% of studies.
gwascat_counts <- gwascat_counts[order(gwascat_counts$assn_count, decreasing=T),]
gwascat_counts$assn_count_cum <- cumsum(gwascat_counts$assn_count)
n_assn_total <- sum(gwascat_counts$assn_count)
i_study_50 <- which.min(abs(gwascat_counts$assn_count_cum - n_assn_total/2))
n_assn_cum_50 <- gwascat_counts$assn_count_cum[i_study_50]
i_study_90 <- which.min(abs(gwascat_counts$assn_count_cum - n_assn_total*.9))
n_assn_cum_90 <- gwascat_counts$assn_count_cum[i_study_90]

p5 <- plot_ly() %>%
  add_trace(x = 1:nrow(gwascat_counts), y = gwascat_counts$assn_count_cum, 
	 type = 'scatter', mode = 'lines', line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
	add_annotations(x = n_gwas, y = n_assn_total*.8, xanchor = "right", showarrow = F, text = paste0("N_gwas = ", n_gwas)) %>%
  add_trace(x = c(0,n_gwas), y = c(n_assn_total/2, n_assn_total/2), type = 'scatter', mode = "lines", line = list(dash = "dot")) %>%
  add_annotations(x = i_study_50, y = n_assn_cum_50, ax = 80, ay = 50, showarrow = T, text = paste0("50%_assn<br>(",i_study_50,",",n_assn_cum_50,")")) %>%
  add_trace(x = c(0,n_gwas), y = c(n_assn_total*.9, n_assn_total*.9), type = 'scatter', mode = "lines", line = list(dash = "dot")) %>%
  add_annotations(x = i_study_90, y = n_assn_cum_90, ax = 50, ay = 50, showarrow = T, text = paste0("90%_assn<br>(",i_study_90,",",n_assn_cum_90,")")) %>%
  add_trace(x = c(0,n_gwas), y = c(n_assn_total, n_assn_total), type = 'scatter', mode = "lines", line = list(dash = "dot")) %>%
  add_annotations(x = n_gwas/5, y = n_assn_total, xanchor = "left", yanchor = "top", showarrow = F, text = paste0("N_assn = ",n_assn_total)) %>%
  add_trace(x = c(n_gwas, n_gwas), y = c(0,n_assn_total), type = 'scatter', mode = "lines", line = list(dash = "dot")) %>%
  	layout(title = "GWAS Catalog: Cumulative associations<br>(studies ordered by per-study associations)",
	xaxis = list(title = "i_study"),
	yaxis = list (title = "N_assn_cum"),
	margin = list(t = 100, l = 60, r = 60),
	font = list(family = "Arial", size = 18),
	showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d\n%H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
#
p5
export(p=p5, file="data/gwascat_counts_p5.png")
#
###


###
#Below not counts, maybe separate file?
gwascat_snp2gene <- gwascat_snp2gene[!(gwascat_snp2gene$gsymb %in% c("intergenic", "Unknown")), ]

t1 <- table(gwascat_trait$trait_uri)
t1 <- t1[order(t1, decreasing = T)]
print("MOST COMMON GWAS TRAITS:\n")
for (i in 1:20)
{
  uri <- names(t1)[i]
  trait <- gwascat_trait$trait[gwascat_trait$trait_uri == uri][1]
  print(sprintf("%2d. (N_study = %3d) [%s] \"%s\"\n", i, t1[i], sub("^.*/", "", uri), substr(trait, 1, 44)))
}
gwascat_trait$ontology <- sub('^.*/([^_]+).*$','\\1',gwascat_trait$trait_uri)
#Pie chart of trait ontologies.
tbl <- table(gwascat_trait$ontology)
ax0 <- list(showline = F, zeroline = F, showticklabels = F, showgrid = F)
p6 <- plot_ly(type = "pie", hole = 0.5, values = tbl[names(tbl)],
       labels = paste0(names(tbl), " (", tbl[names(tbl)], ")") ) %>%
  layout(title = paste0("GWAS Catalog:<br>Traits by ontology"), 
         xaxis = ax0, yaxis = ax0, showlegend = T, margin = list(t=120), legend = list(x=0.4,y=0.5),
         font = list(family = "Arial", size = 12))
p6
export(p=p6, file="data/gwascat_counts_p6.png")
###
dbcon <- dbConnect(MySQL(), host="juniper.health.unm.edu", dbname="tcrd")
sql <- "SELECT
  t.id AS \"tcrd_tid\",  t.name,  t.tdl,  t.fam,  t.idg2,
  p.id AS \"tcrd_pid\",  p.sym AS \"protein_sym\",  p.geneid AS \"protein_geneid\"
FROM
  target t
JOIN t2tc ON t.id = t2tc.target_id
JOIN protein p ON t2tc.protein_id = p.id
WHERE
  p.sym IS NOT NULL"
tcrd <- dbGetQuery(dbcon,sql)
ok <- dbDisconnect(dbcon)

gsyms_tcrd <- unique(tcrd$protein_sym)
print(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))

gsyms_gwascat <- unique(gwascat_snp2gene$gsymb)
writeLines(sprintf("GWASCat Unique gene symbols: %d",length(gsyms_gwascat)))
writeLines(sprintf("GWASCat gene symbols in TCRD: %d", length(intersect(gsyms_gwascat, gsyms_tcrd))))
writeLines(sprintf("GWASCat gene symbols NOT in TCRD: %d", length(setdiff(gsyms_gwascat, gsyms_tcrd))))

tcrd <- merge(tcrd, data.frame(gsym=gsyms_gwascat, in_gwascat=rep(T, length(gsyms_gwascat))),
              by.x="protein_sym", by.y="gsym", all.x=T, all.y=F)
tcrd$in_gwascat[is.na(tcrd$in_gwascat)] <- F
tcrd$idg2 <- as.logical(tcrd$idg2)
t2 <- table(tcrd$tdl[tcrd$in_gwascat])
print(sprintf("%s: %d\n", names(t2), t2))

###

gwascat_gene <- unique(gwascat_snp2gene[,c("study_accession", "gsymb")])
gwascat_gene <- merge(gwascat_gene, tcrd[,c("protein_sym","name","idg2")], all.x=T, all.y=F,
	by.x="gsymb", by.y="protein_sym")
t1 <- table(gwascat_gene$gsymb)
t1 <- t1[order(t1, decreasing = T)]
print("MOST COMMON GWAS GENES:\n")
for (i in 1:20)
{
  gsymb <- names(t1)[i]
  name <- gwascat_gene$name[gwascat_gene$gsymb == gsymb][1]
  idg2 <- tcrd$idg2[tcrd$protein_sym == gsymb][1]
  print(sprintf("%2d. (N_study = %3d) [%5s] \"%s\" (idg2=%s)\n", i, t1[i], gsymb, substr(name, 1, 44), idg2))
}

print("MOST COMMON GWAS GENES (IDG2):\n")
t1 <- table(gwascat_gene$gsymb[gwascat_gene$idg2])
t1 <- t1[order(t1, decreasing = T)]
for (i in 1:10)
{
  gsymb <- names(t1)[i]
  name <- gwascat_gene$name[gwascat_gene$gsymb == gsymb][1]
  idg2 <- tcrd$idg2[tcrd$protein_sym == gsymb][1]
  print(sprintf("%2d. (N_study = %3d) [%5s] \"%s\" (idg2=%s)\n", i, t1[i], gsymb, substr(name, 1, 44), idg2))
}

###
print(sprintf("elapsed time (total): %.2fs",(proc.time()-t0)[3]))

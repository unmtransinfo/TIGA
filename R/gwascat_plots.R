#!/usr/bin/env Rscript
#############################################################################
### gwascat has one row per study, identified by study_accession.
### gwascat_assn has one row per snp-association, linked to gwascat by study_accession.
### http://www.ebi.ac.uk/gwas/docs/fileheaders
#############################################################################
#iCite annotations via iCite API
#############################################################################
library(readr)
library(data.table)
library(RMySQL, quietly = T)
library(dplyr, quietly = T)
library(plotly, quietly = T)
library(webshot, quietly=T)

t0 <- proc.time()

###
gwas <- read_delim("data/gwascat_gwas.tsv", "\t", col_types=cols(.default=col_character(), DATE=col_date(), ASSOCIATION_COUNT=col_integer(), DATE_ADDED_TO_CATALOG=col_date(), study_N=col_integer()))
assn <- read_delim("data/gwascat_assn.tsv", "\t", col_types=cols(.default=col_character(), DATE=col_date(), DATE_ADDED_TO_CATALOG=col_date()), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(), oddsratio=col_double(), beta=col_double())
snp2gene <- read_delim("data/gwascat_snp2gene.tsv", "\t", col_types=cols(.default=col_character(), reported_or_mapped = col_factor()))
trait <- read_delim("data/gwascat_trait.tsv", "\t", col_types=cols(.default=col_character())
icite <- read_delim("data/gwascat_icite.tsv", col_types=cols(.default=col_character())
###
setDT(gwas)
setDT(assn)
setDT(snp2gene)
setDT(trait)
setDT(icite)
###
n_total <- nrow(gwas)
n_gwas <- length(unique(gwas$study_accession))
print(sprintf("Total GWAS Catalog rows: %d\n", n_total))
print(sprintf("Total GWAS Catalog unique study accessions: %d\n", n_gwas))
print(sprintf("Total GWAS Catalog unique PubMedIDs: %d\n", length(unique(gwas$PUBMEDID))))
#
###
print(sprintf("Highest RCR GWAS publications:\n"))
for (i in 1:30)
{
  print(sprintf("%d. RCR = %.1f ; C/yr = %.1f ; C_count = %d ; nih_pctl = %.1f, (%d) %s: %s\n",
	i, icite$relative_citation_ratio[i], icite$citations_per_year[i],
	icite$citation_count[i], icite$nih_percentile[i],
	icite$year[i], icite$journal[i], icite$title[i]))
}

###
#gwas_counts table via SQL in Go_gwascat_DbCreate.sh (~30min):
# study_accession
# trait_count: traits per study
# assn_count: associations per study
# snp_count: SNPs per study
# gene_r_count: genes-reported per study
# gene_m_count: genes-mapped per study
# study_count: studies per publication
###
dbcon <- dbConnect(MySQL(), host="localhost", dbname="gwascatalog")
gwas_counts <- dbGetQuery(dbcon, "SELECT * FROM gwas_counts")
ok <- dbDisconnect(dbcon)
###
#
gwas_counts <- merge(gwas_counts, gwas[,c("study_accession", "pubmedid", "sample_size")], all.x = T, all.y = F,  by = "study_accession")

t <- table(gwas$journal)
v <- as.vector(t)
names(v) <- names(t)
v <- sort(v, decreasing = T)

p1 <- subplot(nrows=2, margin=0.1,
  plot_ly(x=as.numeric(format(gwas$date, "%Y")), type="histogram"),
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
export(p=p1, file="data/gwas_counts.png")
###
#Platform and year.  Fractional counts for multi-platform studies.
plat_year <- data.frame(year = as.numeric(format(gwas$date, "%Y")), 
                platform = as.character(sub(' \\[.*$', '', gwas$`platform_.snps_passing_qc.`)))

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
export(p=p2, file="data/gwas_counts.png")
###

max_snp_count <- max(gwas_counts$snp_count)
max_gene_count <- max(gwas_counts$gene_r_count)
max_assn_count <- max(gwas_counts$assn_count)
max_trait_count <- max(gwas_counts$trait_count)
max_study_perpmid_count <- max(gwas_counts$study_perpmid_count)

#Merge icite metadata
gwas_counts <- merge(gwas_counts, icite, all.x=T, all.y=F, by.x="pubmedid", by.y="pmid")

p3 <- subplot(nrows = 2,
  plot_ly(type = "histogram", x = gwas_counts$snp_count),
  plot_ly(type = "histogram", x = gwas_counts$gene_r_count),
  plot_ly(type = "histogram", x = gwas_counts$assn_count),
  plot_ly(type = "histogram", x = gwas_counts$trait_count)) %>%
  layout(title = paste0("GWAS Catalog: per-study counts<BR> (N_gwas = ", n_gwas, ")"),
	annotations = list(x = c(0.1, 0.6, 0.1, 0.6), y = c(0.9, 0.9, 0.4, 0.4),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("snp_count", "gene_count", "assn_count", "trait_count")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d\n%H:%M:%S"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p3
export(p=p3, file="data/gwas_counts.png")
#
p4 <- subplot(nrows = 2, margin = 0.08,
  plot_ly(type = "histogram", x = icite$citation_count),
  plot_ly(type = "histogram", x = icite$nih_percentile),
  plot_ly(type = "histogram", x = icite$citations_per_year),
  plot_ly(type = "histogram", x = icite$relative_citation_ratio)) %>%
  layout(title = paste0("GWAS Catalog: publication iCite stats<BR> (N_pubs = ",
nrow(icite), " ; N_gwas = ", n_gwas, ")"),
	annotations = list(x = c(0.1, 0.7, 0.1, 0.7), y = c(0.9, 0.9, 0.3, 0.3),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("citation_count", "nih_percentile", "citations_per_year", "relative_citation_ratio")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d %H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
#
p4
export(p=p4, file="data/gwas_counts.png")
#
### ROC showing cumulative associations 
### Add hline for top 1% of studies.
gwas_counts <- gwas_counts[order(gwas_counts$assn_count, decreasing=T),]
gwas_counts$assn_count_cum <- cumsum(gwas_counts$assn_count)
n_assn_total <- sum(gwas_counts$assn_count)
i_study_50 <- which.min(abs(gwas_counts$assn_count_cum - n_assn_total/2))
n_assn_cum_50 <- gwas_counts$assn_count_cum[i_study_50]
i_study_90 <- which.min(abs(gwas_counts$assn_count_cum - n_assn_total*.9))
n_assn_cum_90 <- gwas_counts$assn_count_cum[i_study_90]

p5 <- plot_ly() %>%
  add_trace(x = 1:nrow(gwas_counts), y = gwas_counts$assn_count_cum, 
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
export(p=p5, file="data/gwas_counts.png")
#
###


###
#Below not counts, maybe separate file?
snp2gene <- snp2gene[!(snp2gene$gsymb %in% c("intergenic", "Unknown")), ]

t1 <- table(trait$trait_uri)
t1 <- t1[order(t1, decreasing = T)]
print("MOST COMMON GWAS TRAITS:\n")
for (i in 1:20)
{
  uri <- names(t1)[i]
  trait_this <- trait$trait[trait$trait_uri == uri][1]
  print(sprintf("%2d. (N_study = %3d) [%s] \"%s\"\n", i, t1[i], sub("^.*/", "", uri), substr(trait_this, 1, 44)))
}
trait$ontology <- sub('^.*/([^_]+).*$','\\1',trait$trait_uri)
#Pie chart of trait ontologies.
tbl <- table(trait$ontology)
ax0 <- list(showline = F, zeroline = F, showticklabels = F, showgrid = F)
p6 <- plot_ly(type = "pie", hole = 0.5, values = tbl[names(tbl)],
       labels = paste0(names(tbl), " (", tbl[names(tbl)], ")") ) %>%
  layout(title = paste0("GWAS Catalog:<br>Traits by ontology"), 
         xaxis = ax0, yaxis = ax0, showlegend = T, margin = list(t=120), legend = list(x=0.4,y=0.5),
         font = list(family = "Arial", size = 12))
p6
export(p=p6, file="data/gwas_counts.png")
###
tcrd <- read_delim("data/tcrd_targets.tsv", "\t")
setDT(tcrd)

gsyms_tcrd <- unique(tcrd$protein_sym)
print(sprintf("TCRD targets: %d ; geneSymbols: %d", nrow(tcrd), length(gsyms_tcrd)))

gsyms_gwascat <- unique(snp2gene$gsymb)
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

gwascat_gene <- unique(snp2gene[,c("study_accession", "gsymb")])
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

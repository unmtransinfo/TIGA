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
library(plotly, quietly = T)

message(paste(commandArgs(), collapse=" "))
###
gwas <- read_delim("data/gwascat_gwas.tsv", "\t", col_types=cols(.default=col_character(), DATE=col_date(), ASSOCIATION_COUNT=col_integer(), DATE_ADDED_TO_CATALOG=col_date(), study_N=col_integer()))
assn <- read_delim("data/gwascat_assn.tsv", "\t", col_types=cols(.default=col_character(), DATE=col_date(), DATE_ADDED_TO_CATALOG=col_date(), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(), oddsratio=col_double(), beta=col_double()))
snp2gene <- read_delim("data/gwascat_snp2gene.tsv", "\t", col_types=cols(.default=col_character(), REPORTED_OR_MAPPED = col_factor()))
trait <- read_delim("data/gwascat_trait.tsv", "\t", col_types=cols(.default=col_character()))
icite <- read_delim("data/gwascat_icite.tsv", "\t", col_types=cols(pmid=col_character()))
###
setDT(gwas)
setDT(assn)
setDT(snp2gene)
setDT(trait)
setDT(icite)
###
message(sprintf("Total GWAS Catalog rows: %d", nrow(gwas)))
message(sprintf("Total GWAS Catalog unique study accessions: %d", uniqueN(gwas$STUDY_ACCESSION)))
message(sprintf("Total GWAS Catalog unique PubMedIDs: %d", uniqueN(gwas$PUBMEDID)))
#
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
gwas_counts <- read_delim("data/gwascat_counts.tsv", "\t", col_types = cols(.default=col_double(), study_accession=col_character()))
setDT(gwas_counts)
###
#
gwas_counts <- merge(gwas_counts, gwas[, .(STUDY_ACCESSION, PUBMEDID, study_N)], all.x = T, all.y = F,  by.x = "study_accession", by.y="STUDY_ACCESSION")

###
print(sprintf("Highest RCR GWAS publications:\n"))
icite <- icite[order(-relative_citation_ratio)]
myrange <- 1:10
writeLines(sprintf("%d. RCR = %.1f ; C/yr = %.1f ; C_count = %d ; nih_pctl = %.1f, (%d) %s: %s",
	myrange, icite$relative_citation_ratio[myrange], icite$citations_per_year[myrange],
	icite$citation_count[myrange], icite$nih_percentile[myrange],
	icite$year[myrange], icite$journal[myrange], icite$title[myrange]))
###

max_snp_count <- max(gwas_counts$snp_count)
max_gene_count <- max(gwas_counts$gene_r_count)
max_assn_count <- max(gwas_counts$assn_count)
max_trait_count <- max(gwas_counts$trait_count)
max_study_perpmid_count <- max(gwas_counts$study_perpmid_count)

#Merge icite metadata
gwas_counts <- merge(gwas_counts, icite, all.x=T, all.y=F, by.x="PUBMEDID", by.y="pmid")

subplot(nrows = 2,
  plot_ly(type = "histogram", x = gwas_counts$snp_count),
  plot_ly(type = "histogram", x = gwas_counts$gene_r_count),
  plot_ly(type = "histogram", x = gwas_counts$assn_count),
  plot_ly(type = "histogram", x = gwas_counts$trait_count)) %>%
  layout(title = paste0("GWAS Catalog: per-study counts<BR> (N_gwas = ", uniqueN(gwas$STUDY_ACCESSION), ")"),
	annotations = list(x = c(0.1, 0.6, 0.1, 0.6), y = c(0.9, 0.9, 0.4, 0.4),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("snp_count", "gene_count", "assn_count", "trait_count")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d\n%H:%M:%S"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
subplot(nrows = 2, margin = 0.08,
  plot_ly(type = "histogram", x = icite$citation_count),
  plot_ly(type = "histogram", x = icite$nih_percentile),
  plot_ly(type = "histogram", x = icite$citations_per_year),
  plot_ly(type = "histogram", x = icite$relative_citation_ratio)) %>%
  layout(title = paste0("GWAS Catalog: publication iCite stats<BR> (N_pubs = ", nrow(icite)),
	annotations = list(x = c(0.1, 0.7, 0.1, 0.7), y = c(0.9, 0.9, 0.3, 0.3),
          xanchor = "left", xref = "paper", yref = "paper", showarrow = F,
          text = c("citation_count", "nih_percentile", "citations_per_year", "relative_citation_ratio")),
	margin = list(t = 100, l = 60, r = 60),
         font = list(family = "Arial", size = 18), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d %H:%M:%S"), showarrow=F, x=1.0, y=1.2, xref="paper", yref="paper")
#
### ROC showing cumulative associations 
### Add hline for top 1% of studies.
gwas_counts <- gwas_counts[order(-assn_count)]
gwas_counts$assn_count_cum <- cumsum(gwas_counts$assn_count)
n_assn_total <- sum(gwas_counts$assn_count)
i_study_50 <- which.min(abs(gwas_counts$assn_count_cum - n_assn_total/2))
n_assn_cum_50 <- gwas_counts$assn_count_cum[i_study_50]
i_study_90 <- which.min(abs(gwas_counts$assn_count_cum - n_assn_total*.9))
n_assn_cum_90 <- gwas_counts$assn_count_cum[i_study_90]
n_gwas <- uniqueN(gwas$STUDY_ACCESSION)

plot_ly() %>%
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
###


###
#Below not counts, maybe separate file?
snp2gene <- snp2gene[!(GSYMB %in% c("intergenic", "Unknown"))]

trait_counts <- trait[, .(N_study = uniqueN(STUDY_ACCESSION), MAPPED_TRAIT=first(MAPPED_TRAIT)), by=MAPPED_TRAIT_URI]
setorder(trait_counts, -N_study)
print("MOST COMMON GWAS TRAITS:\n")
for (i in 1:20)
{
  print(sprintf("%2d. (N_study = %3d) [%s] \"%s\"\n", i, trait_counts[i, N_study], trait_counts[i, sub("^.*/", "", MAPPED_TRAIT_URI)], trait_counts[i, substr(MAPPED_TRAIT, 1, 44)]))
}

###
tcrd <- read_delim("data/tcrd_targets.tsv", "\t")
setDT(tcrd)

print(sprintf("TCRD targets: %d ; geneSymbols: %d; ENSGs: %d", nrow(tcrd), uniqueN(tcrd$tcrdGeneSymbol), uniqueN(tcrd$ensemblGeneId)))

writeLines(sprintf("GWASCat Unique gene symbols: %d", uniqueN(snp2gene$GSYMB)))
writeLines(sprintf("GWASCat gene symbols in TCRD: %d", length(intersect(snp2gene[, unique(GSYMB)], tcrd[, unique(tcrdGeneSymbol)]))))
writeLines(sprintf("GWASCat gene symbols NOT in TCRD: %d", length(setdiff(snp2gene[, unique(GSYMB)], tcrd[, unique(tcrdGeneSymbol)]))))

###

tcrd <- merge(tcrd, data.frame(gsym=snp2gene[, unique(GSYMB)], in_gwascat=T), by.x="tcrdGeneSymbol", by.y="gsym", all.x=T, all.y=F)
tcrd[is.na(in_gwascat), in_gwascat := F]
tcrd[, idgList := as.logical(idgList)]


gene_counts <- snp2gene[, .(N_study = uniqueN(STUDY_ACCESSION)), by=GSYMB]
gene_counts <- merge(gene_counts, tcrd[, .(tcrdGeneSymbol, tcrdTargetName, idgList)], all.x=T, all.y=F, by.x="GSYMB", by.y="tcrdGeneSymbol")
setorder(gene_counts, -N_study)
print("MOST COMMON GWAS GENES:\n")
for (i in 1:20)
{
  print(sprintf("%2d. (N_study = %3d) [%5s] \"%s\" (idgList=%s)\n", i, gene_counts[i, N_study], gene_counts[i, GSYMB], gene_counts[i, substr(tcrdTargetName, 1, 44)], gene_counts[i, idgList]))
}

###
print("MOST COMMON GWAS GENES (idgList):\n")
gene_counts <- snp2gene[, .(N_study = uniqueN(STUDY_ACCESSION)), by=GSYMB]
gene_counts <- merge(gene_counts, tcrd[(idgList), .(tcrdGeneSymbol, tcrdTargetName, idgList)], all.x=F, all.y=F, by.x="GSYMB", by.y="tcrdGeneSymbol")
setorder(gene_counts, -N_study)
print("MOST COMMON GWAS GENES:\n")
for (i in 1:20)
{
  print(sprintf("%2d. (N_study = %3d) [%5s] \"%s\" (idgList=%s)\n", i, gene_counts[i, N_study], gene_counts[i, GSYMB], gene_counts[i, substr(tcrdTargetName, 1, 44)], gene_counts[i, idgList]))
}

###

tdl_counts <- tcrd[(in_gwascat), .N, by=TDL]
tdl_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_counts <- tdl_counts[order(TDL)]
print(tdl_counts)
#Pie chart of TDL counts.
ax0 <- list(showline=F, zeroline=F, showticklabels=F, showgrid=F)
plot_ly(tdl_counts, type="pie", hole=0.2, values=~N, labels=~factor(TDL), marker = list(colors = c("green", "blue", "red", "#222222")),
  textposition = 'inside', textinfo = 'label+value+percent',
        #insidetextfont = list(color = '#FFFFFF'),
        hoverinfo = 'text', text = ~paste0(TDL, "\n", N, " genes")) %>%
  layout(title=sprintf("MAPPED GENES:<br>Target Development Levels (TDLs)<br>N_total = %d", sum(tdl_counts$N)), 
         xaxis=ax0, yaxis=ax0, showlegend=F, margin=list(t=120),
         font=list(family="Arial", size=12))

fam_counts <- tcrd[(in_gwascat), .N, by=tcrdTargetFamily]
setorder(fam_counts, -N)
print(fam_counts)
plot_ly(fam_counts, type="pie", hole=0.5, values=~N, labels=~factor(tcrdTargetFamily),
  textposition = 'inside', textinfo = 'label+value+percent',
        hoverinfo = 'text', text = ~paste0(tcrdTargetFamily, "\n", N, " genes")) %>%
  layout(title=sprintf("MAPPED GENES:<br>Target family<br>N_total = %d", sum(fam_counts$N)), 
         xaxis=ax0, yaxis=ax0, showlegend=T, margin=list(t=120), legend=list(x=0.4, y=0.5),
         font=list(family="Arial", size=12))

###
# Table of all TDL, Family combos:

tdl_fam_counts <- tcrd[(in_gwascat), .N, by=c("TDL", "tcrdTargetFamily")]
tdl_fam_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_fam_counts <- dcast(tdl_fam_counts, formula = tcrdTargetFamily ~ TDL, value.var = "N", fill = 0)
setorder(tdl_fam_counts, -Tclin, -Tchem, -Tbio, -Tdark)
setnames(tdl_fam_counts, old="tcrdTargetFamily", new="Family")
print(tdl_fam_counts)
write_delim(tdl_fam_counts, "data/tdl_fam_counts.tsv", "\t")



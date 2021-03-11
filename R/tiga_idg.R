#!/usr/bin/env Rscript
#############################################################################
### IDG counts and plots
#############################################################################
library(readr)
library(data.table)
library(plotly, quietly = T)
library(RMySQL, quietly = T)

message(paste(commandArgs(), collapse=" "))

dbcon <- dbConnect(MySQL(), host="tcrd.kmc.io", port=3306, dbname="tcrd660", user="tcrd", password="")

sql="SELECT DISTINCT
	target.id AS \"tcrdTargetId\",
	target.name AS \"tcrdTargetName\",
	target.fam AS \"tcrdTargetFamily\",
	target.tdl AS \"TDL\",
	target.idg AS \"idgList\",
	protein.sym AS \"geneSymbol\",
	protein.geneid AS \"ncbiGeneId\",
  protein.uniprot,
  protein.stringid AS \"ensemblProteinId\",
  protein.dtoid,
  protein.dtoclass,
	xref.value AS \"ensemblGeneId\"
FROM
	target
JOIN
	t2tc ON t2tc.target_id = target.id
JOIN
	protein ON protein.id = t2tc.protein_id
JOIN
	xref ON xref.protein_id = protein.id
WHERE
	xref.xtype = 'Ensembl' AND xref.value REGEXP '^ENSG'
ORDER BY
	protein.sym"

tcrd <- dbGetQuery(dbcon, sql)
setDT(tcrd)

print(sprintf("TCRD targets: %d ; geneSymbols: %d; ENSGs: %d; ENSPs: %d; UniProts: %d", uniqueN(tcrd$tcrdTargetId), 
  uniqueN(tcrd$geneSymbol), uniqueN(tcrd$ensemblGeneId),  uniqueN(tcrd$ensemblProteinId),
  uniqueN(tcrd$uniprot)))
###
tcrd_dto <- read_delim("data/TCRDv6.4_DTO.tsv", "\t")
setDT(tcrd_dto)
setnames(tcrd_dto, old=c("STRING ID"), new=c("ensemblProteinId"))
tcrd_dto <- unique(tcrd_dto[, .(ensemblProteinId, DTO_Lvl1, DTO_Lvl2, DTO_Lvl3, DTO_Lvl4, DTO_Lvl5, DTO_Lvl6, DTO_Class)])
###
#
tiga <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), 
	n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
	geneNtrait=col_integer(), geneNstudy=col_integer(),
	traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(),
	meanRank=col_double(), meanRankScore=col_double()))
setDT(tiga)
print(sprintf("TIGA geneSymbols: %d; ENSGs: %d", uniqueN(tiga$geneSymbol), uniqueN(tiga$ensemblId)))
###

tcrd <- merge(tcrd, data.frame(ensemblId=tiga[, unique(ensemblId)], in_tiga=T), by.x="ensemblGeneId", by.y="ensemblId", all.x=T, all.y=F)
tcrd[is.na(in_tiga), in_tiga := F]
tcrd[, idgList := as.logical(idgList)]

###

tdl_counts <- tcrd[(in_tiga), .N, by=TDL]
tdl_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_counts <- tdl_counts[order(TDL)]
print(tdl_counts)
#Pie chart of TDL counts.
ax0 <- list(showline=F, zeroline=F, showticklabels=F, showgrid=F)
plot_ly(tdl_counts, type="pie", hole=0.2, values=~N, labels=~factor(TDL), marker = list(colors = c("green", "blue", "red", "#222222")),
  textposition = 'inside', textinfo = 'label+value+percent',
        #insidetextfont = list(color = '#FFFFFF'),
        hoverinfo = 'text', text = ~paste0(TDL, "\n", N, " genes")) %>%
  layout(title=sprintf("TIGA MAPPED GENES:<br>Target Development Levels (TDLs)<br>N_total = %d", sum(tdl_counts$N)), 
         xaxis=ax0, yaxis=ax0, showlegend=F, margin=list(t=120),
         font=list(family="Arial", size=12))

fam_counts <- tcrd[(in_tiga), .N, by=tcrdTargetFamily]
setorder(fam_counts, -N)
print(fam_counts)
plot_ly(fam_counts, type="pie", hole=0.5, values=~N, labels=~factor(tcrdTargetFamily),
  textposition = 'inside', textinfo = 'label+value+percent',
        hoverinfo = 'text', text = ~paste0(tcrdTargetFamily, "\n", N, " genes")) %>%
  layout(title=sprintf("TIGA MAPPED GENES:<br>Target family<br>N_total = %d", sum(fam_counts$N)), 
         xaxis=ax0, yaxis=ax0, showlegend=T, margin=list(t=120), legend=list(x=0.4, y=0.5),
         font=list(family="Arial", size=12))
###

tcrd <- merge(tcrd, tcrd_dto, by="ensemblProteinId", allow.cartesian=T, all.x=T, all.y=F)

dto_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=DTO_Lvl2]
setorder(dto_counts, -N)
print(dto_counts)
plot_ly(dto_counts, type="pie", hole=0.1, values=~N, labels=~factor(DTO_Lvl2),
  textposition = 'inside', textinfo = 'label+value+percent',
        hoverinfo = 'text', text = ~paste0(DTO_Lvl2, "\n", N, " genes")) %>%
  layout(title=sprintf("TIGA MAPPED GENES:<br>DTO_Lvl2<br>N_total = %d", sum(dto_counts$N)), 
         xaxis=ax0, yaxis=ax0, showlegend=T, margin=list(t=120),
         font=list(family="Arial", size=12))

###
# Table of all gene counts for TDL, TargetFamily combos:
tcrd[tcrdTargetFamily=="IC", tcrdTargetFamily := "Ion channel"]
tcrd[tcrdTargetFamily=="NR", tcrdTargetFamily := "Nuclear receptor"]
FAMS = c("GPCR", "Ion channel", "Kinase", "Enzyme", "Transporter", "Nuclear receptor")
writeLines(sprintf("TCRD family being merged into Other: \"%s\"", 
  tcrd[(!tcrdTargetFamily %in% FAMS), unique(tcrdTargetFamily)]))
#
tdl_fam_counts <- tcrd[, .(N = uniqueN(ensemblProteinId)), by=c("TDL", "tcrdTargetFamily")]
tdl_fam_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_fam_counts <- dcast(tdl_fam_counts, formula = tcrdTargetFamily ~ TDL, value.var = "N", fill = 0)
setnames(tdl_fam_counts, old="tcrdTargetFamily", new="Family")
tdl_fam_counts[!(Family %in% FAMS), Family := "Other"]
tdl_fam_counts <- tdl_fam_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tdl_fam_counts <- rbindlist(list(tdl_fam_counts, data.table(Family="Total", Tclin=sum(tdl_fam_counts$Tclin), Tchem=sum(tdl_fam_counts$Tchem), 
  Tbio=sum(tdl_fam_counts$Tbio), Tdark=sum(tdl_fam_counts$Tdark))))
tdl_fam_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tdl_fam_counts[, Family := factor(Family, levels=c(FAMS, "Other", "Total"), ordered=T)]
setorder(tdl_fam_counts, Family)
print(tdl_fam_counts)
write_delim(tdl_fam_counts, "data/tdl_fam_counts.tsv", "\t")
#
# Table of all gene counts (in TIGA) for TDL, TargetFamily combos:
tiga_tdl_fam_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=c("TDL", "tcrdTargetFamily")]
tiga_tdl_fam_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tiga_tdl_fam_counts <- dcast(tiga_tdl_fam_counts, formula = tcrdTargetFamily ~ TDL, value.var = "N", fill = 0)
setnames(tiga_tdl_fam_counts, old="tcrdTargetFamily", new="Family")
tiga_tdl_fam_counts[!(Family %in% FAMS), Family := "Other"]
tiga_tdl_fam_counts <- tiga_tdl_fam_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tiga_tdl_fam_counts <- rbindlist(list(tiga_tdl_fam_counts, data.table(Family="Total", Tclin=sum(tiga_tdl_fam_counts$Tclin), 
  Tchem=sum(tiga_tdl_fam_counts$Tchem), Tbio=sum(tiga_tdl_fam_counts$Tbio), Tdark=sum(tiga_tdl_fam_counts$Tdark))))
tiga_tdl_fam_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tiga_tdl_fam_counts[, Family := factor(Family, levels=c(FAMS, "Other", "Total"), ordered=T)]
setorder(tiga_tdl_fam_counts, Family)
print(tiga_tdl_fam_counts)
write_delim(tiga_tdl_fam_counts, "data/tdl_fam_counts_TIGA.tsv", "\t")
#
write_delim(data.table(tiga_tdl_fam_counts)[, `:=`(
	Tclin = sprintf("%d / %d", Tclin, tdl_fam_counts$Tclin),
	Tchem = sprintf("%d / %d", Tchem, tdl_fam_counts$Tchem),
	Tbio = sprintf("%d / %d", Tbio, tdl_fam_counts$Tbio),
	Tdark = sprintf("%d / %d", Tdark, tdl_fam_counts$Tdark),
  Total = sprintf("%d / %d", Total, tdl_fam_counts$Total))], "data/tdl_fam_counts_MERGED.tsv", "\t")

###
# Same but with DTO Level2 classes.
# Table of all gene counts for TDL, DTO_Lvl2 combos:
FAMS_DTOL2 <- c(
  "G-protein coupled receptor",
  "Ion channel",
  "Kinase",
	"Calcium-binding protein",
	"Cell-cell junction",
	"Cell adhesion",
	"Cellular structure",
	"Chaperone",
	"Enzyme modulator",
	"Enzyme",
	"Epigenetic regulator",
	"Extracellular structure",
	"Immune response",
	"Nuclear receptor",
	"Nucleic acid binding",
	"Transcription factor",
	"Transporter",
	"Receptor",
	"Signaling",
	"Storage",
	"Surfactant")

tdl_dto_counts <- tcrd[, .(N = uniqueN(ensemblProteinId)), by=c("TDL", "DTO_Lvl2")]
tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_dto_counts <- dcast(tdl_dto_counts, formula = DTO_Lvl2 ~ TDL, value.var = "N", fill = 0)
setnames(tdl_dto_counts, old="DTO_Lvl2", new="Family")
tdl_dto_counts[!(Family %in% FAMS_DTOL2), Family := "Other"]
tdl_dto_counts <- tdl_dto_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tdl_dto_counts <- rbindlist(list(tdl_dto_counts, data.table(Family="Total", Tclin=sum(tdl_dto_counts$Tclin), Tchem=sum(tdl_dto_counts$Tchem), 
  Tbio=sum(tdl_dto_counts$Tbio), Tdark=sum(tdl_dto_counts$Tdark))))
tdl_dto_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTOL2, "Other", "Total"), ordered=T)]
setorder(tdl_dto_counts, Family)
print(tdl_dto_counts)
write_delim(tdl_dto_counts, "data/tdl_dto_counts.tsv", "\t")
#
# Table of all gene counts (in TIGA) for TDL, DTO_Lvl2 combos:
tiga_tdl_dto_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=c("TDL", "DTO_Lvl2")]
tiga_tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tiga_tdl_dto_counts <- dcast(tiga_tdl_dto_counts, formula = DTO_Lvl2 ~ TDL, value.var = "N", fill = 0)
setnames(tiga_tdl_dto_counts, old="DTO_Lvl2", new="Family")
tiga_tdl_dto_counts[!(Family %in% FAMS_DTOL2), Family := "Other"]
tiga_tdl_dto_counts <- tiga_tdl_dto_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tiga_tdl_dto_counts <- rbindlist(list(tiga_tdl_dto_counts, data.table(Family="Total", Tclin=sum(tiga_tdl_dto_counts$Tclin), 
  Tchem=sum(tiga_tdl_dto_counts$Tchem), Tbio=sum(tiga_tdl_dto_counts$Tbio), Tdark=sum(tiga_tdl_dto_counts$Tdark))))
tiga_tdl_dto_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tiga_tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTOL2, "Other", "Total"), ordered=T)]
setorder(tiga_tdl_dto_counts, Family)
print(tiga_tdl_dto_counts)
write_delim(tiga_tdl_dto_counts, "data/tdl_dto_counts_TIGA.tsv", "\t")
#
write_delim(data.table(tiga_tdl_dto_counts)[, `:=`(
	Tclin = sprintf("%d / %d", Tclin, tdl_dto_counts$Tclin),
	Tchem = sprintf("%d / %d", Tchem, tdl_dto_counts$Tchem),
	Tbio = sprintf("%d / %d", Tbio, tdl_dto_counts$Tbio),
	Tdark = sprintf("%d / %d", Tdark, tdl_dto_counts$Tdark),
  Total = sprintf("%d / %d", Total, tdl_dto_counts$Total))], "data/tdl_dto_counts_MERGED.tsv", "\t")


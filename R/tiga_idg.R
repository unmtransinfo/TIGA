#!/usr/bin/env Rscript
#############################################################################
### IDG counts - Generates table as in TIGA paper.
#############################################################################
library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)

ODIR <- "data/20210212"
#
ifile_tcrd <- ifelse((length(args)>0), args[1], paste0(ODIR, "/tcrd_targets.tsv")) #BioClients.idg API
ifile_tcrd2dto <- ifelse((length(args)>1), args[2], paste0(ODIR, "/tcrd2dto.tsv")) #From BioClients.idg.tcrd.Client listTargetsByDTO
#ifile_tcrd2dto <- ifelse((length(args)>1), args[2], paste0(ODIR, "/TCRDv6.4_DTO.tsv")) #From ???
ifile_tiga <- ifelse((length(args)>2), args[3], paste0(ODIR, "/gt_stats.tsv.gz"))
ofile_fam_counts <- ifelse((length(args)>3), args[4], paste0(ODIR, "/tdl_fam_counts_TIGA.tsv"))
ofile_fam_counts_merged <- ifelse((length(args)>4), args[5], paste0(ODIR, "/tdl_fam_counts_MERGED.tsv", "\t"))
ofile_dto_counts <- ifelse((length(args)>5), args[6], paste0(ODIR, "/tdl_dto_counts_TIGA.tsv"))
ofile_dto_counts_merged <- ifelse((length(args)>6), args[7], paste0(ODIR, "/tdl_dto_counts_MERGED.tsv", "\t"))
#
message(sprintf("INPUT TCRD file: %s", ifile_tcrd))
message(sprintf("INPUT TCRD2DTO file: %s", ifile_tcrd2dto))
message(sprintf("INPUT TIGA GT_STATS file: %s", ifile_tiga))
message(sprintf("OUTPUT FAM counts file: %s", ofile_fam_counts))
message(sprintf("OUTPUT FAM MERGED (with totals) counts file: %s", ofile_fam_counts_merged))
message(sprintf("OUTPUT DTO counts file: %s", ofile_dto_counts))
message(sprintf("OUTPUT DTO MERGED (with totals) counts file: %s", ofile_dto_counts_merged))
#
tcrd <- read_delim(ifile_tcrd, "\t", na=c("", "NA", "NULL"), col_types=cols(.default=col_character(), idgList=col_logical()))
setDT(tcrd)
tcrd[, `:=`(dtoId=NULL, dtoClass=NULL)]

message(sprintf("TCRD targets: %d; geneSymbols: %d; ENSGs: %d; ENSPs: %d; UniProts: %d", uniqueN(tcrd$tcrdTargetId), 
  uniqueN(tcrd$tcrdGeneSymbol), uniqueN(tcrd$ensemblGeneId),  uniqueN(tcrd$ensemblProteinId), uniqueN(tcrd$uniprotId)))
###
tcrd_dto <- read_delim(ifile_tcrd2dto, "\t")
setDT(tcrd_dto)
#setnames(tcrd_dto, old=c("STRING ID"), new=c("ensemblProteinId"))
#tcrd_dto <- unique(tcrd_dto[, .(ensemblProteinId, DTO_Lvl1, DTO_Lvl2, DTO_Lvl3, DTO_Lvl4, DTO_Lvl5, DTO_Lvl6, DTO_Class)])
###
#
tiga <- read_delim(ifile_tiga, '\t', col_types=cols(.default=col_character(), 
	n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
	geneNtrait=col_integer(), geneNstudy=col_integer(),
	traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(),
	meanRank=col_double(), meanRankScore=col_double()))
setDT(tiga)
message(sprintf("TIGA geneSymbols: %d; ENSGs: %d", uniqueN(tiga$geneSymbol), uniqueN(tiga$ensemblId)))
###

tcrd <- merge(tcrd, data.frame(ensemblId=tiga[, unique(ensemblId)], in_tiga=T), by.x="ensemblGeneId", by.y="ensemblId", all.x=T, all.y=F)
tcrd[is.na(in_tiga), in_tiga := F]
tcrd[, idgList := as.logical(idgList)]
#
tdl_counts <- tcrd[(in_tiga), .N, by=TDL]
tdl_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_counts <- tdl_counts[order(TDL)]
print(tdl_counts)
#
fam_counts <- tcrd[(in_tiga), .N, by=tcrdTargetFamily]
setorder(fam_counts, -N)
print(fam_counts)
#
#tcrd <- merge(tcrd, tcrd_dto, by="ensemblProteinId", allow.cartesian=T, all.x=T, all.y=F)
tcrd <- merge(tcrd, tcrd_dto[dtoGeneration==4 & dtoClass!="protein", .(uniprotId, dtoId, dtoClass, dtoGeneration)], by="uniprotId", allow.cartesian=T, all.x=T, all.y=F)
#
#dto_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=DTO_Lvl2]
dto_counts <- tcrd[(in_tiga), .(N = uniqueN(uniprotId)), by=c("dtoId", "dtoClass")]
setorder(dto_counts, -N)
print(dto_counts)
###
# Table of all gene counts for TDL, TargetFamily combos:
tcrd[tcrdTargetFamily=="IC", tcrdTargetFamily := "Ion channel"]
tcrd[tcrdTargetFamily=="NR", tcrdTargetFamily := "Nuclear receptor"]
FAMS = c("GPCR", "Ion channel", "Kinase", "Enzyme", "Transporter", "Nuclear receptor")
message(sprintf("TCRD family being merged into Other: \"%s\"", tcrd[(!tcrdTargetFamily %in% FAMS), unique(tcrdTargetFamily)]))
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
write_delim(tdl_fam_counts, ofile_fam_counts, "\t")
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
write_delim(tiga_tdl_fam_counts, ofile_fam_counts, "\t")
#
write_delim(data.table(tiga_tdl_fam_counts)[, `:=`(
	Tclin = sprintf("%d / %d", Tclin, tdl_fam_counts$Tclin),
	Tchem = sprintf("%d / %d", Tchem, tdl_fam_counts$Tchem),
	Tbio = sprintf("%d / %d", Tbio, tdl_fam_counts$Tbio),
	Tdark = sprintf("%d / %d", Tdark, tdl_fam_counts$Tdark),
  Total = sprintf("%d / %d", Total, tdl_fam_counts$Total))], ofile_fam_counts_merged, "\t")

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

FAMS_DTO <- c(
	"G-protein coupled receptor",
	"Ion channel",
	"Kinase",
	"Protein kinase",
	"Enzyme",
	"Transporter",
	"Transcription factor",
	"Nucleic acid binding",
	"Epigenetic regulator",
	"Transferase")

#tdl_dto_counts <- tcrd[, .(N = uniqueN(ensemblProteinId)), by=c("TDL", "DTO_Lvl2")]
tdl_dto_counts <- tcrd[, .(N = uniqueN(ensemblProteinId)), by=c("TDL", "dtoClass")]

tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
#tdl_dto_counts <- dcast(tdl_dto_counts, formula = DTO_Lvl2 ~ TDL, value.var = "N", fill = 0)
#setnames(tdl_dto_counts, old="DTO_Lvl2", new="Family")
#tdl_dto_counts[!(Family %in% FAMS_DTOL2), Family := "Other"]
#
tdl_dto_counts <- dcast(tdl_dto_counts, formula = dtoClass ~ TDL, value.var = "N", fill = 0)
setnames(tdl_dto_counts, old="dtoClass", new="Family")
tdl_dto_counts[!(Family %in% FAMS_DTO), Family := "Other"]
#
tdl_dto_counts <- tdl_dto_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tdl_dto_counts <- rbindlist(list(tdl_dto_counts, data.table(Family="Total", Tclin=sum(tdl_dto_counts$Tclin), Tchem=sum(tdl_dto_counts$Tchem), 
  Tbio=sum(tdl_dto_counts$Tbio), Tdark=sum(tdl_dto_counts$Tdark))))
tdl_dto_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
#tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTOL2, "Other", "Total"), ordered=T)]
tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTO, "Other", "Total"), ordered=T)]

setorder(tdl_dto_counts, Family)
print(tdl_dto_counts)
write_delim(tdl_dto_counts, ofile_dto_counts, "\t")
#
# Table of all gene counts (in TIGA) for TDL, DTO_Lvl2 combos:
#tiga_tdl_dto_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=c("TDL", "DTO_Lvl2")]
tiga_tdl_dto_counts <- tcrd[(in_tiga), .(N = uniqueN(ensemblProteinId)), by=c("TDL", "dtoClass")]

tiga_tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
#tiga_tdl_dto_counts <- dcast(tiga_tdl_dto_counts, formula = DTO_Lvl2 ~ TDL, value.var = "N", fill = 0)
#setnames(tiga_tdl_dto_counts, old="DTO_Lvl2", new="Family")
tiga_tdl_dto_counts <- dcast(tiga_tdl_dto_counts, formula = dtoClass ~ TDL, value.var = "N", fill = 0)
setnames(tiga_tdl_dto_counts, old="dtoClass", new="Family")
#
tiga_tdl_dto_counts[!(Family %in% FAMS_DTO), Family := "Other"]
tiga_tdl_dto_counts <- tiga_tdl_dto_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tiga_tdl_dto_counts <- rbindlist(list(tiga_tdl_dto_counts, data.table(Family="Total", Tclin=sum(tiga_tdl_dto_counts$Tclin), 
  Tchem=sum(tiga_tdl_dto_counts$Tchem), Tbio=sum(tiga_tdl_dto_counts$Tbio), Tdark=sum(tiga_tdl_dto_counts$Tdark))))
tiga_tdl_dto_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tiga_tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTO, "Other", "Total"), ordered=T)]
setorder(tiga_tdl_dto_counts, Family)
print(tiga_tdl_dto_counts)
write_delim(tiga_tdl_dto_counts, ofile_fam_counts, "\t")
#
dto_counts_merged <- data.table(tiga_tdl_dto_counts)[, `:=`(
	Tclin = sprintf("%d / %d", Tclin, tdl_dto_counts$Tclin),
	Tchem = sprintf("%d / %d", Tchem, tdl_dto_counts$Tchem),
	Tbio = sprintf("%d / %d", Tbio, tdl_dto_counts$Tbio),
	Tdark = sprintf("%d / %d", Tdark, tdl_dto_counts$Tdark),
  Total = sprintf("%d / %d", Total, tdl_dto_counts$Total))]
write_delim(dto_counts_merged, ofile_dto_counts_merged, "\t")
#

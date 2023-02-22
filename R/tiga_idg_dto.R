#!/usr/bin/env Rscript
#############################################################################
### IDG counts - tables of TIGA and total gene counts by DTO and TDL.
### Count by ENSG (not ENSP) to agree with gene counts elsewhere.
#############################################################################
library(readr)
library(data.table)

message(paste(commandArgs(), collapse=" "))
args <- commandArgs(trailingOnly=TRUE)
#
if (interactive()) {
  rel_y <- as.integer(readline(prompt="Enter RELEASE_YEAR: "))
  rel_m <- as.integer(readline(prompt="Enter RELEASE_MONTH: "))
  rel_d <- as.integer(readline(prompt="Enter RELEASE_DAY: "))
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (length(args)==3) {
  rel_y <- as.integer(args[1])
  rel_m <- as.integer(args[2])
  rel_d <- as.integer(args[3])
  ODIR <- sprintf("data/%d%02d%02d", rel_y, rel_m, rel_d)
} else if (file.exists("LATEST_RELEASE.txt")) {
  GC_REL <- trimws(read_file("LATEST_RELEASE.txt"))
  rel_y <- as.integer(sub("\\-.*$", "", GC_REL))
  rel_m <- as.integer(sub("\\d+\\-(\\d+)\\-.*$", "\\1", GC_REL))
  rel_d <- as.integer(sub("\\d+\\-\\d+\\-(\\d+).*$", "\\1", GC_REL))
  message(sprintf("LATEST_RELEASE: %s", GC_REL))
  ODIR <- sprintf("data/%s", gsub("\\-", "", GC_REL))
} else {
  message("ERROR: Syntax: tiga_idg_dto.R RELEASE_YEAR RELEASE_MONTH RELEASE_DAY")
  quit()
}
#
ifile_tcrd <- paste0(ODIR, "/tcrd_targets.tsv") #BioClients.idg API
ifile_tcrd2dto <- paste0(ODIR, "/tcrd2dto.tsv") #From BioClients.idg.tcrd.Client listTargetsByDTO
ifile_tiga <- paste0(ODIR, "/gt_stats.tsv.gz")
ofile_dto_counts <- paste0(ODIR, "/tdl_dto_counts_TIGA.tsv")
ofile_dto_counts_merged <- paste0(ODIR, "/tdl_dto_counts_MERGED.tsv")
#
message(sprintf("INPUT TCRD file: %s", ifile_tcrd))
message(sprintf("INPUT TCRD2DTO file: %s", ifile_tcrd2dto))
message(sprintf("INPUT TIGA GT_STATS file: %s", ifile_tiga))
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
print(fam_counts)
#
tcrd <- merge(tcrd, tcrd_dto[dtoGeneration==4 & dtoClass!="protein", .(uniprotId, dtoId, dtoClass, dtoGeneration)], by="uniprotId", allow.cartesian=T, all.x=T, all.y=F)
#
dto_counts <- tcrd[(in_tiga), .(N = uniqueN(uniprotId)), by=c("dtoId", "dtoClass")]
setorder(dto_counts, -N)
print(dto_counts)
###
# Table of ALL TCRD gene counts for TDL, DTO-Gen4 combos:
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
#
tdl_dto_counts <- tcrd[, .(N = uniqueN(ensemblGeneId)), by=c("TDL", "dtoClass")]
tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
#
tdl_dto_counts <- dcast(tdl_dto_counts, formula = dtoClass ~ TDL, value.var = "N", fill = 0)
setnames(tdl_dto_counts, old="dtoClass", new="Family")
tdl_dto_counts[!(Family %in% FAMS_DTO), Family := "Other"]
#
tdl_dto_counts <- tdl_dto_counts[, .(Tclin=sum(Tclin), Tchem=sum(Tchem), Tbio=sum(Tbio), Tdark=sum(Tdark)), by=Family]
tdl_dto_counts <- rbindlist(list(tdl_dto_counts, data.table(Family="Total", Tclin=sum(tdl_dto_counts$Tclin), Tchem=sum(tdl_dto_counts$Tchem), 
  Tbio=sum(tdl_dto_counts$Tbio), Tdark=sum(tdl_dto_counts$Tdark))))
tdl_dto_counts[, Total := sum(Tclin, Tchem, Tbio, Tdark), by=Family]
tdl_dto_counts[, Family := factor(Family, levels=c(FAMS_DTO, "Other", "Total"), ordered=T)]
setorder(tdl_dto_counts, Family)
print(tdl_dto_counts)
write_delim(tdl_dto_counts, ofile_dto_counts, "\t")
#
###
# Table of TIGA gene counts for TDL, DTO-Gen4 combos:
message(sprintf("Unique ENSGs in TIGA: %d", length(unique(tcrd[(in_tiga), ensemblGeneId]))))
# Need to ensure DTO-Gen4 and TDL are disjoint sets, no multi-membership for counts.
tcrd_tiga <- unique(tcrd[(in_tiga), .(TDL=first(TDL), dtoClass=first(dtoClass)), by=ensemblGeneId])
#
tiga_tdl_dto_counts <- tcrd_tiga[, .(N = uniqueN(ensemblGeneId)), by=c("TDL", "dtoClass")]
tiga_tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
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
write_delim(tiga_tdl_dto_counts, ofile_dto_counts, "\t")
#
# Merging TIGA and ALL tables for concise presentation.
dto_counts_merged <- data.table(tiga_tdl_dto_counts)[, `:=`(
	Tclin = sprintf("%d / %d", Tclin, tdl_dto_counts$Tclin),
	Tchem = sprintf("%d / %d", Tchem, tdl_dto_counts$Tchem),
	Tbio = sprintf("%d / %d", Tbio, tdl_dto_counts$Tbio),
	Tdark = sprintf("%d / %d", Tdark, tdl_dto_counts$Tdark),
  Total = sprintf("%d / %d", Total, tdl_dto_counts$Total))]
write_delim(dto_counts_merged, ofile_dto_counts_merged, "\t")
#
###
# Now try instead with new DTO top-level class membership input file (from iu_idsl_jena).
tcrd <- read_delim(ifile_tcrd, "\t", na=c("", "NA", "NULL"), col_types=cols(.default=col_character(), idgList=col_logical()))
setDT(tcrd)
dto_tlcm <- read_delim("data/dto_complete_merged_toplevelsuperclassmembership.tsv", "\t", na=c("", "NA", "na", "NULL", "null"))
setDT(dto_tlcm)
#dto_tlcm <- dto_tlcm[!(is.na(SuperClassUriLev_0) | is.na(SuperClassLabelLev_0))]
message(sprintf("dtoIds in TCRD: %d; in DTO file: %d; in common: %d", tcrd[!is.na(dtoId), uniqueN(dtoId)], dto_tlcm[, uniqueN(id)],
                length(intersect(tcrd[!is.na(dtoId), unique(dtoId)], dto_tlcm[, unique(id)]))))
tcrd_tlcm <- merge(tcrd, dto_tlcm, by.x="dtoId", by.y="id", all.x=T, all.y=F)
#z <- tcrd_tlcm[(!is.na(dtoId) & !is.na(uri) & !is.na(SuperClassUriLev_0))]
tdl_dto_counts <- tcrd_tlcm[, .(N = uniqueN(ensemblGeneId)), by=c("TDL", "SuperClassLabelLev_2")]
tdl_dto_counts[, TDL := factor(TDL, levels=c("Tclin", "Tchem", "Tbio", "Tdark"), ordered=T)]
tdl_dto_counts <- dcast(tdl_dto_counts, formula = SuperClassLabelLev_2 ~ TDL, value.var = "N", fill = 0)
setnames(tdl_dto_counts, old="SuperClassLabelLev_2", new="Family")
#tdl_dto_counts[!(Family %in% FAMS_DTO), Family := "Other"]
print(tdl_dto_counts)

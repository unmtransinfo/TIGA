#!/usr/bin/env Rscript
###
library(readr)
library(data.table)


###
# JensenLab DISEASES:
#Experiments means (1) DistiLD (GWAS) and (2) COSMIC (somatic mutations in cancer).
diseases_exp <- read_delim(paste0(Sys.getenv("HOME"), "/../data/JensenLab/DISEASES/human_disease_experiments_full.tsv"), "\t", col_names=c("geneEnsp", "geneSymbol", "doId", "doName", "DISEASES_source", "DISEASES_evidence", "DISEASES_confidence"))
setDT(diseases_exp)
diseases_exp <- diseases_exp[order(-DISEASES_confidence), .SD, by=c("doId", "doName")]

###
# Selected DO diseases which map to EFO, in both datasets, and associate with numerous genes:
# Mapping from OxO! https://www.ebi.ac.uk/spot/oxo/
do2efo <- read_delim("data/oxo_do2efo.tsv", "\t")
setDT(do2efo)
do2efo <- do2efo[distance == 1]
setnames(do2efo, old=c("curie_id", "label", "mapped_curie", "mapped_label"), new=c("doId", "doName", "efoId", "efoName"))
do2efo[, `:=`(mapping_source_prefix=NULL, mapping_target_prefix=NULL, distance=NULL)] #distance all 1
do2efo[, efoId := sub(":", "_", efoId)]

#
tiga <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), 
n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
geneNtrait=col_integer(), geneNstudy=col_integer(),
traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(),
geneMuScore=col_double(), geneMuRank=col_integer(),
traitMuScore=col_double(), traitMuRank=col_integer()))
setDT(tiga)                                                           
#
tiga <- merge(tiga, do2efo[, .(efoId, doId, doName)], by="efoId")

###
I_max <- 100L
results_diseases <- data.table(
	efoId=rep(as.character(NA), I_max),
	efoName=rep(as.character(NA)),
	doId=rep(as.character(NA), I_max),
	doName=rep(as.character(NA), I_max),
	tiga_Ngenes=rep(as.integer(NA), I_max),
	diseases_Ngenes=rep(as.integer(NA), I_max),
	inCommon_Ngenes=rep(as.integer(NA), I_max),
	inCommon_pct=rep(as.integer(NA), I_max)
)
Ngenes_min <- 100L
i <- 0L
cat("i\tefoId\tefoName\tdoId\tdoName\ttiga_Ngenes\tdiseases_Ngenes\tNgenes_InCommon\tNgenes_InCommonPct\n")
for (efoId_this in tiga[, unique(efoId)]) {
  if (tiga[efoId==efoId_this, uniqueN(ensemblId)] < Ngenes_min) { next }
  tiga_this <- tiga[efoId==efoId_this]
  doId_this <- tiga_this[, first(doId)]
  if (diseases_exp[doId==doId_this, uniqueN(geneSymbol)] < Ngenes_min) { next }
  #if (!(doId_this %in% diseases_exp$doId)) { next }
  i <- i + 1L
  efoName_this <- tiga_this[, first(trait)]
  doName_this <- tiga_this[, first(doName)]
  tiga_Ngenes <- tiga_this[, uniqueN(ensemblId)]
  diseases_exp_this <- diseases_exp[doId == doId_this]
  diseases_Ngenes <- diseases_exp_this[, uniqueN(geneSymbol)]
  genes_in_common <- intersect(diseases_exp_this$geneSymbol, tiga_this$geneSymbol)
  cat(sprintf("%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f%%\n", i, efoId_this, efoName_this, doId_this, doName_this, tiga_Ngenes, diseases_Ngenes, length(genes_in_common), 100*length(genes_in_common)/min(diseases_Ngenes, tiga_Ngenes)))
  set(results_diseases, i, names(results_diseases), list(efoId_this, efoName_this, doId_this, doName_this, tiga_Ngenes, diseases_Ngenes, length(genes_in_common), 100*length(genes_in_common)/min(diseases_Ngenes, tiga_Ngenes)))
  if (i==I_max) { break }
}
setorder(results_diseases, -tiga_Ngenes, na.last=T)
results_diseases <- results_diseases[!is.na(efoId)]
write_delim(results_diseases, "data/benchmarks_results_diseases.tsv", "\t")
write_delim(results_diseases[, .(efoId)], "data/benchmarks.efoId", col_names=F)
#
###
# Totals:
message(sprintf("TOTAL TIGA_Ngenes: %d; DISEASES_Ngenes: %d; InCommon: %d (median %.1f%%)\n",
	results_diseases[, sum(tiga_Ngenes)],
	results_diseases[, sum(diseases_Ngenes)],
	results_diseases[, sum(inCommon_Ngenes)],
	results_diseases[, median(inCommon_pct)]))
#
###
# OpenTargets (via OpenTargets API)
# python3 -m BioClients.opentargets.Client searchAssociations -v --idtype 'disease' --i data/benchmarks.efoId --o data/benchmarks_opentargets.tsv
#
# id, is_direct, target_id, gene_name, gene_symbol, disease_id, disease_efo_label, score_overall
opentargets <- read_delim("data/benchmarks_opentargets.tsv", "\t", col_types=cols(.default=col_character(), is_direct=col_logical()))
setDT(opentargets)

for (name in names(opentargets)) {
  if ((grepl("^assn_score_.*\\.", name) | grepl("^evidence_count_.*\\.", name)) & !grepl("\\.ot_genetics_portal$", name)) {
    opentargets[[name]] <- NULL
  }
}
opentargets$assn_score_source.ot_genetics_portal <- as.numeric(opentargets$assn_score_source.ot_genetics_portal)
opentargets$evidence_count_source.ot_genetics_portal <- as.numeric(opentargets$evidence_count_source.ot_genetics_portal)
opentargets <- opentargets[assn_score_source.ot_genetics_portal>0.0 & evidence_count_source.ot_genetics_portal>0.0]
opentargets[, `:=`(id=NULL, is_direct=NULL)]
setnames(opentargets, old=c("target_id", "gene_name", "gene_symbol", "disease_id", "disease_efo_label", "score_overall"), new=c("ensgId", "geneName", "geneSymbol", "efoId", "efoName", "otOverallScore"))

###
results_ot <- data.table(
	efoId=rep(as.character(NA), I_max),
	efoName=rep(as.character(NA)),
	doId=rep(as.character(NA), I_max),
	doName=rep(as.character(NA), I_max),
	tiga_Ngenes=rep(as.integer(NA), I_max),
	ot_Ngenes=rep(as.integer(NA), I_max),
	inCommon_Ngenes=rep(as.integer(NA), I_max),
	inCommon_pct=rep(as.integer(NA), I_max)
)
i <- 0L
cat("i\tefoId\tefoName\tdoId\tdoName\ttiga_Ngenes\topentargets_Ngenes\tNgenes_InCommon\tNgenes_InCommonPct\n")
for (efoId_this in tiga[, unique(efoId)]) {
  if (tiga[efoId==efoId_this, uniqueN(ensemblId)] < Ngenes_min) { next }
  tiga_this <- tiga[efoId==efoId_this]
  doId_this <- tiga_this[, first(doId)]
  if (opentargets[efoId==efoId_this, uniqueN(geneSymbol)] < Ngenes_min) { next }
  i <- i + 1L
  efoName_this <- tiga_this[, first(trait)]
  doName_this <- tiga_this[, first(doName)]
  tiga_Ngenes <- tiga_this[, uniqueN(ensemblId)]
  opentargets_this <- opentargets[efoId == efoId_this]
  opentargets_this <- opentargets_this[assn_score_source.ot_genetics_portal>=.8]
  opentargets_Ngenes <- opentargets_this[, uniqueN(geneSymbol)]
  genes_in_common <- intersect(opentargets_this$geneSymbol, tiga_this$geneSymbol)
  cat(sprintf("%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f%%\n", i, efoId_this, efoName_this, doId_this, doName_this, tiga_Ngenes, opentargets_Ngenes, length(genes_in_common), 100*length(genes_in_common)/min(opentargets_Ngenes, tiga_Ngenes)))
  set(results_ot, i, names(results_ot), list(efoId_this, efoName_this, doId_this, doName_this, tiga_Ngenes, opentargets_Ngenes, length(genes_in_common), 100*length(genes_in_common)/min(opentargets_Ngenes, tiga_Ngenes)))
  if (i==I_max) { break }
}
setorder(results_ot, -tiga_Ngenes, na.last=T)
results_ot <- results_ot[!is.na(efoId)]
write_delim(results_ot, "data/benchmarks_results_opentargets.tsv", "\t")

###
# Totals:
message(sprintf("TOTAL TIGA_Ngenes: %d; OPENTARGETS_Ngenes: %d; InCommon: %d (median %.1f%%)\n",
	results_ot[, sum(tiga_Ngenes)],
	results_ot[, sum(ot_Ngenes)],
	results_ot[, sum(inCommon_Ngenes)],
	results_ot[, median(inCommon_pct)]))
#

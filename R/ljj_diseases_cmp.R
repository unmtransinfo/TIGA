#!/usr/bin/env Rscript
###
library(readr)
library(data.table)

###

#Experiments means (1) DistiLD (GWAS) and (2) COSMIC (somatic mutations in cancer).
diseases_exp <- read_delim("/home/data/JensenLab/DISEASES/human_disease_experiments_full.tsv", "\t", 
                             col_names=c("geneEnsp", "geneSymbol", "doId", "doName", "source", "evidence", "confidence"))
setDT(diseases_exp)
diseases_exp <- diseases_exp[order(-confidence), .SD, by=c("doId", "doName")]

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
results <- data.table(
	efoId=rep(as.character(NA), I_max),
	efoName=rep(as.character(NA)),
	doId=rep(as.character(NA), I_max),
	doName=rep(as.character(NA), I_max),
	diseases_Ngenes=rep(as.integer(NA), I_max),
	tiga_Ngenes=rep(as.integer(NA), I_max),
	inCommon_Ngenes=rep(as.integer(NA), I_max),
	inCommon_pct=rep(as.integer(NA), I_max)
)
Ngenes_min <- 100L
i <- 0L
for (efoId_this in tiga[, unique(efoId)]) {
  if (tiga[efoId==efoId_this, uniqueN(ensemblId)] < Ngenes_min) { next }
  tiga_this <- tiga[efoId==efoId_this]
  doId_this <- tiga_this[, first(doId)]
  if (!(doId_this %in% diseases_exp$doId)) { next }
  i <- i + 1L
  efoName_this <- tiga_this[, first(trait)]
  doName_this <- tiga_this[, first(doName)]
  tiga_Ngenes <- tiga_this[, uniqueN(ensemblId)]
  diseases_exp_this <- diseases_exp[doId == doId_this]
  diseases_Ngenes <- diseases_exp_this[, uniqueN(geneSymbol)]
  genes_in_common <- intersect(diseases_exp_this$geneSymbol, tiga_this$geneSymbol)
  print(sprintf("%d. %s (\"%s\") ; %s (\"%s\"); TIGA_Ngenes: %d; DISEASES_Ngenes: %d; InCommon: %d (%.1f%%)", i, efoId_this, efoName_this, 
                  doId_this, doName_this, tiga_Ngenes, diseases_Ngenes, length(genes_in_common), 
                  100*length(genes_in_common)/min(diseases_Ngenes, tiga_Ngenes)))

  set(results, i, names(results), list(efoId_this, efoName_this, doId_this, doName_this, tiga_Ngenes, diseases_Ngenes, length(genes_in_common), 100*length(genes_in_common)/min(diseases_Ngenes, tiga_Ngenes)))
  if (i==I_max) { break }
}
setorder(results, -tiga_Ngenes, na.last=T)
results <- results[!is.na(efoId)]
write_delim(results, "data/ljj_diseases_cmp_results.tsv", "\t")

###
# Totals:
message(sprintf("TOTAL TIGA_Ngenes: %d; DISEASES_Ngenes: %d; InCommon: %d (median %.1f%%)",
	results[, sum(tiga_Ngenes)],
	results[, sum(diseases_Ngenes)],
	results[, sum(inCommon_Ngenes)],
	results[, median(inCommon_pct)]))
#
# library(plotly)

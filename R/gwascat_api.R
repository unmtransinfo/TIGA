# Process files from GWAS Catalog API.
# What data add value relative to the downloads?

library(readr)
library(data.table)

studyAssns <- read_delim("data/gwascat_StudyAssociations.tsv.gz", "\t", col_types = cols(.default=col_character(), standardError = col_double(), pvalue = col_double(), betaNum = col_double(), orPerCopyNum = col_double(), multiSnpHaplotype = col_logical(), snpInteraction = col_logical(), pvalueMantissa = col_double(), pvalueExponent = col_double(), locus_haplotypeSnpCount = col_integer(), allele_genomeWide = col_logical(), allele_limitedList = col_logical()))
setDT(studyAssns)
message(sprintf("associations: %d", nrow(studyAssns)))

studyAssns[, rsId := sub("-.", "", allele_riskAlleleName)]
#studyAssns <- studyAssns[!is.na(rsId)] #Why missing? See code: ignore reported genes. 

studyAssns[, betaUnit := sub("^\\(?[Zz][- ][Ss]core\\)?$", "Z-score", betaUnit)]
betaUnits <- studyAssns[, .(.N), by="betaUnit"][order(-N)]

###
# "LINC" = Long Intergenic Non-protein-Coding
studySnps <- read_delim("data/gwascat_Snps.tsv.gz", "\t", col_types=cols(.default=col_character(), merged=col_logical(), lastUpdateDate=col_datetime(),   genomicContext_isIntergenic = col_logical(), genomicContext_isUpstream = col_logical(), genomicContext_isDownstream = col_logical(), genomicContext_distance = col_double(), genomicContext_isClosestGene = col_logical(), loc_chromosomePosition = col_double()))
setDT(studySnps)
message(sprintf("SNPs: %d", nrow(studySnps)))

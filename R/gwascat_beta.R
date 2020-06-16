#!/usr/bin/env Rscript
#############################################################################
### Input file from gwascat_assn.R,
### which separated OR from BETA, via heuristic: If all OR_or_BETA values
### for a study are >=1, assume OR.
###
### From Catalog help:
###`OR or BETA`: Reported odds ratio or beta-coefficient associated with
### strongest SNP risk allele. Note that if an OR <1 is reported this is
### inverted, along with the reported allele, so that all ORs included in
### the Catalog are >1. Appropriate unit and increase/decrease are included
### for beta coefficients.
#############################################################################
library(readr)
library(data.table, quietly=T)
library(plotly, quietly=T)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
  (ifile <- args[1])
  (ofile <- args[2])
} else if (length(args)==0) {
  ifile <- "data/gwascat_assn.tsv"
  ofile <- "data/gwascat_beta.tsv"
} else {
  message("ERROR: Syntax: gwascat_assn.R ASSNFILE OFILE\n\t...or no args for defaults.")
  quit()
}
writeLines(sprintf("Input: %s", ifile))
writeLines(sprintf("Output: %s", ofile))

assn <- read_delim(ifile, "\t", col_types=cols(.default=col_character(), oddsratio=col_double(), beta=col_double()))
setDT(assn)
setnames(assn, old=c("MAPPED_TRAIT_URI", "MAPPED_TRAIT"), new=c("TRAIT_URI", "TRAIT"))

writeLines(sprintf("Studies with OR: %d", assn[!is.na(oddsratio), uniqueN(STUDY_ACCESSION)]))
writeLines(sprintf("Studies with BETA: %d", assn[!is.na(beta), uniqueN(STUDY_ACCESSION)]))

###
qs <- quantile(assn[["beta"]][!is.na(assn[["beta"]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("beta %4s-ile: %9.1f", names(qs), qs))

qs <- quantile(assn[["oddsratio"]][!is.na(assn[["oddsratio"]])], c(0, .25, .5, .75, seq(.9, 1, .01)))
writeLines(sprintf("oddsratio %4s-ile: %9.1f", names(qs), qs))

#Cleaning
assn[, ci_95_text := `95%_CI_(TEXT)`]
assn[, ci_95_text := sub("[\u0080\u0093\u0094]+", "-", ci_95_text)]
assn[, ci_95_text := gsub("[Ââ]", "", ci_95_text)]
assn[, ci_95_text := sub("^\\((.*)\\)\\s*$", "[\\1]", ci_95_text)] #both parens to brackets
assn[, ci_95_text := sub("^\\((.*)\\]\\s*$", "[\\1]", ci_95_text)] #open paren to bracket
assn[, ci_95_text := sub("^\\[(.*)\\)\\s*$", "[\\1]", ci_95_text)] #close paren to bracket
assn[, ci_95_text := sub("^([0-9\\-][0-9\\.]*)\\-([0-9\\.]+)\\]$", "[\\1-\\2]", ci_95_text)] #add missing open bracket
assn[, ci_95_text := sub("^\\[([0-9\\-][0-9\\.]*)\\-([0-9\\.]+)$", "[\\1-\\2]", ci_95_text)] #add missing close bracket
assn[, ci_95_text := sub("^([0-9\\-][0-9\\.]*)\\-([0-9\\.]+)$", "[\\1-\\2]", ci_95_text)] #add both missing brackets
assn[, ci_95_text := gsub("^\\[(.*)\\s*\\-\\s*(.*)\\](.*)$", "[\\1-\\2]\\3", ci_95_text)] #remove spaces in ci

assn[, ci_95_text := sub("^([0-9\\-]+\\.[0-9]+)\\.([0-9]+\\.[0-9]+)$", "[\\1-\\2]", ci_95_text)] #weirdo
#
###
# Heuristic: all units include "(in|de)crease"
assn[, beta_units := sub("^.*\\]\\s*", "", ci_95_text)]
assn[, ci_min := as.numeric(sub("^\\[(.+)\\-(.+)\\].*$", "\\1", ci_95_text))]
assn[, ci_max := as.numeric(sub("^\\[(.+)\\-(.+)\\].*$", "\\2", ci_95_text))]
#
#Flag unparsable beta units
assn[, bad_beta := (!is.na(ci_95_text) & !grepl("^\\[.*\\]$", ci_95_text)  & !grepl("(in|de)crease", ci_95_text))]
print(unique(assn[(bad_beta==TRUE), .(bad = sub("^.*\\] *", "", ci_95_text))]))
#
#Flag NR CIs (not recorded)
assn[, ci_is_nr := (grepl("^\\[NR\\]", ci_95_text))]
#
assn$beta_units[!grepl("(in|de)crease", assn$beta_units)] <- NA
#
#Top BETA units
tbl <- data.table(table(assn$beta_units))
setnames(tbl, c("beta_units", "Freq"))
setorder(tbl, -Freq)
writeLines(sprintf("%2d. %8d: %s", 1:20, tbl[1:20, Freq], tbl[1:20, beta_units]))
writeLines(sprintf("Pct of associations either \"unit increase\" or \"unit decrease\": %.1f%% (%d/%d)", 
    100*tbl[grepl("unit (in|de)crease", beta_units), sum(Freq)]/tbl[, sum(Freq)], tbl[grepl("unit (in|de)crease", beta_units), sum(Freq)], tbl[, sum(Freq)]))
###
#Can we generate z-scores? Some?
#Maybe for some traits. Assume same traits have same units?
#Another approach: beta value counts for mu scoring/ranking.
xxx <- assn[(grepl("^unit (in|de)crease", beta_units) & ci_is_nr==F), .(TRAIT_URI, TRAIT, `95%_CI_(TEXT)`, ci_95_text, beta, beta_units, ci_min, ci_max)]
setorder(xxx, -ci_min, na.last=T)

anns <- c(sprintf("median:%.3f<br>mean:%.3f<br>range:[%.3f-%.3f]", median(xxx$ci_min, na.rm=T), mean(xxx$ci_min, na.rm=T), min(xxx$ci_min, na.rm=T), max(xxx$ci_min, na.rm=T)),
  sprintf("median:%.3f<br>mean:%.3f<br>range:[%.3f-%.3f]", median(xxx$ci_max, na.rm=T), mean(xxx$ci_max, na.rm=T), min(xxx$ci_max, na.rm=T), max(xxx$ci_max, na.rm=T))
          )
plot_ly() %>%
  add_boxplot(name="CI_MIN", data=xxx, y=~ci_min) %>%
  add_boxplot(name="CI_MAX", data=xxx, y=~ci_max) %>%
  layout(title="Effect size beta 95%_CI",
         xaxis=list(title="", tickangle=45), yaxis=list(title="", type="log"),
         margin = list(t=100, l=100),
         font = list(family="monospace", size=16)) %>%
  add_annotations(text=anns, 
        showarrow=F, x=c(0, 1), y=c(1, 1), xref="paper", yref="paper")

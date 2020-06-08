#!/usr/bin/env Rscript
###
# Higher level groupings could provide UI visual channel for trait classification.
# General groupings can include child trait-gene associations as needed for expected behavior.
###
library(readr)
library(data.table, quietly=T)

efo_groups <- read_delim("/Users/data/gwascatalog/tiga_data/efo_groups.tsv", "\t")
setDT(efo_groups)

setorder(efo_groups, level, -N_sub)

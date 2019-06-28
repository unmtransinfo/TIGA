#!/usr/bin/env Rscript
#############################################################################
### UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE for mapping
### confidence scoring. Also get CONTEXT, functionalClass in API. 
### INTERGENIC? "RISK ALLELE FREQUENCY"?
### MERGED?
#############################################################################
library(readr)
library(data.table)
library(plotly)

ifile_assn <- "data/gwascat_assn.tsv"

assn <- read_delim(ifile_assn, "\t", col_types=cols(.default=col_character(), INTERGENIC=col_logical(), `P-VALUE`=col_double(), PVALUE_MLOG=col_double(), `OR_or_BETA`=col_double(), UPSTREAM_GENE_DISTANCE=col_integer(), DOWNSTREAM_GENE_DISTANCE=col_integer(), DATE=col_date(format="%Y-%m-%d"), DATE_ADDED_TO_CATALOG=col_date(format="%Y-%m-%d")))
setDT(assn)

message(sprintf("Associations with UPSTREAM_GENE_DISTANCE: %d (%.1f%%)", assn[!is.na(UPSTREAM_GENE_DISTANCE), .N], 100*assn[!is.na(UPSTREAM_GENE_DISTANCE), .N]/nrow(assn)))
message(sprintf("Associations with DOWNSTREAM_GENE_DISTANCE: %d (%.1f%%)", assn[!is.na(DOWNSTREAM_GENE_DISTANCE), .N], 100*assn[!is.na(DOWNSTREAM_GENE_DISTANCE), .N]/nrow(assn)))
plots <- list(
  plot_ly(name="Histogram", type="histogram", x=pmin(assn$UPSTREAM_GENE_DISTANCE, assn$DOWNSTREAM_GENE_DISTANCE, na.rm=T),
              histnorm="probability") %>%
  add_annotations(text=sprintf("N_[BP>200k] = %d", assn[pmin(UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, na.rm=T)>2e5, .N]), x=1, xanchor="right", y=.1, xref="paper", yref="paper", showarrow=F),
  plot_ly() %>%
    add_trace(name="ReLU", type="scatter", x=seq(0, 1e5, by=1000),  mode="lines+markers", marker=list(size=5), y=1-seq(0, 1e5, by=1000)/1e5) %>%
    add_trace(name="Sigmoid", x=seq(0, 1e5, by=1000), y = exp(-seq(-.5e5, .5e5, by=1000)/1e4)/(1+exp(-seq(-.5e5, .5e5, by=1000)/1e4)))
  )
subplot(plots, nrows=2, shareX=T, titleX=T, titleY=T, margin=.05) %>%
  layout(title="Association/SNP (UP|DOWN)STREAM_GENE_DISTANCE<br>(closest gene)",
         margin = list(t=100, r=50),
         xaxis=list(title="BPs", range=c(0,2e5)), yaxis=list(title="N_snp"),
         font=list(family="monospace", size=18))


qtl <- quantile(pmin(assn$UPSTREAM_GENE_DISTANCE, assn$DOWNSTREAM_GENE_DISTANCE, na.rm=T), na.rm=T)
message(sprintf("%4sile: %d\n", names(qtl), qtl))

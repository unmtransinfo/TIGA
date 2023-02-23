#!/usr/bin/env Rscript
###
library(readr)
library(data.table)
library(plotly)

tigastats <- read_delim(paste0(Sys.getenv("HOME"), "/Dropbox/UNM/TIGA/TIGA_Statistics.tsv"), delim="\t", skip=2,
  col_names = c("Release", "Studies", "Publications", "Traits", "Genes", "GC_Studies", "GC_Traits", "Comment"))
setDT(tigastats)

tigastats[, Comment := NULL]
tigastats[, Release := as.Date(Release)]

fig <- plot_ly(tigastats, x=~Release, mode="markers")
fig <- fig %>% add_lines(y=~Studies, name="Studies", mode="markers", line=list(shape="linear"), marker=list(size=15))
fig <- fig %>% add_lines(y=~Publications, name="Publications", mode="markers", line=list(shape="linear"), marker=list(size=15)) 
fig <- fig %>% add_lines(y=~Traits, name="Traits", mode="markers", line=list(shape="linear"), marker=list(size=15))
fig <- fig %>% add_lines(y=~Genes, name="Genes", mode="markers", line=list(shape="linear"), marker=list(size=15))
fig <- fig %>% layout(title = "TIGA Statistics",
                      paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
                      xaxis = list(title = "Release Date (GWAS Catalog and TIGA)",
                                   gridcolor = 'rgb(255,255,255)',
                                   showgrid = TRUE,
                                   showline = FALSE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = FALSE),
                      yaxis = list(title = "Counts",
                                   gridcolor = 'rgb(255,255,255)',
                                   showgrid = TRUE,
                                   showline = FALSE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = FALSE))
ann_pubs_min <- list(
  xref = 'x', yref = 'y',
  x = min(tigastats[, Release]),
  y = tigastats[Release==max(tigastats[, Release]), Publications],
  xanchor = 'center', yanchor = 'top',
  text = sprintf("%d", tigastats[Release==min(tigastats[, Release]), Publications]),
  font = list(size=16, color="blue"), showarrow=F)
ann_pubs_max <- list(
  xref = 'x', yref = 'y',
  x = max(tigastats[, Release]),
  y = tigastats[Release==max(tigastats[, Release]), Publications],
  xanchor = 'center', yanchor = 'top',
  text = sprintf("%d", tigastats[Release==max(tigastats[, Release]), Publications]),
  font = list(size=16, color="blue"), showarrow=F)
fig <- fig %>% layout(annotations = ann_pubs_min)
fig <- fig %>% layout(annotations = ann_pubs_max)
fig

---
title: "GWAX association scoring"
author: "Jeremy Yang"
output:
  html_document
---

```{r echo=FALSE}
knitr::opts_chunk$set(echo=FALSE, message=FALSE)
base::date()
```

# Introduction

Multivariable non-parametric ranking via &mu; scores.

## &mu; scores

 * Non-dominated solutions are cases which are not inferior to any other case at any variable.
 * A &mu; score is defined as the number of lower cases minus the number of higher.
 * The resulting ranking is the useful result, not so much the score itself.
 * Using muStat package.
 * What about ties?


```{r message=FALSE, warning=FALSE}
library(readr)
library(data.table)
library(plotly)
library(muStat)
```

Read gene-trait file.

```{r}
gt <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_traits_g=col_integer(), n_genes_t=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), rcras=col_double()))
setDT(gt)
gt <- gt[!is.na(or_median) & !is.na(gsymb) &!is.na(tdl)]
MIN_ASSN <- 20
traits_df <- gt[, .(.N),  by=c("trait", "trait_uri")]
traits_df <- traits_df[N>=MIN_ASSN]
gt <- gt[trait_uri %in% traits_df$trait_uri]
message(sprintf("Genes in dataset: %d", uniqueN(gt$gsymb)))
message(sprintf("Traits in dataset: %d", uniqueN(gt$trait_uri)))
message(sprintf("G-T associations in dataset: %d", nrow(gt)))
```

Normalize data so bigger is better. I.e. use inverse of n\_traits\_g and n\_genes\_t.

```{r}
gt[, n_traits_g_inv := 1 / n_traits_g]
gt[, n_genes_t_inv := 1 / n_genes_t]
```

We are interested in rankings __for a given trait__. So for each trait,
convert to matrix for muStat::mu.GE().
The (i,j) entry of GE matrix is 1 if \code{x_i >= y_j}, 0 otherwise. The square matrix GE is stored by column in a vector. Thus nrow(GE_matrix) = nrow(x)^2.

```{r}
for (trait_this in unique(gt$trait)) {
  message(sprintf("\"%s\"", trait_this))
  gtmat <- as.matrix(gt[trait==trait_this, .(n_study, n_snp, n_traits_g_inv, n_genes_t_inv, pvalue_mlog_median, or_median, rcras)])
  #rownames(gtmat) <- paste(gt[trait==trait_this, .(gsymb)], sub("^.*/", "", gt[trait==trait_this, .(trait)]), sep=":")
  ge <- mu.GE(gtmat)
  mu_scores <- mu.score(ge)
  #"Error: cannot allocate vector of size 54.2 Gb"!!
  break
}
```


#############################################################################
### Plot a single gene with associated traits.
### Prototype GWA-X in gene-mode.
### Note that gt_stats.csv must exist, produced by gwascat_trait_stats.R.
#############################################################################
library(dplyr, quietly = T)
library(plotly, quietly = T)
library(webshot, quietly=T)

###
gt_stats <- read.csv(file = "data/gt_stats.csv")
###
# ... OR from db:
#library(RMySQL, quietly = T)
#dbcon <- dbConnect(MySQL(), host="localhost", dbname="gwascatalog")
#gt_stats <- dbGetQuery(dbcon, "SELECT * FROM gt_stats")
#ok <- dbDisconnect(dbcon)
###
gt_stats$tdl_color <- NA
gt_stats$tdl_color[gt_stats$tdl == "Tdark"] <- "gray"
gt_stats$tdl_color[gt_stats$tdl == "Tbio"] <- "red"
gt_stats$tdl_color[gt_stats$tdl == "Tchem"] <- "green"
gt_stats$tdl_color[gt_stats$tdl == "Tclin"] <- "blue"
#Color unmapped:
gt_stats$tdl_color[is.na(gt_stats$tdl_color)] <- "#CCCCCC"

#Focus on genes with more data (rows).
t <- table(gt_stats$gsymb)
g <- as.vector(t)
names(g) <- names(t)
g <- sort(g, decreasing = T)
g_top <- g[1:length(g)/20] #top 5%
  
#Plot a random gene:
#gsymb <- sample(names(g_top), 1)

#gsymb <- "APOE" #Apolipoprotein E
#gsymb <- "CDKN2B" #Cyclin-dependent kinase 4 inhibitor B
#gsymb <- "GCKR" #Glucokinase regulatory protein
#gsymb <- "CDKN2A" #Tumor suppressor ARF
#gsymb <- "ABO" #Histo-blood group ABO system transferase
#gsymb <- "FADS1" #Fatty acid desaturase 1
#gsymb <- "FTO" #Alpha-ketoglutarate-dependent dioxygenase FT
#gsymb <- "APOC1" #Apolipoprotein C-I
#gsymb <- "CETP" #Cholesteryl ester transfer protein
#gsymb <- "CSMD1" #CUB and sushi domain-containing protein 1
gsymb <- "CELSR2" #Cadherin EGF LAG seven-pass G-type receptor [IDG2]
#gsymb <- "CAMK1D" #Calcium/calmodulin-dependent protein kinase

g2t <- gt_stats[gt_stats$gsymb == gsymb,]
gname <- g2t$name[g2t$gsymb == gsymb][1]
fam <- g2t$fam[g2t$gsymb == gsymb][1]
print(sprintf("%s: %s", gsymb, gname))
g2t <- g2t[order(g2t$trait),]

for (i in 1:nrow(g2t))
{
  print(sprintf("%2d. n_genes_t = %2d ;n_study = %2d ; n_snp = %2d ; p_median_nlog = %6.1f  (%s) %s\n", i,  
          g2t$n_genes_t[i], g2t$n_study[i], g2t$n_snp[i], g2t$p_median_nlog[i],
          sub('^.*/', '',g2t$trait_uri[i]), g2t$trait[i]))
}

p1 <- plot_ly() %>%
  add_trace(x = 1 / g2t$n_genes_t, y = g2t$p_median_nlog, 
	type = 'scatter', mode = 'markers',
	marker = list(symbol = "circle", color = g2t$tdl_color, size = 4*g2t$n_study),
	text = paste0(sub('^.*/', '', g2t$trait_uri), "<br>(", g2t$trait, ")<br>",
		"n_snp = ",g2t$n_snp, " ; n_study = ", g2t$n_study)
	) %>%
  layout(yaxis = list(title = "median(-log(p))", range=c(0,150)), 
         xaxis = list(type = "log", title = "Specificity (log(1/n_genes_t))"),
        title = paste0("GENE: ", gsymb, " (", fam, ")", " (", gname, ") [",g2t$tdl[1], "]",
                       "<br>(N_trait = ", nrow(g2t), ")"),
        margin = list(t=100,r=50,b=60,l=60),
        font = list(family = "monospace", size=18)) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p1
export(p=p1, file="data/gwax_gene_p1.png")
#
#4 genes via subplot()
gsymbs <- c("APOE", "GCKR", "FTO", "CELSR2")
gnames <- c()
plots <- as.list(rep(NA, 4))
for (i in 1:length(gsymbs))
{
  gsymb <- gsymbs[i]
  g2t <- gt_stats[gt_stats$gsymb == gsymb,]
  gname <- as.character(g2t$name[g2t$gsymb == gsymb][1])
  gnames <- c(gnames, gname)
  fam <- g2t$fam[g2t$gsymb == gsymb][1]
  print(sprintf("%s: %s", gsymb, gname))
  g2t <- g2t[order(g2t$trait),]
  plots[[i]] <- plot_ly(name = gsymb, x = 1 / g2t$n_genes_t, y = g2t$p_median_nlog, 
              type = 'scatter', mode = 'markers',
              marker = list(symbol = "circle", color = g2t$tdl_color, size = 4*g2t$n_study),
              text = paste0(sub('^.*/', '', g2t$trait_uri), "<br>(", g2t$trait, ")<br>",
                            "n_snp = ",g2t$n_snp, " ; n_study = ", g2t$n_study)) %>%
    layout(xaxis = list(type = "log", range = c(log10(.0001),log10(2)), title = "Specificity (log(1/n_genes_t))"), 
           yaxis = list(range = c(0,150), title = "median(-log(p))"))
}

p2 <- subplot(nrows = 2, shareX = T, shareY = T, titleX = T, titleY = T, margin = 0.03,
  plots[[1]],
  plots[[2]],
  plots[[3]],
  plots[[4]]
)  %>%
  layout(title = paste0("GENES: (N = ", length(gsymbs), ")<br>", paste(gsymbs, collapse = ", ")),
         annotations = list(x = c(0, 0.5, 0, 0.5), y = c(0.9, 0.9, 0.4, 0.4), 
                            xanchor = "left",
                            text = paste0(gsymbs, "<br>", gnames), showarrow = F,
                            xref = "paper", yref = "paper"),
         margin = list(t = 60),
         font = list(family = "Arial", size = 14), showlegend = F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p2
export(p=p2, file="data/gwax_gene_p2.png")


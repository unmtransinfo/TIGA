#############################################################################
### Plot single trait with all associated genes.
### Prototype GWA-X in trait-mode.
### Note that gt_stats.csv must exist, produced by gwax_gt_stats.R.
#############################################################################
library(readr, quietly=T)
library(data.table, quietly=T)
library(plotly, quietly=T)
#library(webshot, quietly=T)

###
gt <- read_delim("data/gt_stats.tsv", "\t", col_types=cols(.default=col_character(), 
            n_study=col_integer(), n_snp=col_integer(), n_traits_g=col_integer(), n_genes_t=col_integer(), 
            pvalue_mlog_median=col_double(), or_median=col_double()))
setDT(gt)
###
#
gt$tdl_color <- NA
gt$tdl_color[gt$tdl=="Tdark"] <- "gray"
gt$tdl_color[gt$tdl=="Tbio"] <- "red"
gt$tdl_color[gt$tdl=="Tchem"] <- "green"
gt$tdl_color[gt$tdl=="Tclin"] <- "blue"
#Color unmapped:
gt$tdl_color[is.na(gt$tdl_color)] <- "#CCCCCC"

#Focus on traits with more data (rows).
t <- table(gt$trait_uri)
x <- as.vector(t)
names(x) <- names(t)
x <- sort(x, decreasing=T)
x_top <- x[1:length(x)/20] #top 5%

#Plot a random trait:
#trait_uri <- sample(names(x_top), 1)

trait_id <- "EFO_0001360" #T2DM
#trait_id <- "EFO_0000249" #Alzheimer
#trait_id <- "EFO_0000289" #bipolar disorder
#trait_id <- "EFO_0000249" #Alzheimers disease
#trait_id <- "EFO_0000305" #breast carcinoma
#trait_id <- "EFO_0000270" #asthma
#trait_id <- "EFO_0000692" #schizophrenia
#trait_id <- "EFO_0004340" #body mass index
#trait_id <- "EFO_0005842" #colorectal cancer
#trait_id <- "EFO_0001663" #prostate carcinoma
#trait_id <- "EFO_0004612" #high density lipoprotein cholesterol
#trait_id <- "EFO_0004530" #triglyceride measurement
#trait_id <- "EFO_0000685" #rheumatoid arthritis
#trait_id <- "EFO_0004611" #low density lipoprotein cholesterol
#trait_id <- "EFO_0003761" #unipolar depression
#trait_id <- "EFO_0002508" #Parkinsons disease
#trait_id <- "EFO_0001645" #coronary heart disease
#trait_id <- "EFO_0003885" #multiple sclerosis
#trait_id <- "EFO_0004339" #body height
#trait_id <- "EFO_0000384" #Crohns disease

trait_uri <- paste0("http://www.ebi.ac.uk/efo/", trait_id)

trait <- gt$trait[gt$trait_uri==trait_uri][1]
print(sprintf("%s: %s", trait_uri, trait))
t2g <- gt[trait_uri==trait_uri]
t2g <- t2g[!is.na(gsymb)]

t2g <- t2g[order(gsymb)]

#for (i in 1:nrow(t2g))
#{
#  print(sprintf("%2d. n_traits_g=%2d ; n_study=%2d ; n_snp=%2d ; or_median=%3.1f  (%s) %s\n", i,  
#          t2g$n_traits_g[i], t2g$n_study[i], t2g$n_snp[i], t2g$or_median[i], t2g$gsymb[i], t2g$name[i]))
#}

t2g[, specificity := (1 / n_traits_g)]
t2g[, tdl := as.character(tdl)]
t2g[is.na(tdl), tdl := "NA"]
t2g[, markersize := 20*log(n_study+1)]
t2g[, text := paste0(gsymb, "<br>", name, "<br>", "Family: ", fam, "<br>", "TDL: ", tdl, "<br>", "n_traits_this_gene=", n_traits_g, " ; n_snp=", n_snp, " ; n_study=", n_study)]
#
p1 <- plot_ly(type='scatter', mode='markers') %>%
  add_markers(x=t2g[tdl=="Tclin", specificity], y=t2g[tdl=="Tclin", or_median], name="Tclin",
    marker=list(symbol="circle", color=t2g[tdl=="Tclin", tdl_color], size=t2g[tdl=="Tclin", markersize]), text=t2g[tdl=="Tclin", text]
  ) %>%
  add_markers(x=t2g[tdl=="Tchem", specificity], y=t2g[tdl=="Tchem", or_median], name="Tchem",
    marker=list(symbol="circle", color=t2g[tdl=="Tchem", tdl_color], size=t2g[tdl=="Tchem", markersize]), text=t2g[tdl=="Tchem", text]
  ) %>%
  add_markers(x=t2g[tdl=="Tbio", specificity], y=t2g[tdl=="Tbio", or_median], name="Tbio",
    marker=list(symbol="circle", color=t2g[tdl=="Tbio", tdl_color], size=t2g[tdl=="Tbio", markersize]), text=t2g[tdl=="Tbio", text]
  ) %>%
  add_markers(x=t2g[tdl=="Tdark", specificity], y=t2g[tdl=="Tdark", or_median], name="Tdark",
    marker=list(symbol="circle", color=t2g[tdl=="Tdark", tdl_color], size=t2g[tdl=="Tdark", markersize]), text=t2g[tdl=="Tdark", text]
  ) %>%
  add_markers(x=t2g[tdl=="NA", specificity], y=t2g[tdl=="NA", or_median], name="NA",
    marker=list(symbol="circle", color=t2g[tdl=="NA", tdl_color], size=t2g[tdl=="NA", markersize]), text=t2g[tdl=="NA", text]
  ) %>%
  layout(xaxis=list(type="log", title="Specificity (1/n_trait)"), 
         yaxis=list(type="log", title="Effect (median(OR))"), 
        title=paste0(trait, "<br>(GWAS Catalog TRAIT: ", sub('^.*/', '', trait_uri), ")"),
        margin=list(t=100,r=50,b=60,l=60), legend=list(x=.95, y=.9), showlegend=T, font=list(family="monospace", size=18)) %>%
  add_annotations(text=paste0("(N_gene=", nrow(t2g), ")"), showarrow=F, x=0.1, y=1.0, xref="paper", yref="paper") %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p1
#export(p=p1, file="data/gwax_trait_p1.png")
#
#4 traits via subplot()
trait_ids <- c("EFO_0001360", "EFO_0000249", "EFO_0000289", "EFO_0005842")
trait_uris <- paste0("http://www.ebi.ac.uk/efo/", trait_ids)
traits <- c()
plots <- list()
for (i in 1:length(trait_uris))
{
  trait_uri <- trait_uris[i]
  t2g <- gt[trait_uri==trait_uri]
  trait <- as.character(t2g[t2g$trait_uri==trait_uri, trait][1])
  traits <- c(traits, trait)
  print(sprintf("%s: %s", trait_uri, trait))

  plots[[i]] <- plot_ly(name=trait_uri, x=1/t2g$n_traits_g, y=t2g$or_median, 
	type='scatter', mode='markers',
	marker=list(symbol="circle", color=t2g$tdl_color, size=4*t2g$n_study),
	text=paste0(t2g$gsymb, "<br>", t2g$name, "<br>",
	              "Family: ", t2g$fam, "<br>",
	              "TDL: ", t2g$tdl, "<br>",
	              "n_traits_this_gene=", t2g$n_traits_g, " ; n_snp=",t2g$n_snp, " ; n_study=", t2g$n_study)) %>%
    layout(xaxis=list(type="log", title="Specificity (log(1/n_traits_g))"), 
           yaxis=list(type="log", title="Effect (median(OR))"))
}

p2 <- subplot(plots, nrows=2, shareX=T, shareY=T, titleX=T, titleY=T, margin=0.03)  %>%
  layout(title=paste0("GWAS Catalog TRAITS:<br>",
        paste(traits[1:2], collapse=", "), "<br>", paste(traits[3:4], collapse=", ")),
	annotations=list(x=c(0, 0.5, 0, 0.5), y=c(0.9, 0.9, 0.4, 0.4), 
	text=paste0(trait_ids, "<br>", traits), showarrow=F, xanchor="left", xref="paper", yref="paper"),
	margin=list(t=60), font=list(family="Arial", size=14),  showlegend=F) %>%
  add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper")
#
p2
#export(p=p2, file="data/gwax_trait_p2.png")
#

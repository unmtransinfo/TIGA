########################################################################################
### TIGA: Target Illumination by GWAS Analytics
### gt = gene-trait data
### Dataset gt_stats.tsv.gz from tiga_gt_stats.R.
### Jeremy Yang
########################################################################################
### Requires dqshiny dev version late 2019, via https://github.com/daqana/dqshiny
### remotes::install_github("daqana/dqshiny")
########################################################################################
### Autosuggest is dependent on dqshiny and shiny versions and is fussy. shiny 1.4.0 is
### ok (but possibly not 1.4.0.2) with dqshiny 0.0.3.9000 and shinyBS 0.61.
########################################################################################
library(readr)
library(data.table)
library(igraph, quietly=T)
library(shiny, quietly=T)
library(shinyBS, quietly=T)
library(dqshiny, quietly=T) #https://github.com/daqana/dqshiny
library(DT, quietly=T)
library(plotly, quietly=T)
#
pkgs <- names(sessionInfo()$otherPkgs)
pkgVerTxt <- paste(sprintf("%s %s", pkgs, sapply(pkgs, function(p){paste(packageVersion(p), collapse=".")})), collapse="; ")
message(pkgVerTxt)
###
efoId2Uri <- function(efoId) { #(non-vector)
  if (is.null(efoId)) { return(NA) }
  if (grepl("^EFO_", efoId)) {
    return(sprintf("http://www.ebi.ac.uk/efo/%s", efoId))
  } else if (grepl("^Orphanet_", efoId)) {
    return(sprintf("http://www.orpha.net/ORDO/%s", efoId))
  } else if (grepl("(^HP_|^GO_|^CHEBI_|^UBERON_|^NCBITaxon_|^MONDO_|^UO_|^CL_|^PATO_|^HANCESTRO_)", efoId)) {
    return(sprintf("http://purl.obolibrary.org/obo/%s", efoId))
  } else {
    return(NA)
  }
}
# Return all subclass efoIds for input efoId, including self.
efoId2Subclasses <- function(G, id_this) {
  if (!(id_this %in% V(G)$efoId)) { return(NULL) }
  v_this <- V(G)[V(G)$efoId == id_this]
  bfs_this <- igraph::bfs(G, v_this, neimode="out", unreachable=F)
  subG <- induced_subgraph(G, bfs_this$order[1:sum(!is.na(bfs_this$order))])
  return(V(subG)$efoId)
}
########################################################################################
APPNAME <- "TIGA"
MIN_ASSN <- 1
#
t0 <- proc.time()
DEBUG <- F
if (!file.exists("tiga.Rdata") | DEBUG) {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  gt <- read_delim("gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(), geneNtrait=col_integer(), geneNstudy=col_integer(),	traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(), geneMuScore=col_double(), geneMuRank=col_integer(),	traitMuScore=col_double(), traitMuRank=col_integer()))
  setDT(gt)
  setnames(gt, old=c("geneIdgTdl"), new=c("TDL"))
  #
  gt_prov <- read_delim("gt_provenance.tsv.gz", "\t", col_types=cols(.default=col_character()))
  setDT(gt_prov)
  #
  filtered_studies <- read_delim("filtered_studies.tsv.gz", "\t")
  setDT(filtered_studies)
  filtered_studies_trait <- read_delim("filtered_studies_trait.tsv.gz", "\t")
  setDT(filtered_studies_trait)
  filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_trait))
  filtered_studies[, type := "study"]
  filtered_traits <- read_delim("filtered_traits.tsv.gz", "\t")
  setDT(filtered_traits)
  filtered_traits[, type := "trait"]
  filtered_genes <- read_delim("filtered_genes.tsv.gz", "\t")
  setDT(filtered_genes)
  filtered_genes[, type := "gene"]
  filtered_gene_menu <- paste0("gene:", filtered_genes$ensemblId) #named vector
  names(filtered_gene_menu) <- sprintf("%s:%s", filtered_genes$ensemblSymb, filtered_genes$geneName)
  filtered <- rbindlist(list(filtered_studies[, .(type, id=STUDY_ACCESSION, reason)], filtered_traits[, .(type, id=TRAIT_URI, reason)], filtered_genes[, .(type, id=ensemblId, reason)]))
  #
  trait_table <- gt[, .(trait=first(trait), N_study=first(traitNstudy), N_gene=uniqueN(ensemblId)),  by="efoId"]
  message(sprintf("Traits with N_assn<%d: %d", MIN_ASSN, trait_table[N_gene<MIN_ASSN, uniqueN(efoId)]))
  trait_table <- trait_table[N_gene>=MIN_ASSN]
  trait_table[, trait_uri := sapply(efoId, efoId2Uri)]
  trait_table <- trait_table[, .(trait_uri, efoId, trait, N_study, N_gene)][order(trait)]
  trait_menu <- paste0("trait:", trait_table$efoId) #named vector
  names(trait_menu) <- trait_table$trait
  #
  gene_table <- gt[, .(geneSymbol, geneName, geneFamily, TDL, N_study = geneNstudy, N_trait = uniqueN(efoId), filtered=F),  by=c("ensemblId")]
  gene_table <- rbindlist(list(gene_table, filtered_genes[, .(ensemblId, geneSymbol=ensemblSymb, geneName, geneFamily, TDL, N_study=0, N_trait=0, filtered=T)]))
  gene_table <- unique(gene_table)[order(geneSymbol)]
  #
  gene_menu <- paste0("gene:", gene_table$ensemblId) #named vector
  names(gene_menu) <- sprintf("%s:%s", gene_table$geneSymbol, gene_table$geneName)
  #
  qry_menu <- as.list(c(trait_menu, gene_menu, filtered_gene_menu))
  trait_menu <- as.list(trait_menu)
  gene_menu <- as.list(gene_menu)
  #
  study_table <- read_delim("gwascat_gwas.tsv.gz", "\t", col_types = cols(.default = col_character(), DATE=col_date(), DATE_ADDED_TO_CATALOG=col_date()))
  setDT(study_table)
  study_table <- study_table[, .(STUDY_ACCESSION, STUDY, PUBMEDID, DATE_PUBLISHED = DATE, DATE_ADDED_TO_CATALOG)][order(DATE_PUBLISHED)]
  #
  system("if [ -f \"efo_graph.graphml.gz\" ]; then gunzip -f efo_graph.graphml.gz ; fi")
  efoGraph <- read_graph("efo_graph.graphml", format="graphml")
  #
  save(gt, trait_table, gene_table, trait_menu, gene_menu, filtered_gene_menu, qry_menu, gt_prov, filtered, study_table, efoGraph, file="tiga.Rdata")
} else {
  message(sprintf("Loading tiga.Rdata..."))
  load("tiga.Rdata")
}
#
t_elapsed <- (proc.time()-t0)[3]
#
message(sprintf("Gene count, IDs: %d; symbols: %d", uniqueN(gt$ensemblId), uniqueN(gt$geneSymbol)))
message(sprintf("Trait count (total): %d", uniqueN(gt$efoId)))
#trait_table[['ontology']] <- as.factor(sub("_.*$","", trait_table$efoId))
#trait_counts <- trait_table[, .N, by="ontology"][order(-N)]
#message(sprintf("traits (%10s): %4d / %4d (%4.1f%%)\n", trait_counts$ontology, trait_counts$N, sum(trait_counts$N), 100*trait_counts$N/sum(trait_counts$N)))
message(sprintf("Trait count (filtered; n_assn<%d): %d", MIN_ASSN, uniqueN(gt$efoId)-length(trait_menu)))
message(sprintf("Trait count (menu; n_assn>=%d): %d", MIN_ASSN, length(trait_menu)))
#message(sprintf("DEBUG: COUNT or_median: %d", sum(!is.na(gt$or_median))))
#message(sprintf("DEBUG: COUNT pvalue_mlog_median: %d", sum(!is.na(gt$pvalue_mlog_median))))
#message(sprintf("DEBUG: COUNT rcras: %d", sum(!is.na(gt$rcras))))
#
message(sprintf("Provenance: %d PUBMEDIDS and %d studies for %d gene-trait pairs.", gt_prov[, uniqueN(PUBMEDID)], gt_prov[, uniqueN(STUDY_ACCESSION)], nrow(unique(gt_prov[, .(ensemblId, TRAIT_URI)]))))
flt <- data.table(table(filtered$type, filtered$reason))[N>0][order(V1, V2)]
message(sprintf("FILTERED entities: (%5s) %5d REASON: %s\n", flt$V1, flt$N, flt$V2))
#
message(sprintf("Graph \"%s\": vertices: %d; edges: %d", graph_attr(efoGraph, "name"), vcount(efoGraph), ecount(efoGraph)))
#
dbHtm <- sprintf("<B>Dataset:</B> genes: %d; traits: %d ; studies: %d; publications: %d (t_load: %.1fs)", 
                 uniqueN(gt$ensemblId), uniqueN(gt$efoId), uniqueN(gt_prov$STUDY_ACCESSION), uniqueN(gt_prov$PUBMEDID), t_elapsed)
#
###
idgfams <- c("GPCR", "Kinase", "IC", "NR", "Other")
axes <- c("Effect", "Evidence")
#
#############################################################################
HelpHtm <- function() {
  htm <- ("<P><B>TIGA</B>, Target Illumination by GWAS Analytics, is designed to facilitate drug target illumination by 
scoring and ranking of protein-coding genes associated with traits from genome-wide association studies
(GWAS). Similarly, <B>TIGA</B> can score and rank traits with the same gene-trait association metrics. 
Rather than a comprehensive analysis of GWAS for all biological implications and insights, this
more focused application provides a rational method by which GWAS findings can be 
aggregated and filtered for applicable, actionable intelligence, 
evidence usable by drug discovery scientists to enrich prioritization of target hypotheses. 
Data from the <A HREF=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">NHGRI-EBI GWAS Catalog</A>.
<UL>
<LI>Traits are mapped to EFO, HPO, Orphanet, PATO or GO, with ~90%% mapped to EFO.
<LI>Mapped genes via Ensembl pipeline as per GWAS Catalog documentation. Reported genes ignored for consistency and accountable confidence assessment in this app and downstream.
<LI>In this version, effect size measure (1) odds ratio (OR) or (2) BETA required.
<LI>Due to difficulties parsing and harmonizing beta units, simple count N_beta is employed.
</UL>
<B>Datatypes:</B>
<UL>
  <LI><B>pVal_mLog<SUP>*</SUP></B>: median(-Log(pValue)) supporting trait-gene association.
  <LI><B>OR<SUP>*</SUP></B>: median(odds ratio, inverted if &lt;1) supporting trait-gene association (computed as one if missing).
  <LI><B>N_beta<SUP>*</SUP></B>: simple count of beta values with 95% confidence intervals supporting trait-gene association.
  <LI><B>N_snp<SUP>*</SUP></B>: SNPs involved with trait-gene association.
  <LI><B>N_snpw<SUP>*</SUP></B>: N_snp weighted by distance inverse exponential.
  <LI><B>N_study<SUP>*</SUP></B>: studies supporting trait-gene association.
  <LI><B>study_N<SUP>*</SUP></B>: mean(SAMPLE_SIZE) supporting trait-gene association.
  <LI><B>RCRAS<SUP>*</SUP></B>: Relative Citation Ratio (RCR) Aggregated Score (iCite-RCR-based)
  <LI><B>geneNtrait<SUP>**</SUP></B>: total traits associated with gene.
  <LI><B>traitNgene<SUP>**</SUP></B>: total genes associated with trait.
  <LI><B>muScore</B>: For a given query trait, genes are scored and ranked, or for a given gene, traits are scored and ranked. 
</UL>
<SUP>*</SUP>Variable used in <B>muScore</B>.
<SUP>**</SUP>Variable inverse used in <B>muScore</B>.
<BR/>
Hits are filtered and ranked based on non-parametric multivariate &mu; scores
(Wittkowski, 2008).
<BR/>
<B>UI:</B>
Scatterplot axes are Effect (OR) vs. Evidence as measured by <B>muScore</B>.
<UL>
<LI>Genes markers colored by TDL, and may be be hidden via legend.
<LI>Trait markers colored by EFO top-level class, and may be hidden via legend (TO DO).
<LI>Plot markers may be sized by <B>N_study</B> or <B>RCRAS</B>.
<LI>Note that this app will accept query parameter <B>trait</B> via URL, e.g.
<B><TT>?trait=EFO_0000341</TT></B>.
</UL>
<B>More documentation:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\">GWAS Catalog</a>
<LI><a href=\"https://www.ebi.ac.uk/gwas/docs/fileheaders\">GWAS Catalog data dictionary</a>
<LI><a href=\"https://www.ebi.ac.uk/efo/\">Experimental Factor Ontology (EFO)</a>
</UL>
<B>Issues:</B>
<UL>
<LI>Query traits with apostrophes must be typed in full. Autocomplete autosuggests but does not complete.
<LI>Plot markers may obscure others.
</UL>
<B>Authors:</B> Jeremy Yang<SUP>1</SUP>, Stephen Mathias<SUP>1</SUP>, Cristian
Bologa<SUP>1</SUP>, Lars Juhl Jensen<SUP>2</SUP>, Christophe Lambert<SUP>1</SUP>, David Wild<SUP>3</SUP> and Tudor
Oprea<SUP>1</SUP>.<BR/>
<I><SUP>1</SUP>University of New Mexico, Translational Informatics Division, Dept. of
Internal Medicine; <SUP>2</SUP>Novo Nordisk Center for Protein Research, Copenhagen,
Denmark; <SUP>3</SUP>Indiana University, School of Informatics, Computing and Engineering, Integrative Data Science Lab.</I>
<BR/>
<B>Feedback welcome</B> to corresponding author  
<a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>")
  htm <- paste(htm, sprintf("<hr>\nBuilt with:\n<pre>%s; %s</pre>", R.version.string, pkgVerTxt), sep="\n")
  return(htm)
}

##########################################################################################
ui <- fluidPage(
  tags$style(".green_class {color:#00ff00} .blue_class {color:#0000ff} .red_class {color:#ff0000} .black_class {color:black}"),
  titlePanel(h2("IDG", tags$img(height="50", valign="bottom", src="IDG_logo_only.png"), sprintf("%s: Target Illumination by GWAS Analytics", APPNAME),tags$img(height="40", valign="bottom", src="GWAS_Catalog_logo.png"), span(style="font-size:18px", "GWAS Catalog based app (BETA)")), windowTitle=APPNAME),
  fluidRow(
    column(3, 
      wellPanel(

	dqshiny::autocomplete_input("traitQry", div("Trait",
		actionButton("randTraitQry", tags$img(height="28", valign="bottom", src="dice.png"), style='padding:0px; background-color:#DDDDDD'),
		actionButton("goReset", tags$img(height="28", valign="bottom", src="refresh_icon.png"), style='padding:0px;background-color:#DDDDDD')),
		options=trait_menu, max_options=1000, placeholder="Query trait..."),

	dqshiny::autocomplete_input("geneQry", div("Gene"), options=as.list(c(gene_menu, filtered_gene_menu)), max_options=10000, placeholder="Query gene..."),

        sliderInput("maxHits", "MaxHits", 25, 200, 50, step=25),
	checkboxGroupInput("logaxes", "LogAxes", choices=axes, selected=NULL, inline=T),
        radioButtons("markerSizeBy", "MarkerSizeBy", choiceNames=c("N_study", "RCRAS", "None"), choiceValues=c("n_study", "rcras", NA), selected="n_study", inline=T)
      ),
	wellPanel(htmlOutput(outputId="logHtm")),
	wellPanel(htmlOutput(outputId="resultHtm"))
	),
    column(9,
	tabsetPanel(id="tabset", type="tabs",
		tabPanel(value="plot", title=textOutput("plotTabTxt"), plotlyOutput("tigaPlot", height = "500px")),
		tabPanel(id="hits", title=textOutput("hitsTabTxt"), DT::dataTableOutput("hitrows"), br(), downloadButton("hits_file", label="Download Hits")),
		tabPanel(value="association_detail", title="Detail",
			htmlOutput("association_detail"),
			DT::dataTableOutput("association_detail_studies"),
			DT::dataTableOutput("association_detail_publications")
			),
		tabPanel(value="traits", title="Traits (all)", DT::dataTableOutput("traits")),
		tabPanel(value="genes", title="Genes (all)", DT::dataTableOutput("genes")),
		tabPanel(value="studies", title="Studies (all)", DT::dataTableOutput("studies")),
		tabPanel(value="download", title="Download",
			h1("Downloads"),
			p(em("NOTE: gene-trait file contains all data in separate genes and traits files.")),
			p(downloadButton("gt_file", label="Gene-Trait Associations (all)"), textOutput("gtFileInfoTxt")),
			p(downloadButton("traits_file", label="Traits (all)"), textOutput("traitFileInfoTxt")),
			p(downloadButton("genes_file", label="Genes (all)"), textOutput("geneFileInfoTxt"))
		),
		tabPanel(value="help", title="Help", htmlOutput("helpHtm"))
	))),
  hr(),
  fluidRow(
    column(12, tags$em(strong(sprintf("%s", APPNAME)), " web app from ", 
        tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
        " and ",
        tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
        " built upon ",
        tags$a(href="https://www.ebi.ac.uk/gwas/", target="_blank", span("GWAS Catalog", tags$img(id="gwas_catalog_logo", height="50", valign="bottom", src="GWAS_Catalog_logo.png"))),
        " and ",
        tags$a(href="https://www.ebi.ac.uk/efo/", target="_blank", span("EFO", tags$img(id="efo_logo", height="50", valign="bottom", src="EFO_logo.png")))
        ))),
  bsTooltip("randTraitQry", "Random query trait or gene", "right"),
  bsTooltip("jitter", "Jitter can un-hide plot markers.", "right"),
  bsTooltip("goReset", "Reset.", "right"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

##########################################################################################
server <- function(input, output, session) {
  observeEvent(input$showHelp, {
    showModal(modalDialog(easyClose=T, footer=tagList(modalButton("Dismiss")),
	title=HTML(sprintf("<H2>%s Help</H2>", APPNAME)),
	HTML(HelpHtm())))
  })
  
  output$helpHtm <- reactive({ paste(sprintf("<H2>%s Help</H2>", APPNAME), HelpHtm()) })
  #output$downloadHtm <- reactive({ paste(sprintf("<H2>%s download</H2>", APPNAME), DownloadHtm()) })

  Sys.sleep(1) #Needed?
  traitQryRand_count <- 0 # initialize once per session
  i_query <- 0 # initialize once per session
  
  # ?trait=EFO_0000341
  # ?trait=EFO_1000654&gene=ENSG00000094914
  httpQstr <- reactive({
    qStr <- getQueryString(session) #named list
    if (length(qStr)>0) {
      for (key in names(qStr)) {
        message(sprintf("DEBUG: qStr[[\"%s\"]]=\"%s\"", key, qStr[[key]]))
      }
    }
    return(qStr)
  })

  urlBase <- reactive({
    sprintf("%s//%s:%s%s",
      session$clientData$url_protocol, session$clientData$url_hostname,
      session$clientData$url_port, session$clientData$url_pathname)
  })
  urlText <- reactive({
    sprintf("%s%s", urlBase(), session$clientData$url_search)
  })
 
  qryType <- reactive({
    if (input$traitQry!="" & input$geneQry!="") { return("genetrait") }
    else if (input$traitQry!="" & input$geneQry=="") { return("trait") }
    else if (input$traitQry=="" & input$geneQry!="") { return("gene") }
    else { return(NA) }
  })
  hitType <- reactive({
    ifelse(is.na(qryType()), NA,
    ifelse(qryType()=="genetrait", "genetrait",
    ifelse(qryType()=="trait", "gene",
    ifelse(qryType()=="gene", "trait", NA))))
  })
  
  observeEvent(input$goReset, {
    session$reload()
  })
  
  efoId2Name <- function(efoId_this) {
    name <- trait_table[efoId==efoId_this, first(trait)]
    #message(sprintf("DEBUG: efoId: \"%s\"; trait name: \"%s\"", efoId_this, name))
    return(name)
  }
  ensemblId2Symbol <- function(ensemblId_this) {
    gsymb <- gene_table[ensemblId==ensemblId_this, first(geneSymbol)]
    #message(sprintf("DEBUG: ensemblId: \"%s\"; gene symbol: %s", ensemblId_this, gsymb))
    return(gsymb)
  }
  ensemblId2Name <- function(ensemblId_this) {
    name <- gene_table[ensemblId==ensemblId_this, first(geneName)]
    #message(sprintf("DEBUG: ensemblId: \"%s\"; gene name: %s", ensemblId_this, name))
    return(name)
  }
  
  AssociationDetailHtm <- function(efoId_this, ensemblId_this) {
    htm <- h2("Gene-trait association detail")
    message(sprintf("DEBUG: AssociationDetailHtm: %s; %s", efoId_this, ensemblId_this))
    traitName <- efoId2Name(efoId_this)
    geneName <- ensemblId2Name(ensemblId_this)
    message(sprintf("DEBUG: AssociationDetailHtm TRAIT: %s; GENE: %s", traitName, geneName))
    #htm <- sprintf("TRAIT: %s; GENE: %s", traitName, geneName)
    htm <- c(htm, h3("TRAIT: ", traitName, "; GENE: ", geneName))
    return(htm)
  }

  processUrlParams <- function() {
    message(sprintf("DEBUG: url: \"%s\"", urlText()))
    qStr <- httpQstr()
    if ("trait" %in% names(qStr)) {
      dqshiny::update_autocomplete_input(session, "traitQry", value=efoId2Name(qStr[["trait"]]))
    }
    if ("gene" %in% names(qStr)) {
      dqshiny::update_autocomplete_input(session, "geneQry", value=sprintf("%s:%s", ensemblId2Symbol(qStr[["gene"]]), ensemblId2Name(qStr[["gene"]])))
    }
  }

  traitQryId <- reactive({
    if (i_query==0 & length(names(httpQstr()))>0) { processUrlParams() } #1st query may be via URL http param.
    message(sprintf("DEBUG: input$traitQry: \"%s\"", input$traitQry))
    if (input$randTraitQry>traitQryRand_count) {
      traitQryRand_count <<- input$randTraitQry # Must assign to up-scoped variable.
      traitQryRand_new <- sample(trait_menu, 1)
      dqshiny::update_autocomplete_input(session, "traitQry", value=as.character(names(traitQryRand_new)))
    }
    i_query <<- i_query + 1  # Must assign to up-scoped variable.
    updateTabsetPanel(session, "tabset", selected="plot")
    if (grepl("^\\s*$", input$traitQry)) { return("") }
    return(sub("^.*:", "", input$traitQry))
  })

  geneQryId <- reactive({
    message(sprintf("DEBUG: input$geneQry: \"%s\"", input$geneQry))
    if (grepl("^\\s*$", input$geneQry)) { return("") }
    return(sub("^.*:", "", input$geneQry))
  })

  traitQryName <- reactive({
    if (traitQryId()=="") { return(NULL) }
    if (!(traitQryId() %in% trait_table$efoId)) { return(NULL) }
    return(trait_table[efoId==traitQryId(), trait])
  })
  geneQryName <- reactive({
    if (geneQryId()=="") { return(NULL) }
    if (!(geneQryId() %in% gene_table$ensemblId)) { return(NULL) }
    return(sprintf("%s:%s", gene_table[ensemblId==geneQryId(), geneSymbol], gene_table[ensemblId==geneQryId(), geneName]))
  })


  Hits <- reactive({
    message(sprintf("DEBUG: qryType(): %s", qryType()))
    if (traitQryId()=="" & geneQryId()=="") { return(NULL) }
    if (hitType()=="trait") {
      gt_this <- gt[ensemblId==geneQryId()]
      if (nrow(gt_this)==0) { return(NULL) }
      gt_this <- gt_this[, .(efoId, trait, n_study, study_N_mean, n_snp, n_snpw, traitNgene, pvalue_mlog_median, or_median, betaN, rcras, traitMuScore, traitMuRank)]
      setnames(gt_this, old=c("traitMuScore", "traitMuRank"), new=c("muScore", "muRank"))
      gt_this[, ok := (is.na(muRank) | as.logical(muRank<=input$maxHits))] #1-hit situation problematic.
      setorder(gt_this, muRank)
    } else if (hitType()=="gene") {
      gt_this <- gt[efoId==traitQryId()]
      if (nrow(gt_this)==0) { return(NULL) }
      gt_this$TDL <- factor(gt_this$TDL, levels=c("NA", "Tdark", "Tbio", "Tchem", "Tclin"), ordered=T)
      gt_this <- gt_this[, .(ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, geneNtrait, pvalue_mlog_median, or_median, betaN,rcras, geneMuScore, geneMuRank)]
      setnames(gt_this, old=c("betaN", "geneMuScore", "geneMuRank"), new=c("n_beta", "muScore", "muRank"))
      gt_this[, ok := (is.na(muRank) | as.logical(muRank<=input$maxHits))] #1-hit situation problematic.
      setorder(gt_this, muRank)
    } else { #hitType=="genetrait"
      gt_this <- gt[ensemblId==geneQryId() & efoId==traitQryId()]
      gt_this <- gt_this[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, traitNgene, pvalue_mlog_median, or_median, betaN, rcras)]
      gt_this[, ok := T]
    }
    message(sprintf("DEBUG: nrow(gt_this)=%d", nrow(gt_this)))
    return(gt_this)
  })

  output$hitCount <- renderText({ as.character(nrow(Hits())) })
  
  output$plotTabTxt <- renderText({ ifelse(is.null(Hits()), "Plot", sprintf("Plot (%ss)", hitType())) })
  output$hitsTabTxt <- renderText({ ifelse(!is.null(Hits()), sprintf("Hits (%d %ss)", nrow(Hits()), hitType()), "Hits (0)") })

  output$traitFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(trait_table), ncol(trait_table)) })
  output$geneFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(gene_table), ncol(gene_table)) })
  output$gtFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(gt), ncol(gt)) })
  output$association_detail <- reactive({ AssociationDetailHtm(traitQryId(), geneQryId()) })

  #Hits table has links to tiga:trait+gene
  HitsWithHtm <- reactive({
    hh <- data.table(Hits()) #copy
    if (hitType() == "trait") {
      hh <- hh[, efoId := sprintf("<a href=\"%s\" target=\"_blank\">EFO:%s</a><br><a href=\"%s?trait=%s&gene=%s\">TIGA:%s+%s</a>", sapply(efoId, efoId2Uri), efoId, urlBase(), efoId, geneQryId(), efoId, geneQryId())]
    } else if (hitType() == "gene") {
      hh <- hh[, geneSymbol := sprintf("<a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\">IDG:%s</a><br><a href=\"%s?trait=%s&gene=%s\">TIGA:%s+%s</a>", geneSymbol, geneSymbol, urlBase(), traitQryId(), ensemblId, traitQryId(), geneSymbol)]
    } else if (hitType() == "genetrait") {
      hh <- hh[, efoId := sprintf("<a href=\"%s\" target=\"_blank\">EFO:%s</a>", sapply(efoId, efoId2Uri), efoId)]
      hh <- hh[, geneSymbol := sprintf("<a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\">IDG:%s</a>", geneSymbol, geneSymbol)]
    }
    return(hh)
  })

  output$resultHtm <- reactive({
    message(sprintf("TraitQuery: \"%s\" (%s)", traitQryName(), traitQryId()))
    message(sprintf("GeneQuery: \"%s\" (%s)", geneQryName(), geneQryId()))
    if (traitQryId()=="" & geneQryId()=="") { return("No query. Search? Browse? Or roll dice?") }
    htm <- sprintf("<B>Results:</B>")
    if (qryType() %in% c("trait", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", traitQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", efoId2Uri(traitQryId()), traitQryId()))
    }
    if (qryType() %in% c("gene", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", geneQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"https://pharos.nih.gov/targets/%s\">%s</a>)", geneQryId(), geneQryId()))
    }
    if (!is.null(Hits())) {
      htm <- paste0(htm, sprintf("; N_%s: %d", hitType(), nrow(Hits())))
      if (!is.null(Hits()[["or_median"]]) & !is.null(Hits()[["n_beta"]]))
        htm <- paste0(htm, sprintf("; ORs: %d; betas: %d", Hits()[!is.na(or_median), .N], Hits()[!is.na(n_beta), .N]))
    } else if (geneQryId() %in% filtered$id) {
      htm <- paste0(htm, sprintf("; %s <B>%s: %s</B> filtered; reason: %s.", hitType(), geneQryId(), geneQryName(), filtered[id==geneQryId(), reason]))
    } else {
      htm <- paste0(htm, sprintf("; No %ss found.", hitType()))
    }
    return(htm)
  })

  output$logHtm <- reactive({
    htm <- dbHtm
    return(htm)
  })
  
  #ERRORS here to be fixed.
  markerSize <- reactive({
    if (input$markerSizeBy=="n_study") {
      message(sprintf("DEBUG: markerSizeBy: %s", input$markerSizeBy))
      size <- 10*Hits()[(ok), n_study]
    } else if (input$markerSizeBy=="rcras") {
      message(sprintf("DEBUG: markerSizeBy: %s", input$markerSizeBy))
      size <- 10*Hits()[(ok), rcras]
    } else { #NA
      message(sprintf("DEBUG: markerSizeBy: %s", as.character(input$markerSizeBy)))
      size <- rep(10, nrow(Hits()[(ok)]))
    }
    message(sprintf("DEBUG: nrow(Hits()): %d", nrow(Hits()[(ok)])))
    message(sprintf("DEBUG: length(markerSize): %d", length(size)))
    size <- pmax(size, rep(10, nrow(Hits()[(ok)]))) #min
    size <- pmin(size, rep(80, nrow(Hits()[(ok)]))) #max
    message(sprintf("DEBUG: length(markerSize): %d", length(size)))
    return(size)
  })

  markerTextGenes <- reactive({
    text=paste0(
    paste0("<b>", Hits()[(ok)]$geneSymbol, "</b> (", Hits()[(ok)]$ensemblId, ")"), 
    paste0("<br><b>", Hits()[(ok)]$geneName, "</b>"),
    paste0("<br>Fam:", Hits()[(ok)]$geneFamily),
    paste0(", TDL:", Hits()[(ok)]$TDL),
    paste0("; N_trait = ", Hits()[(ok)]$geneNtrait),
    ";<br>muScore = ", Hits()[(ok)]$muScore,
    "; muRank = ", Hits()[(ok)]$muRank,
    ";<br>N_study = ", Hits()[(ok)]$n_study,
    "; study_N = ", Hits()[(ok)]$study_N_mean, 
    "; N_snp = ", Hits()[(ok)]$n_snp,
    ";<br>N_snpw = ", Hits()[(ok)]$n_snpw,
    "; OR = ", round(Hits()[(ok)]$or_median, digits=2), 
    "; N_beta = ", Hits()[(ok)]$n_beta,
    "; pVal = ", sprintf("%.2g", 10^(-Hits()[(ok)]$pvalue_mlog_median)),
    "; RCRAS = ", round(Hits()[(ok)]$rcras, digits=2))
    return(text)
  })
  markerTextTraits <- reactive({
    text=paste0(
    paste0("<b>", Hits()[(ok)]$efoId, "</b>"), 
    paste0("<br><b>", Hits()[(ok)]$trait, "</b>"),
    paste0("; N_gene = ", Hits()[(ok)]$traitNgene),
    ";<br>muScore = ", Hits()[(ok)]$muScore,
    "; muRank = ", Hits()[(ok)]$muRank,
    ";<br>N_study = ", Hits()[(ok)]$n_study,
    "; study_N = ", Hits()[(ok)]$study_N_mean, 
    "; N_snp = ", Hits()[(ok)]$n_snp,
    ";<br>N_snpw = ", Hits()[(ok)]$n_snpw,
    "; OR = ", round(Hits()[(ok)]$or_median, digits=2), 
    "; N_beta = ", Hits()[(ok)]$n_beta,
    "; pVal = ", sprintf("%.2g", 10^(-Hits()[(ok)]$pvalue_mlog_median)),
    "; RCRAS = ", round(Hits()[(ok)]$rcras, digits=2))
    return(text)
  })

  output$tigaPlot <- renderPlotly({
    xaxis <- list(title="Evidence (muScore)", type="normal", zeroline=F, showline=F)
    yaxis <- list(title="Effect (OddsRatio)", type=ifelse("Effect" %in% input$logaxes, "log", "normal"))
    axis_none <- list(zeroline=F, showline=F, showgrid=F, showticklabels=F)
    #message(sprintf("DEBUG: length(traitQryId()=%d; length(geneQryId()=%d", length(traitQryId()), length(geneQryId())))
    if (traitQryId()=="" & geneQryId()=="") {
      title <- "<I>(No query.)</I>"
      return(plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>%
               layout(title=title, xaxis=xaxis, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (traitQryId()=="" & geneQryId() %in% filtered$id) {
      title <- sprintf("<I>(No hits.)</I><br>Gene %s:<br>%s<br><I>Filtered by TIGA preprocessing<br>Reason: %s</I>", geneQryId(), geneQryName(), filtered[id==geneQryId(), reason])
      return(plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>%
               layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (is.null(Hits())) {
      title <- "<I>(No hits.)</I>"
      return(plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>%
               layout(title=title, xaxis=xaxis, yaxis=axis_none, margin=list(t=120,b=20)))
    }
    
    if (qryType()=="trait") {
      message(sprintf("DEBUG: %s", paste(collapse=",", paste(Hits()[(ok)]$geneSymbol, as.character(markerSize()), sep=":"))))
    } else if (qryType()=="gene") {
      message(sprintf("DEBUG: %s", paste(collapse=",", paste(Hits()[(ok)]$efoId, as.character(markerSize()), sep=":"))))
    }
    
    if (hitType()=="gene") {
      p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok)],
        x=~muScore,
        y=~or_median,
        color=~TDL, colors=c("gray", "black", "red", "green", "blue"),
        marker=list(symbol="circle", size=markerSize()),
        text=markerTextGenes()
      ) %>% config(displayModeBar=F) %>%
      layout(xaxis=xaxis, yaxis=yaxis, 
        title=paste0(traitQryName(), "<br>", "(", qryType(), ":", traitQryId(), ")"),
        margin=list(t=80,r=50,b=60,l=60), showlegend=T,
	legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
        font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits()[(ok)]), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else if (hitType()=="trait") {
      p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok)],
        x=~muScore,
        y=~or_median,
        marker=list(symbol="circle", size=markerSize()),
        text=markerTextTraits()
      ) %>% config(displayModeBar=F) %>%
      layout(xaxis=xaxis, yaxis=yaxis, 
        title=paste0(geneQryName(), "<br>", "(", qryType(), ":", geneQryId(), ")"),
        margin=list(t=80,r=50,b=60,l=60), showlegend=F,
	legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
        font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits()[(ok)]), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else { # "genetrait"
      p <- plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>%
               layout(title=title, xaxis=xaxis, yaxis=axis_none, margin=list(t=120,b=20))
    }
    return(p)
  })
  
  #DT numbers cols from 0.
  #"ensemblId","geneSymbol","geneName","geneFamily","TDL","n_study","study_N_mean","n_snp","n_snpw","geneNtrait","pvalue_mlog_median","or_median","rcras","muScore","muRank"
  output$hitrows <- DT::renderDataTable({
    if (is.null(Hits())) {
      return(NULL)
    } else if (hitType()=="gene") {
      return(DT::datatable(data=HitsWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap",
          	selection=list(target="row", mode="multiple", selected=NULL),
		colnames=c("ENSG", "GSYMB", "GeneName", "idgFam", "idgTDL", "N_study", "study_N", "N_snp", "N_snpw", "N_trait", "pVal_mlog", "OR", "N_beta", "RCRAS", "muScore", "muRank", "ok"),
          	options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
          	               list(className='dt-center', targets=c(0, 3:(ncol(HitsWithHtm())-2))),
          	               list(visible=F, targets=c(0, ncol(HitsWithHtm())-1))
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_median", "or_median", "rcras", "n_snpw"), digits=2))
  } else if (hitType()=="trait") {
      return(DT::datatable(data=HitsWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap",
		selection=list(target="row", mode="multiple", selected=NULL),
		colnames=c("efoId", "trait", "N_study", "study_N", "N_snp", "N_snpw", "N_gene", "pVal_mlog", "OR", "N_beta", "RCRAS", "muScore", "muRank", "ok"),
		options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
				list(className='dt-center', targets=c(0, 2:(ncol(HitsWithHtm())-2))),
				list(visible=F, targets=c(ncol(HitsWithHtm())-1))
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_median", "or_median", "rcras", "n_snpw"), digits=2)) 
  } else if (hitType()=="genetrait") {
      return(DT::datatable(data=HitsWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap",
		selection=list(target="row", mode="multiple", selected=NULL),
      		# efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, traitNgene, pvalue_mlog_median, or_median, betaN, rcras
		colnames=c("efoId", "trait", "ENSG", "GSYMB", "GeneName", "idgFam", "idgTDL", "N_study", "study_N", "N_snp", "N_snpw", "N_gene", "pVal_mlog", "OR", "N_beta", "RCRAS", "ok"),
		options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
				list(className='dt-center', targets=c(0, 2:(ncol(HitsWithHtm())-2))),
				list(visible=F, targets=c(2, ncol(HitsWithHtm())-1))
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_median", "or_median", "rcras", "n_snpw"), digits=2)) 
  } else {
    message(sprintf("ERROR: hitType()=%s; must be 'gene' or 'trait' or 'genetrait'.", hitType()))
    return(NULL)
  }}, server=T)

  #All-traits table has tiga-trait links.
  trait_tableHtm <- reactive({
    dt <- data.table(trait_table) #copy
    #dt[, idHtm := sprintf("<a href=\"%s\" target=\"_blank\">%s</a>", trait_uri, efoId)][, .(ID = idHtm, trait, N_study, N_gene)][order(trait)]
    dt[, idHtm := sprintf("<a href=\"%s?trait=%s\">%s</a>", urlBase(), efoId, efoId)]
    dt[, .(efoId = idHtm, trait, N_study, N_gene)][order(trait)]
  })

  output$traits <- DT::renderDataTable({
    DT::datatable(data=trait_tableHtm()[, .(efoId, trait, N_study, N_gene)], rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)

  #All-genes table has tiga-gene links.
  gene_tableHtm <- reactive({
    dt <- data.table(gene_table)[order(geneSymbol)]
    #dt[, symbHtm := sprintf("<a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\">%s</a>", geneSymbol, geneSymbol)][, .(Symbol = symbHtm, ensemblId, Name = geneName, Family = geneFamily, TDL, N_study, N_trait)][order(Symbol)]
    dt[, symbHtm := sprintf("<a href=\"%s?gene=%s\">%s</a>", urlBase(), ensemblId, geneSymbol)]
    dt[, .(ensemblId, geneSymbol = symbHtm, geneName, geneFamily, TDL, N_study, N_trait, filtered)]
  })

  output$genes <- DT::renderDataTable({
    DT::datatable(data=gene_tableHtm()[, .(geneSymbol, geneName, geneFamily, TDL, N_study, N_trait, filtered)], rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)
  
  study_tableHtm <- reactive({
    dt <- data.table(study_table)
    dt[, gcHtm := sprintf("<a href=\"https://www.ebi.ac.uk/gwas/studies/%s\">%s</a>", STUDY_ACCESSION, STUDY_ACCESSION)]
    dt[, pubmedHtm := sprintf("<a href=\"https://pubmed.ncbi.nlm.nih.gov/%s\">%s</a>", PUBMEDID, PUBMEDID)]
    dt[, .(Accession=gcHtm, Study=STUDY, PMID=pubmedHtm, DatePublished=DATE_PUBLISHED, DateAdded=DATE_ADDED_TO_CATALOG)][order(-DatePublished)]
  })

  output$studies <- DT::renderDataTable({
    DT::datatable(data=study_tableHtm()[, .(Accession, Study, PMID, DatePublished, DateAdded)], rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)

  output$association_detail_studies <- DT::renderDataTable({
    DT::datatable(data=data.table(STUDY_ACCESSION=rep("BLAH", 10), STUDY=rep("GROK", 10)), rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)

  output$association_detail_publications <- DT::renderDataTable({
    DT::datatable(data=data.table(PMID=rep("BLAH", 10), TITLE=rep("GROK", 10)), rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)

  Hits_export <- reactive({
    if (is.null(Hits())) { return(NULL) }
    hits_out <- data.table(Hits()) #copy
    if (hitType()=="gene") {
      hits_out[["efoId"]] <- traitQryId()
      hits_out[["trait"]] <- traitQryName()
      hits_out <- hits_out[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, geneNtrait, pvalue_mlog_median, or_median, n_beta, rcras, muScore, muRank)]
    } else if (hitType()=="trait") {
      hits_out[["ensemblId"]] <- geneQryId()
      hits_out[["geneName"]] <- geneQryName()
      hits_out <- hits_out[, .(ensemblId, geneName, efoId, trait, n_study, study_N_mean, n_snp, n_snpw, traitNgene, pvalue_mlog_median, or_median, n_beta, rcras, muScore, muRank)]
    }
    return(hits_out)
  })

  output$hits_file <- downloadHandler(
    filename = function() {
      if (hitType()=="gene") {
        sprintf("tiga_hits_%s.tsv", traitQryId())
      } else if (hitType()=="trait") {
        sprintf("tiga_hits_%s.tsv", geneQryId())
      } else {
        ("tiga_hits.tsv")
      }
    },
    content = function(file) {
      if (is.null(Hits_export())) { return(NULL) }
      write_delim(Hits_export(), file, delim="\t")
    }
  )
  output$traits_file <- downloadHandler(
    filename = "tiga_traits.tsv",
    content = function(file) {
      write_delim(trait_table, file, delim="\t")
    }
  )
  output$genes_file <- downloadHandler(
    filename = "tiga_genes.tsv",
    content = function(file) {
      write_delim(gene_table, file, delim="\t")
    }
  )
  output$gt_file <- downloadHandler(
    filename = "tiga_gene-trait_stats.tsv",
    content = function(file) {
      write_delim(gt, file, delim="\t")
    }
  )
}
###
shinyApp(ui, server)

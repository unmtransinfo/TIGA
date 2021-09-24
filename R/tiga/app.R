########################################################################################
### TIGA: Target Illumination GWAS Analytics
### gt = gene-trait data
### Input files:
###  gt_stats.tsv.gz            (from tiga_gt_stats.R) 
###  gt_provenance.tsv.gz       (from tiga_gt_provenance.R)
###  filtered_studies.tsv       (from tiga_gt_prepfilter.R)
###  filtered_traits.tsv        (from tiga_gt_prepfilter.R)
###  filtered_genes.tsv         (from tiga_gt_prepfilter.R)
###  gwascat_gwas.tsv           (from gwascat_gwas.R)
###  efo_graph.graphml.gz       (from efo_graph.R)
########################################################################################
### Requires shinysky, devtools::install_github("AnalytixWare/ShinySky")
########################################################################################
library(readr)
library(data.table)
library(igraph, quietly=T)
library(shiny, quietly=T)
library(shinyBS, quietly=T)
library(shinysky, quietly=T)
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
APPNAME_FULL <- "TIGA: Target Illumination GWAS Analytics"
MIN_ASSN <- 1
TDL_COLORS <- list(Tclin = "#B3D8FF", Tchem = "#B3FFCA", Tbio = "#FA838A", Tdark = "#A7A4A4")
#
if (!file.exists("tiga.Rdata")) {
  message(sprintf("ERROR: tiga.Rdata NOT FOUND."))
} else {
  message(sprintf("Loading tiga.Rdata..."))
  load("tiga.Rdata")
}
#
message(sprintf("Gene count, IDs: %d; symbols: %d", uniqueN(gt$ensemblId), uniqueN(gt$geneSymbol)))
message(sprintf("Trait count (total): %d", uniqueN(gt$efoId)))
message(sprintf("Trait count (filtered; n_assn<%d): %d", MIN_ASSN, uniqueN(gt$efoId)-length(trait_menu)))
message(sprintf("Trait count (menu; n_assn>=%d): %d", MIN_ASSN, length(trait_menu)))
#
message(sprintf("Provenance: %d PUBMEDIDS and %d studies for %d gene-trait pairs.", gt_prov[, uniqueN(PUBMEDID)], gt_prov[, uniqueN(STUDY_ACCESSION)], nrow(unique(gt_prov[, .(ensemblId, TRAIT_URI)]))))
flt <- data.table(table(filtered$type, filtered$reason))[N>0][order(V1, V2)]
message(sprintf("FILTERED entities: (%5s) %5d REASON: %s\n", flt$V1, flt$N, flt$V2))
#
message(sprintf("Graph \"%s\": vertices: %d; edges: %d", graph_attr(efoGraph, "name"), vcount(efoGraph), ecount(efoGraph)))
#
dbHtm <- sprintf("<B><i>Dataset</i></B>: GWAS Catalog version: %s; genes: %d; traits: %d; studies: %d; publications: %d", 
	GWASCATALOG_RELEASE, uniqueN(gt$ensemblId), uniqueN(gt$efoId), uniqueN(gt_prov$STUDY_ACCESSION), uniqueN(gt_prov$PUBMEDID))
#
###
idgfams <- c("GPCR", "Kinase", "IC", "NR", "Other")
TDLs <- c("Tclin", "Tchem", "Tbio", "Tdark")
#
#############################################################################
HelpHtm <- function() {
  htm <- sprintf("<P><B>TIGA</B>, Target Illumination GWAS Analytics, facilitates drug target illumination by 
scoring and ranking protein-coding genes associated with traits from genome-wide association studies
(GWAS). Similarly, <B>TIGA</B> can score and rank traits with the same gene-trait association metrics. 
Rather than a comprehensive analysis of GWAS for all biological implications and insights, this
focused application provides a rational method by which GWAS findings can be 
aggregated and filtered for applicable, actionable intelligence, with 
evidence usable by drug discovery scientists to enrich prioritization of target hypotheses. 
<P>
<B>Data sources:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">GWAS Catalog</a> [%s]
<LI><a href=\"https://www.ebi.ac.uk/efo/\" TARGET=\"_blank\">Experimental Factor Ontology (EFO)</a> [%s]
<LI><a href=\"http://juniper.health.unm.edu/tcrd/\" TARGET=\"_blank\">Target Central Resource Db (TCRD)</a> [%s]
</UL>", GWASCATALOG_RELEASE, EFO_RELEASE, TCRD_RELEASE)
  htm <- paste(sep="\n", htm, 
"<B>Notes:</B>
<UL>
<LI>Traits are mapped to EFO, Experimental Factor Ontology.
<LI>Mapped genes via Ensembl pipeline as per 
<a href=\"https://www.ebi.ac.uk/gwas/docs/methods/curation\" TARGET=\"_blank\">GWAS Catalog documentation</a>. 
Reported genes ignored for consistency and accountable confidence assessment in this app and downstream.
<LI>In this version, effect size measure (1) odds ratio (OR) or (2) BETA required.
<LI>Due to lack of reported, extracted, parsed and harmonized beta units, N_beta count is employed
as simple, rational measure of effect evidence and confidence (but not magnitude).
</UL>
<B>Aggregate statistics and scores:</B>
<UL>
  <LI><B>pVal_mLog<SUP>*</SUP></B>: max(-Log(pValue)) supporting gene-trait association.
  <LI><B>RCRAS<SUP>*</SUP></B>: Relative Citation Ratio (iCite RCR) Aggregated Score.
  <LI><B>N_snpw<SUP>*</SUP></B>: N_snp weighted by distance inverse exponential.
  <LI><B>N_study</B>: studies supporting gene-trait association (unique count).
  <LI><B>OR</B>: median(odds ratio, inverted if &lt;1) supporting gene-trait association.
  <LI><B>N_beta</B>: simple count of beta values with 95%% confidence intervals supporting gene-trait association.
  <LI><B>N_snp</B>: SNPs involved with gene-trait association (unique count).
  <LI><B>study_N</B>: mean(SAMPLE_SIZE) supporting gene-trait association.
  <LI><B>geneNtrait</B>: total traits associated with gene (unique count).
  <LI><B>traitNgene</B>: total genes associated with trait (unique count).
  <LI><B>traitNstudy</B>: total studies associated with trait (unique count).
  <LI><B>geneNstudy</B>: total studies associated with gene (unique count).
  <LI><B>meanRankScore</B>: Gene-trait pairs (GTs) are ranked based on selected variables, determined by benchmarking versus gold standard associations.  meanRankScore = 100 - Percentile(meanRank).
</UL>
<SUP>*</SUP>Variable used in <B>meanRankScore</B>.
<BR/>
Hits are ranked based on meanRankScore
<BR/>
<B>Other datatypes:</B>
<UL>
  <LI><B>ENSG</B>: Ensembl Gene ID.
  <LI><B>efoId</B>: EFO trait ID.
  <LI><B>STUDY_ACCESSION</B>: GWAS Catalog study ID.
  <LI><B>PMID</B>: PubMed ID.
  <LI><B>TDL</B>: IDG Target Development Level, a knowledge based classification. Tclin = high-confidence drug targets; Tchem = small-molecule modulator exists; Tbio = biological
function elucidated; Tdark = minimal knowledge.
</UL>
<BR/>
<B>UI:</B>
Scatterplot axes are Effect (OR or beta) vs. Evidence as measured by <B>meanRankScore</B>.
Odds ratio (OR) is the median, beta is a count of non-zero beta values, hence a
measure of effect-evidence but not magnitude. Missing ORs represented as 0 for plotting.
Note that for a given trait or gene, studies and associations may have ORs, betas, or collectively, both,
though the plot displays only the datatype selected or auto-chosen.
By default the axis is auto-chosen to reflect the predominant type.
<BR/>
Note that this app will accept query parameters <B>trait</B> (EFO_ID) and/or <B>gene</B>
(ENSEMBL_ID) via URL, e.g.
<B><TT>?trait=EFO_1000654</TT></B>, <B><TT>?gene=ENSG00000094914</TT></B>,
<B><TT>?trait=EFO_1000654&gene=ENSG00000094914</TT></B>.
<BR/>
<B>Authors:</B>
Jeremy Yang<SUP>1</SUP>, Dhouha Grissa<SUP>2</SUP>, Stephen Mathias<SUP>1</SUP>,
Cristian Bologa<SUP>1</SUP>, Anna Waller<SUP>1</SUP>, David Wild<SUP>3</SUP>,
Christophe Lambert<SUP>1</SUP>, Lars Juhl Jensen<SUP>2</SUP> and Tudor Oprea<SUP>1</SUP>.<BR/>
<I><SUP>1</SUP>University of New Mexico, Translational Informatics Division, Dept. of
Internal Medicine; <SUP>2</SUP>Novo Nordisk Foundation Center for Protein Research, Copenhagen,
Denmark; <SUP>3</SUP>Indiana University, School of Informatics, Computing and Engineering, Integrative Data Science Lab.</I>
<BR/>
<B>Publication:</B> <a href=\"https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab427/6292081\"> Yang et al, Bioinformatics, 04 June 2021</a>
<BR/>
<B>Feedback welcome</B> to <a href=\"mailto:jjyang_AT_salud_DOT_unm_DOT_edu\">Jeremy Yang</a>.<br/>
This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>")
  htm <- paste(htm, sprintf("<hr>\nPowered by: <tt>%s; %s</tt>", R.version.string, pkgVerTxt), sep="\n")
  return(htm)
}

##########################################################################################
ui <- fluidPage(
  tags$style(".green_class {color:#00ff00} .blue_class {color:#0000ff} .red_class {color:#ff0000} .black_class {color:black}"),
  #tags$style("label {display: table-cell; text-align: center; vertical-align: middle;} .form-group {display: table-row;}"),
  titlePanel(h2("IDG", tags$img(height="50", valign="bottom", src="IDG_logo_only.png"), APPNAME_FULL), windowTitle=APPNAME_FULL),
  fluidRow(
    column(4, 
      wellPanel(
        #tokens: A list whose length equal to nrow(local) where each element is array of string tokens.
        shinysky::textInput.typeahead(id="traitQry", placeholder="Trait...",
		local=trait_menu, valueKey="efoId", tokens=trait_menu$trait,
		template=HTML("<p class='repo-name'>{{efoId}}</p> <p class='repo-description'>{{trait}}</p>"), 
		limit=10), p(),
        shinysky::textInput.typeahead(id="geneQry", placeholder="Gene...",
		local=gene_menu, valueKey="ensemblId", tokens=gene_menu$geneSymbol,
		template=HTML("<p class='repo-name'>{{geneSymbol}}</p> <p class='repo-description'>{{geneName}}</p>"),
		limit=10), p(),
        	actionButton("goSubmit", label="Submit", icon=icon("cogs"), style='background-color:#EEEEEE;border-width:2px'),
        	actionButton("goReset", label="Reset", icon=icon("power-off"), style='background-color:#EEEEEE;border-width:2px')),
      wellPanel(
        fluidRow(column(2, HTML("<b><i>Hits</i></b>")),
        column(10, checkboxGroupInput("tdls", "TDLs", choiceValues=TDLs, choiceNames=TDLs, selected=c("Tclin", "Tchem", "Tbio", "Tdark"), inline=T)))),
      wellPanel(
        fluidRow(column(2, HTML("<b><i>Plot</i></b>")),
        column(10, radioButtons("yAxis", "Y-Axis", choiceNames=c("OR", "N_beta", "Auto"), choiceValues=c("or_median", "n_beta", "auto"), selected="auto", inline=T)))),
	wellPanel(htmlOutput(outputId="logHtm")),
	wellPanel(htmlOutput(outputId="resultHtm"))),
    column(8,
	tabsetPanel(id="tabset", type="tabs",
		tabPanel(value="hitsTable", title=htmlOutput("hitsTableTabHtm"), htmlOutput("tigaTableTitleHtm"), DT::dataTableOutput("tigaTable"), br(), downloadButton("hits_file",
			label=htmlOutput("hitsDownloadButtonHtm"))),
		tabPanel(value="hitsPlot", title=htmlOutput("hitsPlotTabHtm"), plotlyOutput("tigaPlot", height = "500px")),
		tabPanel(value="provenance", title=HTML("<center>Provenance<br/>(Gene-Trait)</center>"), htmlOutput("provenance_summary"),
			DT::dataTableOutput("provenance_studies"), br(),
			downloadButton("prov_detail_file", label="Download Provenance (this association)")),
		tabPanel(value="traits", title=HTML("<center><i>Traits<br/>(ALL)</i></center>"), DT::dataTableOutput("traits")),
		tabPanel(value="genes", title=HTML("<center><i>Genes<br/>(ALL)</i></center>"), DT::dataTableOutput("genes")),
		tabPanel(value="studies", title=HTML("<center><i>Studies<br/>(ALL)</i></center>"), DT::dataTableOutput("studies")),
		tabPanel(value="downloads", title=HTML("<i><br/>Downloads</i>"),
			h1("Downloads"),
			wellPanel(p(tags$b("ALL gene-trait associations:")), downloadButton("gt_file", label="Gene-Trait Associations"), htmlOutput("gtFileInfoHtm")),
			wellPanel(p(tags$b("ALL traits involved in gene-trait associations:")), downloadButton("traits_file", label="Traits"), htmlOutput("traitFileInfoHtm")),
			wellPanel(p(tags$b("ALL genes involved in gene-trait associations:")), downloadButton("genes_file", label="Genes"), htmlOutput("geneFileInfoHtm")),
			wellPanel(p(tags$b("Provenance for ALL gene-trait associations:")), downloadButton("provenance_file", label="Provenance"), htmlOutput("provFileInfoHtm")),
			HTML("For latest and archived files also see: <b><a href=\"https://unmtid-shinyapps.net/download/TIGA/\" target=\"_blank\">TIGA downloads directory</a></b>")
		),
		tabPanel(value="help", title=HTML("<i><br/>Help</i>"), htmlOutput("helpHtm"))
	))),
  hr(),
  fluidRow(
    column(12, tags$em(strong(sprintf("%s", APPNAME)), " web app from ", 
        tags$a(href="https://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
        " and ",
        tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
        " built from ",
        tags$a(href="https://www.ebi.ac.uk/gwas/", target="_blank", span(paste0("GWAS Catalog [", GWASCATALOG_RELEASE, "]"), tags$img(id="gwas_catalog_logo", height="50", valign="bottom", src="GWAS_Catalog_logo.png"))),
        " and ",
        tags$a(href="https://www.ebi.ac.uk/efo/", target="_blank", span(paste0("EFO [", EFO_RELEASE, "]"), tags$img(id="efo_logo", height="50", valign="bottom", src="EFO_logo.png")))
        ))),
  bsTooltip("traitQry", "Enter trait name fragment for autosuggest and resolution to EFO ID.", "bottom"),
  bsTooltip("geneQry", "Enter gene symbol or name fragment for autosuggest and resolution to Ensembl ID.", "bottom"),
  bsTooltip("goReset", "Reset.", "bottom"),
  bsTooltip("yAxis", "Auto chooses predominant effect size datatype. Datatype not selected viewable via hits table.", "bottom"),
  bsTooltip("tdls", "Filters hits shown but not inclusive hits download.", "bottom"),
  bsTooltip("tabset", "HitsPlot, HitsTable and Provenance depends on query and results; Traits, Genes, Downloads and Help are static.", "left"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "top"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "top"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "top"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "top")
)

##########################################################################################
server <- function(input, output, session) {

  output$helpHtm <- reactive({ paste(sprintf("<H2>%s Help</H2>", APPNAME), HelpHtm()) })

  Sys.sleep(1) #Needed?
  i_query <- 0 # initialize once per session
  
  # ?trait=EFO_0000341
  # ?trait=EFO_1000654&gene=ENSG00000094914
  httpQstr <- reactive({
    qStr <- getQueryString(session) #named list
    if (length(qStr)>0)
      for (key in names(qStr))
        message(sprintf("DEBUG: qStr[[\"%s\"]]=\"%s\"", key, qStr[[key]]))
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
 
  observeEvent(input$goReset, {
    #dqshiny::update_autocomplete_input(session, "traitQry", value="")
    shiny::updateTextInput(session, "traitQry", value="")

    #dqshiny::update_autocomplete_input(session, "geneQry", value="")
    shiny::updateTextInput(session, "geneQry", value="")

    updateQueryString("?", "push", session)
    session$reload()
  })
  
  efoId2Name <- function(efoId_this) {
    name <- trait_table[efoId==efoId_this, first(trait)]
    return(name)
  }
  ensemblId2Symbol <- function(ensemblId_this) {
    gsymb <- gene_table[ensemblId==ensemblId_this, first(geneSymbol)]
    return(gsymb)
  }
  ensemblId2Name <- function(ensemblId_this) {
    name <- gene_table[ensemblId==ensemblId_this, first(geneName)]
    return(name)
  }
  
  DetailSummaryHtm <- function(efoId_this, ensemblId_this) {
    #message(sprintf("DEBUG: efoId_this: %s; ensemblId_this: %s; efoId2Name(efoId_this): %s; ensemblId2Symbol(ensemblId_this): %s; ensemblId2Name(ensemblId_this): %s", efoId_this, ensemblId_this, efoId2Name(efoId_this), ensemblId2Symbol(ensemblId_this), ensemblId2Name(ensemblId_this)))
    htm <- sprintf("<table width=\"100%%\"><tr><td width=\"45%%\" align=\"right\" valign=\"bottom\"><h3>TRAIT: %s</br><i>%s</i></h3></td><td width=\"5%%\" align=\"center\"><h3>&#8226;</h3></td><td align=\"left\" valign=\"bottom\"><h3>GENE: %s<br/><i>%s (%s)</i></h3></td></tr></table>", efoId_this, efoId2Name(efoId_this), ensemblId_this, ensemblId2Symbol(ensemblId_this), ensemblId2Name(ensemblId_this))
    if (is.null(Hits()) | nrow(Hits())==0) {
      htm <- paste(htm, "<center><h3>NO ASSOCIATIONS FOUND</h3></center>", "\n")
    } else {
      htm <- paste(htm, "<center><h3>Gene-trait association provenance</h3></center>", "\n")
      keys <- sort(setdiff(names(Hits()), c("geneSymbol", "geneName", "trait", "ok2show")))
      for (k in keys)
        htm <- paste(htm, sprintf("<b>%s: </b><tt>%s</tt>", k, Hits()[[k]][1]), sep=" &#8226; ")
    }
    return(htm)
  }

  processUrlParams <- function(qStr) {
    message(sprintf("DEBUG: url: \"%s\"", urlText()))
    if ("trait" %in% names(qStr)) {
      #dqshiny::update_autocomplete_input(session, "traitQry", value=efoId2Name(qStr[["trait"]]))
      shiny::updateTextInput(session, "traitQry", value=efoId2Name(qStr[["trait"]]))
    }
    if ("gene" %in% names(qStr)) {
      #dqshiny::update_autocomplete_input(session, "geneQry", value=ensemblId2Symbol(qStr[["gene"]]))
      shiny::updateTextInput(session, "geneQry", value=ensemblId2Symbol(qStr[["gene"]]))
    }
  }

  # Returns both input fields as a list(gene = ***, trait = ***)
  qryIds <- reactive({
    input$goSubmit #trigger with goSubmit button?
    ids <- list(trait = NA, gene = NA)
    qStr <- httpQstr()
    if (i_query==0 & length(names(qStr))>0) {  #1st query may be via URL http param.
      processUrlParams(qStr)
      if ("trait" %in% names(qStr)) ids$trait <- qStr[["trait"]]
      if ("gene" %in% names(qStr)) ids$gene <- qStr[["gene"]]
    }
    if (!grepl("^\\s*$", input$traitQry)) ids$trait <- input$traitQry
    if (!grepl("^\\s*$", input$geneQry)) ids$gene <- input$geneQry
    message(sprintf("DEBUG: qryIds(): i_query: %d; input$traitQry: %s; ids$trait: %s; input$geneQry: %s; ids$gene: %s", i_query, input$traitQry, ids$trait, input$geneQry, ids$gene))
    i_query <<- i_query + 1  # Must assign to up-scoped variable.
    return(ids)
  })
  qryType <- reactive({
    if (!is.na(qryIds()$trait) & is.na(qryIds()$gene)) { return("trait") }
    else if (is.na(qryIds()$trait) & !is.na(qryIds()$gene)) { return("gene") }
    else if (!is.na(qryIds()$trait) & !is.na(qryIds()$gene)) { return("genetrait") }
    else { return(NA) } #No query
  })
  hitType <- reactive({
    ifelse(is.na(qryType()), NA,
    ifelse(qryType()=="gene", "trait",
    ifelse(qryType()=="trait", "gene",
    ifelse(qryType()=="genetrait", "genetrait", NA))))
  })
  traitQryName <- reactive({
    if (is.na(qryIds()$trait)) { return(NULL) }
    if (!(qryIds()$trait %in% trait_table$efoId)) { return(NULL) }
    return(trait_table[efoId==qryIds()$trait, trait])
  })
  geneQryName <- reactive({
    if (is.na(qryIds()$gene)) { return(NULL) }
    if (!(qryIds()$gene %in% gene_table$ensemblId)) { return(NULL) }
    return(sprintf("%s:%s", gene_table[ensemblId==qryIds()$gene, geneSymbol], gene_table[ensemblId==qryIds()$gene, geneName]))
  })

  Hits <- reactive({
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      return(NULL)
    }
    if (hitType()=="trait") {
      gt_this <- gt[ensemblId==qryIds()$gene]
      if (nrow(gt_this)==0) { return(NULL) }
      gt_this <- gt_this[, .(efoId, trait, pvalue_mlog_max, rcras, n_snpw, meanRankScore, n_study, study_N_mean, n_snp, traitNgene, or_median, n_beta)]
      setorder(gt_this, -meanRankScore)
    } else if (hitType()=="gene") {
      gt_this <- gt[efoId==qryIds()$trait]
      if (nrow(gt_this)==0) { return(NULL) }
      gt_this$TDL <- factor(gt_this$TDL, levels=c("NA", "Tdark", "Tbio", "Tchem", "Tclin"), ordered=T)
      gt_this <- gt_this[, .(ensemblId, geneSymbol, geneName, geneFamily, TDL, pvalue_mlog_max, rcras, n_snpw, meanRankScore, n_study, study_N_mean, n_snp, geneNtrait, or_median, n_beta)]
      setorder(gt_this, -meanRankScore)
    } else { #hitType=="genetrait"
      # These data appear in Provenance tab.
      gt_this <- gt[ensemblId==qryIds()$gene & efoId==qryIds()$trait]
      gt_this <- gt_this[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, pvalue_mlog_max, rcras, n_snpw, n_study, study_N_mean, n_snp, traitNgene, or_median, n_beta)]
    }
    gt_this[, ok2show := T]
    return(gt_this)
  })

  Hits2Show <- reactive({
    if (is.null(Hits())) {
      return(NULL)
    }
    gt_this <- Hits()[(ok2show)]
    if (hitType()=="gene") {
      gt_this <- gt_this[(TDL %in% input$tdls)]
    }
    return(gt_this)
  })
  Hits2Plot <- reactive({
    return(Hits2Show()[, `:=`(or_median = ifelse(is.na(or_median), 0, or_median), n_beta = ifelse(is.na(n_beta), 0, n_beta))]) #NA represented as zero for plotting.
  })

  yAxisAuto <- reactive({
    if (is.null(Hits())) {
      return(input$yAxis)
    } else if (input$yAxis=="auto") {
      n_or <- Hits()[!is.na(or_median), .N]
      n_nbeta <- Hits()[n_beta>0, .N]
      return(ifelse(n_nbeta>n_or, "n_beta", "or_median"))
    } else {
      return(input$yAxis)
    }
  })

  output$hitsPlotTabHtm <- reactive({ HTML(ifelse(!is.null(Hits()), sprintf("<center>HitsPlot<br/>(%d %ss)</center>", Hits2Show()[, .N], hitType()), "HitsPlot")) })
  output$hitsTableTabHtm <- reactive({ HTML(ifelse(!is.null(Hits()), sprintf("<center>HitsTable<br/>(%d %ss)</center>", Hits2Show()[, .N], hitType()), "<center>HitsTable<br/>(0)</center>")) })
  output$hitsDownloadButtonHtm <- reactive({ HTML(sprintf("Download Hits (%d %ss)", Hits()[, .N], hitType())) })

  output$gtFileInfoHtm <- reactive({ 
    paste(sep="\n", sprintf("<tt>rows: %d; cols: %d (%s)</tt>", nrow(gt), ncol(gt), paste(collapse=", ", names(gt))))
    })
  output$traitFileInfoHtm <- reactive({ 
    paste("\n", sprintf("<tt>rows: %d; cols: %d (%s)</tt>", nrow(trait_table), ncol(trait_table), paste(collapse=", ", names(trait_table))))
    })
  output$geneFileInfoHtm <- reactive({
    paste("\n", sprintf("<tt>rows: %d; cols: %d (%s)</tt>", nrow(gene_table[(!filtered)]), ncol(gene_table), paste(collapse=", ", names(gene_table))))
    })
  output$provFileInfoHtm <- reactive({
    paste("\n", sprintf("<tt>rows: %d; cols: %d (%s)</tt>", nrow(gt_prov), ncol(gt_prov), paste(collapse=", ", names(gt_prov))))
    })
  
  output$provenance_summary <- reactive({ DetailSummaryHtm(qryIds()$trait, qryIds()$gene) })

  output$tigaTableTitleHtm <- reactive({
    if (is.na(qryType())) { htm <- "" }
    else if (qryType() == "trait") { htm <- sprintf("<H3>TRAIT: \"%s\" (<a target=\"_blank\" href=\"%s\">%s</a>)</H3>", traitQryName(), efoId2Uri(qryIds()$trait), qryIds()$trait) }
    else if (qryType() == "gene") { htm <- sprintf("<H3>GENE: \"%s\" (<a target=\"_blank\" href=\"https://pharos.nih.gov/targets/%s\">%s</a>)</H3>", geneQryName(), qryIds()$gene, qryIds()$gene) } 
    else { htm <- "" }
    return(htm)
  })

  #Hits table has links to tiga:genetrait, and to external resources EFO and Pharos.
  Hits2ShowWithHtm <- reactive({
    hwh <- data.table(Hits2Show()) #copy
    if (hitType() == "trait") {
      hwh <- hwh[, efoId := sprintf("%s<br/><a href=\"%s?trait=%s&gene=%s\"><i class=\"fa fa-search\"></i></a> &nbsp; <a href=\"%s\" target=\"_blank\"><i class=\"fa fa-external-link\"></i></a>", efoId, urlBase(), efoId, qryIds()$gene, sapply(efoId, efoId2Uri))]
    } else if (hitType() == "gene") {
      hwh <- hwh[, geneSymbol := sprintf("%s<br/><a href=\"%s?trait=%s&gene=%s\"><i class=\"fa fa-search\"></i></a> &nbsp; <a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\"><i class=\"fa fa-external-link\"></i></a>", geneSymbol, urlBase(), qryIds()$trait, ensemblId, geneSymbol)]
    }
    return(hwh)
  })

  # Hide/show/select tabs, based on query inputs.
  observeEvent(qryIds(), {
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      hideTab("tabset", "hitsPlot")
      hideTab("tabset", "hitsTable")
      hideTab("tabset", "provenance")
      updateTabsetPanel(session, "tabset", selected="traits")
    } else if (qryType() %in% c("gene", "trait")) {
      showTab("tabset", "hitsPlot")
      showTab("tabset", "hitsTable")
      hideTab("tabset", "provenance")
      updateTabsetPanel(session, "tabset", selected="hitsTable")
    } else if (qryType()=="genetrait") {
      hideTab("tabset", "hitsPlot")
      hideTab("tabset", "hitsTable")
      showTab("tabset", "provenance")
      updateTabsetPanel(session, "tabset", selected="provenance")
    }
  })

  output$resultHtm <- reactive({
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      return("No query. Search? Browse?")
    }
    message(sprintf("TraitQuery: \"%s\" (%s)", traitQryName(), qryIds()$trait))
    message(sprintf("GeneQuery: \"%s\" (%s)", geneQryName(), qryIds()$gene))
    htm <- sprintf("<B>Results:</B>")
    if (qryType() %in% c("trait", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", traitQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", efoId2Uri(qryIds()$trait), qryIds()$trait))
    }
    if (qryType() %in% c("gene", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", geneQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"https://pharos.nih.gov/targets/%s\">%s</a>)", qryIds()$gene, qryIds()$gene))
    }
    if (!is.null(Hits())) {
      htm <- paste0(htm, sprintf("; N_%s: %d shown (%d total)", hitType(), Hits2Show()[, .N], Hits()[, .N]))
      #if (!is.null(Hits()[["or_median"]]) & !is.null(Hits()[["n_beta"]]))
      #  message(sprintf("ORs: %d; N_betas>0: %d", Hits()[or_median>0, .N], Hits()[n_beta>0, .N]))
    } else if (qryIds()$trait %in% filtered$id) {
      htm <- paste0(htm, sprintf("; %s <B>%s: %s</B> filtered; reason: %s.", hitType(), qryIds()$trait, traitQryName(), filtered[id==qryIds()$trait, reason]))
    } else if (qryIds()$gene %in% filtered$id) {
      htm <- paste0(htm, sprintf("; %s <B>%s: %s</B> filtered; reason: %s.", hitType(), qryIds()$gene, geneQryName(), filtered[id==qryIds()$gene, reason]))
    } else {
      htm <- paste0(htm, sprintf("; No %ss found.", hitType()))
    }
    return(htm)
  })

  output$logHtm <- reactive({
    htm <- dbHtm
    return(htm)
  })
  
  markerSize <- reactive({ return(10) })

  markerTextGenes <- reactive({
    text <- paste0(
    paste0("<b>", Hits2Plot()[, geneSymbol], "</b> (", Hits2Plot()[, ensemblId], ")"), 
    paste0("<br><b>", Hits2Plot()[, geneName], "</b>"),
    paste0("<br>Fam:", Hits2Plot()[, geneFamily]),
    paste0(", TDL:", Hits2Plot()[, TDL]),
    paste0("; geneNtrait = ", Hits2Plot()[, geneNtrait]),
    ";<br>meanRankScore = ", round(Hits2Plot()[, meanRankScore], digits=3),
    ";<br>N_study = ", Hits2Plot()[, n_study],
    "; study_N = ", Hits2Plot()[, study_N_mean], 
    "; N_snp = ", Hits2Plot()[, n_snp],
    ";<br>N_snpw = ", Hits2Plot()[, n_snpw],
    "; OR = ", round(Hits2Plot()[, or_median], digits=2), 
    "; N_beta = ", Hits2Plot()[, n_beta],
    "; pVal = ", sprintf("%.2g", 10^(-Hits2Plot()[, pvalue_mlog_max])),
    "; RCRAS = ", round(Hits2Plot()[, rcras], digits=2)
    # ";<br>DEBUG: .I = ", Hits2Plot()[, .I]
      )
    return(text)
  })
  markerTextTraits <- reactive({
    text <- paste0(
    paste0("<b>", Hits2Plot()[, efoId], "</b>"), 
    paste0("<br><b>", Hits2Plot()[, trait], "</b>"),
    paste0("; traitNgene = ", Hits2Plot()[, traitNgene]),
    ";<br>meanRankScore = ", round(Hits2Plot()[, meanRankScore], digits=3),
    ";<br>N_study = ", Hits2Plot()[, n_study],
    "; study_N = ", Hits2Plot()[, study_N_mean], 
    "; N_snp = ", Hits2Plot()[, n_snp],
    ";<br>N_snpw = ", Hits2Plot()[, n_snpw],
    "; OR = ", round(Hits2Plot()[, or_median], digits=2), 
    "; N_beta = ", Hits2Plot()[, n_beta],
    "; pVal = ", sprintf("%.2g", 10^(-Hits2Plot()[, pvalue_mlog_max])),
    "; RCRAS = ", round(Hits2Plot()[, rcras], digits=2))
    return(text)
  })

  output$tigaPlot <- renderPlotly({
    xaxis <- list(title="Evidence (meanRankScore)", type="normal", zeroline=F, showline=F, dtick=10)
    #yaxis <- list(title=ifelse(yAxisAuto()=="n_beta", "N_Beta", "Effect (OddsRatio)"), type="normal", dtick=ifelse(yAxisAuto()=="n_beta", NA, 1))
    yaxis <- list(title=ifelse(yAxisAuto()=="n_beta", "N_Beta", "Effect (OddsRatio)"), type="normal", dtick=NA)
    axis_none <- list(zeroline=F, showline=F, showgrid=F, showticklabels=F)
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      title <- "<I>(No query.)</I>"
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (is.na(qryIds()$trait) & qryIds()$gene %in% filtered$id) {
      title <- sprintf("<I>(No hits.)</I><br>Gene %s:<br>%s<br><I>Filtered by TIGA preprocessing<br>Reason: %s</I>", qryIds()$gene, geneQryName(), filtered[id==qryIds()$gene, reason])
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (is.null(Hits2Plot())) {
      title <- "<I>(No hits.)</I>"
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    }

    # Fixed x-range, and y-range for a given (all) hitlist, to avoid moving points.
    xaxis[["range"]] <- c(-2, 102) #Extra space for markers.
    if (yAxisAuto()=="n_beta") {
      yaxis[["range"]] <- c(floor(min(Hits()[, n_beta])-.1), ceiling(max(Hits()[, n_beta])+.1))
    } else if (yAxisAuto()=="or_median") {
      yaxis[["range"]] <- c(-0.1, ceiling(max(Hits()[, or_median])+.1))
    }

    if (hitType()=="gene") {
      if (yAxisAuto()=="n_beta")
        p <- plot_ly(type='scatter', mode='markers', data=Hits2Plot(), 
          x = ~meanRankScore, 
          y = ~n_beta,
                marker=list(symbol="circle", size=markerSize()), text=markerTextGenes(),
                color=~TDL, colors=c("gray", "black", "red", "green", "blue"))
      else # or_median
        p <- plot_ly(type='scatter', mode='markers', data=Hits2Plot(), 
          x = ~meanRankScore, 
          y = ~or_median,
                marker=list(symbol="circle", size=markerSize()), text=markerTextGenes(),
                color=~TDL, colors=c("gray", "black", "red", "green", "blue"))
      p <- config(p, displayModeBar=F) %>% layout(xaxis=xaxis, yaxis=yaxis, 
          title=paste0(toupper(qryType()), ":\"", traitQryName(), "\"<br>", "(", qryIds()$trait, ")"),
          margin=list(t=80,r=50,b=60,l=60), showlegend=T,
	  legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
          font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits2Plot()), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else if (hitType()=="trait") {
      if (yAxisAuto()=="n_beta")
        p <- plot_ly(type='scatter', mode='markers', data=Hits2Plot(), x=~meanRankScore, y=~n_beta,
                marker=list(symbol="circle", size=markerSize()), text=markerTextTraits())
      else # or_median
        p <- plot_ly(type='scatter', mode='markers', data=Hits2Plot(), x=~meanRankScore, y=~or_median,
                marker=list(symbol="circle", size=markerSize()), text=markerTextTraits())
      p <- config(p, displayModeBar=F) %>% layout(xaxis=xaxis, yaxis=yaxis, 
          title=paste0(toupper(qryType()), ":\"", geneQryName(), "\"<br>", "(", qryIds()$gene, ")"),
          margin=list(t=80,r=50,b=60,l=60), showlegend=F,
	  legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
          font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits2Plot()), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else { # "genetrait"
      title <- "<I>(No plot in gene+trait mode.)</I>"
      p <- plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20))
    }
    return(p)
  })

  #DT numbers cols from 0.
  #"ensemblId","geneSymbol","geneName","geneFamily","TDL","pvalue_mlog_max","rcras","n_snpw","meanRankScore","n_study","study_N_mean","n_snp","geneNtrait","or_median","n_beta"
  output$tigaTable <- DT::renderDataTable({
    if (is.null(Hits())) return(NULL)
    if (hitType()=="gene") {
      return(DT::datatable(data=Hits2ShowWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap", selection="none",
		colnames=c("ENSG", "GSYMB", "GeneName", "idgFam", "idgTDL", "pVal_mlog", "RCRAS", "N_snpw", "meanRankScore", "N_study", "study_N", "N_snp", "geneNtrait", "OR", "N_beta", "ok2show"),
          	options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
          	               list(className='dt-center', targets=c(0, 1, 3:(ncol(Hits2ShowWithHtm())-2))),
          	               list(visible=F, targets=c(0, ncol(Hits2ShowWithHtm())-1)) #Hide EnsemblId
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_max", "or_median", "rcras", "n_snpw", "meanRankScore"), digits=2)
    %>% DT::formatStyle(c("pvalue_mlog_max", "rcras", "n_snpw"), fontWeight="bold")
    %>% DT::formatStyle(c("meanRankScore"), color="black", fontWeight="bold")
    %>% DT::formatStyle("TDL", backgroundColor=styleEqual(c("Tclin", "Tchem", "Tbio", "Tdark"), c(TDL_COLORS[["Tclin"]], TDL_COLORS[["Tchem"]], TDL_COLORS[["Tbio"]], TDL_COLORS[["Tdark"]]))
  )
        )
  } else if (hitType()=="trait") {
      return(DT::datatable(data=Hits2ShowWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap", selection="none",
		colnames=c("efoId", "trait", "pVal_mlog", "RCRAS", "N_snpw", "meanRankScore", "N_study", "study_N", "N_snp", "traitNgene", "OR", "N_beta", "ok2show"),
		options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
				list(className='dt-center', targets=c(0, 1, 3:(ncol(Hits2ShowWithHtm())-2))),
				list(visible=F, targets=c(ncol(Hits2ShowWithHtm())-1))
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_max", "or_median", "rcras", "n_snpw", "meanRankScore"), digits=2)
    %>% DT::formatStyle(c("pvalue_mlog_max", "rcras", "n_snpw"), fontWeight="bold")
    %>% DT::formatStyle(c("meanRankScore"), color="black", fontWeight="bold")
	) 
  }
  }, server=T)

  #All-traits table has tiga-trait links.
  trait_tableHtm <- reactive({
    dt <- data.table(trait_table) #copy
    dt[, efoHtm := sprintf("<a href=\"%s?trait=%s\">%s</a>", urlBase(), efoId, efoId)]
    dt[, .(efoId = efoHtm, trait, N_study, traitNgene)][order(trait)]
  })

  output$traits <- DT::renderDataTable({
    DT::datatable(data=trait_tableHtm()[, .(efoId, trait, N_study, traitNgene)], rownames=F, width="100%", options=list(autoWidth=F, dom='tipf'), escape=F, selection="none")
  }, server=T)

  #All-genes table has tiga-gene links.
  gene_tableHtm <- reactive({
    dt <- data.table(gene_table[(!filtered)])[order(geneSymbol)]
    dt[, symbHtm := sprintf("<a href=\"%s?gene=%s\">%s</a>", urlBase(), ensemblId, geneSymbol)]
    dt[, .(ensemblId, geneSymbol = symbHtm, geneName, geneFamily, TDL, N_study, geneNtrait)]
  })

  output$genes <- DT::renderDataTable({
    DT::datatable(data=gene_tableHtm()[, .(geneSymbol, geneName, geneFamily, TDL, N_study, geneNtrait)], rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F, selection="none")
  }, server=T)
  
  study_tableHtm <- reactive({
    dt <- data.table(study_table)
    dt[, gcHtm := sprintf("<a target=\"_blank\" href=\"https://www.ebi.ac.uk/gwas/studies/%s\">%s<i class=\"fa fa-external-link\"></a>", STUDY_ACCESSION, STUDY_ACCESSION)]
    dt[, pubmedHtm := sprintf("<a target=\"_blank\" href=\"https://pubmed.ncbi.nlm.nih.gov/%s\">%s<i class=\"fa fa-external-link\"></a>", PUBMEDID, PUBMEDID)]
    dt[, efoId := sub("^.*/", "", MAPPED_TRAIT_URI)]
    dt[, efoHtm := sprintf("<a href=\"%s?trait=%s\">%s</a>", urlBase(), efoId, efoId)]
    dt[, .(Trait = MAPPED_TRAIT, efoId = efoHtm, Study=STUDY, Accession=gcHtm, PMID=pubmedHtm, DatePublished=DATE_PUBLISHED)][order(Trait)]
  })

  output$studies <- DT::renderDataTable({
    DT::datatable(data=study_tableHtm()[, .(Trait, efoId, Study, Accession, PMID, DatePublished)], rownames=F, width="100%", options=list(autoWidth=F, dom='tipf'), escape=F, selection="none")
  }, server=T)

  output$provenance_studies <- DT::renderDataTable({
    DT::datatable(data=provenance_studies_tableHtm(), rownames=F, width="100%", options=list(autoWidth=F, dom='tip'), escape=F, selection="none",
      caption = tags$caption(style='caption-side:top; text-align:center;', 'Table: ', tags$b('Studies with association evidence')))
  }, server=T)

  provenance_studies_tableHtm <- reactive({
    dt <- data.table(merge(gt_prov[efoId == qryIds()$trait & ensemblId == qryIds()$gene, .(STUDY_ACCESSION)], study_table[, .(STUDY_ACCESSION, STUDY, PUBMEDID, DATE_PUBLISHED)], by="STUDY_ACCESSION", all.x=T, all.y=F))
    dt <- unique(dt)
    dt[, gcHtm := sprintf("<a target=\"_blank\" href=\"https://www.ebi.ac.uk/gwas/studies/%s\">%s<i class=\"fa fa-external-link\"></a>", STUDY_ACCESSION, STUDY_ACCESSION)]
    dt[, pubmedHtm := sprintf("<a target=\"_blank\" href=\"https://pubmed.ncbi.nlm.nih.gov/%s\">%s<i class=\"fa fa-external-link\"></a>", PUBMEDID, PUBMEDID)]
    dt[, .(Accession=gcHtm, Study=STUDY, PMID=pubmedHtm, DatePublished=DATE_PUBLISHED)][order(-DatePublished)]
  })

  Provenance_export <- reactive({
    prov_out <- data.table(merge(gt_prov[efoId == qryIds()$trait & ensemblId == qryIds()$gene, .(efoId, ensemblId, STUDY_ACCESSION)], study_table[, .(STUDY_ACCESSION, STUDY, PUBMEDID, DATE_PUBLISHED)], by="STUDY_ACCESSION", all.x=T, all.y=F))
    if (nrow(prov_out)==0) { return(NULL) }
    prov_out <- merge(prov_out, gene_table[, .(ensemblId, geneSymbol, geneName)], by="ensemblId")
    prov_out <- merge(prov_out, trait_table[, .(efoId, trait)], by="efoId")
    prov_out <- prov_out[, .(geneSymbol, ensemblId, geneName, efoId, trait, STUDY_ACCESSION, STUDY, DATE_PUBLISHED, PUBMEDID)]
    prov_out <- unique(prov_out)
    return(prov_out)
  })

  Hits_export <- reactive({
    if (is.null(Hits())) { return(NULL) }
    hits_out <- data.table(Hits()) #copy
    if (hitType()=="gene") {
      hits_out[["efoId"]] <- qryIds()$trait
      hits_out[["trait"]] <- traitQryName()
      hits_out <- hits_out[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, pvalue_mlog_max, rcras, meanRankScore, study_N_mean, n_snp, n_snpw, geneNtrait, or_median, n_beta)]
    } else if (hitType()=="trait") {
      hits_out[["ensemblId"]] <- qryIds()$gene
      hits_out[["geneName"]] <- geneQryName()
      hits_out <- hits_out[, .(ensemblId, geneName, efoId, trait, n_study, pvalue_mlog_max, rcras, meanRankScore, study_N_mean, n_snp, n_snpw, traitNgene, or_median, n_beta)]
    }
    return(hits_out)
  })

  output$hits_file <- downloadHandler(
    filename = function() {
      if (hitType()=="gene") {
        sprintf("tiga_hits_%s.tsv", qryIds()$trait)
      } else if (hitType()=="trait") {
        sprintf("tiga_hits_%s.tsv", qryIds()$gene)
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
      write_delim(gene_table[(!filtered)], file, delim="\t")
    }
  )
  output$gt_file <- downloadHandler(
    filename = "tiga_gene-trait_stats.tsv",
    content = function(file) {
      write_delim(gt, file, delim="\t")
    }
  )
  output$prov_detail_file <- downloadHandler(
    filename = function() {
      sprintf("tiga_provenance_%s_%s.tsv", qryIds()$gene, qryIds()$trait)
    },
    content = function(file) {
      if (is.null(Provenance_export())) { return(NULL) }
      write_delim(Provenance_export(), file, delim="\t")
    }
  )
  output$provenance_file <- downloadHandler(
    filename = "tiga_gene-trait_provenance.tsv",
    content = function(file) {
      write_delim(gt_prov, file, delim="\t")
    }
  )
}
###
shinyApp(ui, server)

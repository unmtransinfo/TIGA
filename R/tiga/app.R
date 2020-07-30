########################################################################################
### TIGA: Target Illumination GWAS Analytics
### gt = gene-trait data
### Input files:
###  gt_stats.tsv.gz               (from tiga_gt_stats.R) 
###  gt_provenance.tsv.gz          (from tiga_gt_stats.R)
###  filtered_studies.tsv.gz       (from tiga_gt_stats.R)
###  filtered_traits.tsv.gz        (from tiga_gt_stats.R)
###  filtered_genes.tsv.gz         (from tiga_gt_stats.R)
###  filtered_studies_trait.tsv.gz (from gwascat_trait.R)
###  gwascat_gwas.tsv.gz           (from gwascat_gwas.R)
###  efo_graph.graphml.gz          (from efo_graph.R)
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
APPNAME_FULL <- "TIGA: Target Illumination GWAS Analytics"
MIN_ASSN <- 1
#
t0 <- proc.time()
DEBUG <- F
if (!file.exists("tiga.Rdata") | DEBUG) {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  gt <- read_delim("data/gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(), geneNtrait=col_integer(), geneNstudy=col_integer(), traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(), geneMeanRank=col_double(), geneMeanRankScore=col_double(), traitMeanRank=col_double(), traitMeanRankScore=col_double()))
  setDT(gt)
  setnames(gt, old=c("geneIdgTdl"), new=c("TDL"))
  #
  gt_prov <- read_delim("data/gt_provenance.tsv.gz", "\t", col_types=cols(.default=col_character()))
  setDT(gt_prov)
  #
  filtered_studies <- read_delim("data/filtered_studies.tsv.gz", "\t")
  setDT(filtered_studies)
  filtered_studies_trait <- read_delim("data/filtered_studies_trait.tsv.gz", "\t")
  setDT(filtered_studies_trait)
  filtered_studies <- rbindlist(list(filtered_studies, filtered_studies_trait))
  filtered_studies[, type := "study"]
  filtered_traits <- read_delim("data/filtered_traits.tsv.gz", "\t")
  setDT(filtered_traits)
  filtered_traits[, type := "trait"]
  filtered_genes <- read_delim("data/filtered_genes.tsv.gz", "\t")
  setDT(filtered_genes)
  filtered_genes[, type := "gene"]
  filtered_gene_menu <- filtered_genes$ensemblId #named vector
  names(filtered_gene_menu) <- sprintf("%s:%s", filtered_genes$ensemblSymb, filtered_genes$geneName)
  filtered <- rbindlist(list(filtered_studies[, .(type, id=STUDY_ACCESSION, reason)], filtered_traits[, .(type, id=TRAIT_URI, reason)], filtered_genes[, .(type, id=ensemblId, reason)]))
  #
  trait_table <- gt[, .(trait=first(trait), N_study=first(traitNstudy), N_gene=uniqueN(ensemblId)),  by="efoId"]
  message(sprintf("Traits with N_assn<%d: %d", MIN_ASSN, trait_table[N_gene<MIN_ASSN, uniqueN(efoId)]))
  trait_table <- trait_table[N_gene>=MIN_ASSN]
  trait_table[, trait_uri := sapply(efoId, efoId2Uri)]
  trait_table <- trait_table[, .(trait_uri, efoId, trait, N_study, N_gene)][order(trait)]
  trait_menu <- trait_table$efoId #named vector
  names(trait_menu) <- trait_table$trait
  #
  gene_table <- gt[, .(geneSymbol, geneName, geneFamily, TDL, N_study = geneNstudy, N_trait = uniqueN(efoId), filtered=F),  by=c("ensemblId")]
  gene_table <- rbindlist(list(gene_table, filtered_genes[, .(ensemblId, geneSymbol=ensemblSymb, geneName, geneFamily, TDL, N_study=0, N_trait=0, filtered=T)]))
  gene_table <- unique(gene_table)[order(geneSymbol)]
  #
  gene_menu <- gene_table$ensemblId #named vector
  names(gene_menu) <- sprintf("%s:%s", gene_table$geneSymbol, gene_table$geneName)
  #
  trait_menu <- as.list(trait_menu)
  gene_menu <- as.list(gene_menu)
  #
  study_table <- read_delim("data/gwascat_gwas.tsv.gz", "\t", col_types = cols(.default = col_character(), DATE=col_date(), DATE_ADDED_TO_CATALOG=col_date()))
  setDT(study_table)
  study_table <- study_table[, .(STUDY_ACCESSION, STUDY, PUBMEDID, DATE_PUBLISHED = DATE, DATE_ADDED_TO_CATALOG)][order(DATE_PUBLISHED)]
  #
  system("if [ -f \"data/efo_graph.graphml.gz\" ]; then gunzip -f data/efo_graph.graphml.gz ; fi")
  efoGraph <- read_graph("data/efo_graph.graphml", format="graphml")
  system("if [ -f \"data/efo_graph.graphml\" ]; then gzip -f data/efo_graph.graphml ; fi")
  #
  save(gt, trait_table, gene_table, trait_menu, gene_menu, filtered_gene_menu, gt_prov, filtered, study_table, efoGraph, file="tiga.Rdata")
} else {
  message(sprintf("Loading tiga.Rdata..."))
  load("tiga.Rdata")
}
if (!("efoId" %in% names(gt_prov))) gt_prov[, efoId := sub("^.*/", "", TRAIT_URI)] #tmp
#
t_elapsed <- (proc.time()-t0)[3]
#
message(sprintf("Gene count, IDs: %d; symbols: %d", uniqueN(gt$ensemblId), uniqueN(gt$geneSymbol)))
message(sprintf("Trait count (total): %d", uniqueN(gt$efoId)))
#trait_table[, ontology := factor(sub("_.*$","", efoId))]
#trait_counts <- trait_table[, .N, by="ontology"][order(-N)]
#message(sprintf("traits (%10s): %4d / %4d (%4.1f%%)\n", trait_counts$ontology, trait_counts$N, sum(trait_counts$N), 100*trait_counts$N/sum(trait_counts$N)))
message(sprintf("Trait count (filtered; n_assn<%d): %d", MIN_ASSN, uniqueN(gt$efoId)-length(trait_menu)))
message(sprintf("Trait count (menu; n_assn>=%d): %d", MIN_ASSN, length(trait_menu)))
#
message(sprintf("Provenance: %d PUBMEDIDS and %d studies for %d gene-trait pairs.", gt_prov[, uniqueN(PUBMEDID)], gt_prov[, uniqueN(STUDY_ACCESSION)], nrow(unique(gt_prov[, .(ensemblId, TRAIT_URI)]))))
flt <- data.table(table(filtered$type, filtered$reason))[N>0][order(V1, V2)]
message(sprintf("FILTERED entities: (%5s) %5d REASON: %s\n", flt$V1, flt$N, flt$V2))
#
message(sprintf("Graph \"%s\": vertices: %d; edges: %d", graph_attr(efoGraph, "name"), vcount(efoGraph), ecount(efoGraph)))
#
dbHtm <- sprintf("<B>Dataset:</B> genes: %d; traits: %d ; studies: %d; publications: %d", 
                 uniqueN(gt$ensemblId), uniqueN(gt$efoId), uniqueN(gt_prov$STUDY_ACCESSION), uniqueN(gt_prov$PUBMEDID))
#
###
idgfams <- c("GPCR", "Kinase", "IC", "NR", "Other")
#axes <- c("Effect", "Evidence")
#
#############################################################################
HelpHtm <- function() {
  htm <- ("<P><B>TIGA</B>, Target Illumination GWAS Analytics, facilitates drug target illumination by 
scoring and ranking protein-coding genes associated with traits from genome-wide association studies
(GWAS). Similarly, <B>TIGA</B> can score and rank traits with the same gene-trait association metrics. 
Rather than a comprehensive analysis of GWAS for all biological implications and insights, this
focused application provides a rational method by which GWAS findings can be 
aggregated and filtered for applicable, actionable intelligence, with 
evidence usable by drug discovery scientists to enrich prioritization of target hypotheses. 
Data from the <A HREF=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">NHGRI-EBI GWAS Catalog</A>.
<UL>
<LI>Traits are mapped to EFO, Experimental Factor Ontology.
<LI>Mapped genes via Ensembl pipeline as per GWAS Catalog documentation. Reported genes ignored for consistency and accountable
confidence assessment in this app and downstream.
<LI>In this version, effect size measure (1) odds ratio (OR) or (2) BETA required.
<LI>Due to lack of reported, extracted, parsed and harmonized beta units, N_beta count is employed
as simple, rational measure of effect evidence and confidence (but not magnitude).
</UL>
<B>Datatypes:</B>
<UL>
  <LI><B>N_study<SUP>*</SUP></B>: studies supporting trait-gene association.
  <LI><B>pVal_mLog<SUP>*</SUP></B>: median(-Log(pValue)) supporting trait-gene association.
  <LI><B>RCRAS<SUP>*</SUP></B>: Relative Citation Ratio (RCR) Aggregated Score (iCite-RCR-based)
  <LI><B>OR</B>: median(odds ratio, inverted if &lt;1) supporting trait-gene association (computed as one if missing).
  <LI><B>N_beta</B>: simple count of beta values with 95% confidence intervals supporting trait-gene association.
  <LI><B>N_snp</B>: SNPs involved with trait-gene association.
  <LI><B>N_snpw</B>: N_snp weighted by distance inverse exponential.
  <LI><B>study_N</B>: mean(SAMPLE_SIZE) supporting trait-gene association.
  <LI><B>geneNtrait</B>: total traits associated with gene.
  <LI><B>traitNgene</B>: total genes associated with trait.
  <LI><B>meanRank</B>: For a given query trait, genes are ranked, or for a given gene, traits are ranked, based on selected variables, determined by benchmarking versus gold standard curated disease-gene associations.
  <LI><B>meanRankScore</B>: 1/meanRank, for normalization (0,1] and visualization.
</UL>
<SUP>*</SUP>Variable used in <B>meanRank</B>.
<BR/>
Hits are ranked based on meanRankScore
<BR/>
<B>UI:</B>
Scatterplot axes are Effect (OR or beta) vs. Evidence as measured by <B>meanRankScore</B>.
Odds ratio (OR) is the median, beta is a count of non-zero beta values, hence a
measure of effect-evidence but not magnitude. Nonexistent ORs plotted as zero.
<UL>
<LI>Genes markers colored by TDL, and may be be hidden via legend.
<LI>Trait markers colored by EFO top-level class, and may be hidden via legend (TO DO).
<LI>Plot markers may be sized by <B>N_study</B> or <B>RCRAS</B>.
<LI>Note that this app will accept query parameters <B>trait</B> (EFO_ID) and/or <B>gene</B>
(ENSEMBL_ID) via URL, e.g.
<B><TT>?trait=EFO_1000654</TT></B>, <B><TT>?gene=ENSG00000094914</TT></B>,
<B><TT>?trait=EFO_1000654&gene=ENSG00000094914</TT></B>.
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
<B>Authors:</B>
Jeremy Yang<SUP>1</SUP>, Stephen Mathias<SUP>1</SUP>, Cristian Bologa<SUP>1</SUP>,
Anna Waller<SUP>1</SUP>, Dhouha Grissa<SUP>2</SUP>,
Christophe Lambert<SUP>1</SUP>, David Wild<SUP>3</SUP>,
Lars Juhl Jensen<SUP>2</SUP> and Tudor Oprea<SUP>1</SUP>.<BR/>
<I><SUP>1</SUP>University of New Mexico, Translational Informatics Division, Dept. of
Internal Medicine; <SUP>2</SUP>Novo Nordisk Center for Protein Research, Copenhagen,
Denmark; <SUP>3</SUP>Indiana University, School of Informatics, Computing and Engineering, Integrative Data Science Lab.</I>
<BR/>
<B>Feedback welcome</B> to corresponding author  
<a href=\"mailto:jjyang_AT_salud_DOT_unm_DOT_edu\">Jeremy Yang</a>.<br/>
This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>")
  htm <- paste(htm, sprintf("<hr>\nBuilt with: <tt>%s; %s</tt>", R.version.string, pkgVerTxt), sep="\n")
  return(htm)
}

##########################################################################################
ui <- fluidPage(
  tags$style(".green_class {color:#00ff00} .blue_class {color:#0000ff} .red_class {color:#ff0000} .black_class {color:black}"),
  titlePanel(h2("IDG", tags$img(height="50", valign="bottom", src="IDG_logo_only.png"), APPNAME_FULL), windowTitle=APPNAME_FULL),
  fluidRow(
    column(3, 
      wellPanel(
	dqshiny::autocomplete_input("traitQry", "Trait", options=trait_menu, max_options=1000, placeholder="Query trait..."),
	dqshiny::autocomplete_input("geneQry", div("Gene"), options=as.list(c(gene_menu, filtered_gene_menu)), max_options=10000, placeholder="Query gene..."),
        	actionButton("goSubmit", label="Submit", icon=icon("cogs"), style='background-color:#EEEEEE;border-width:2px'),
        	actionButton("goReset", label="Reset", icon=icon("power-off"), style='background-color:#EEEEEE;border-width:2px')),
      wellPanel(
        sliderInput("maxHits", "MaxHits", 25, 200, 50, step=25),
	radioButtons("yAxis", "Y-Axis", choiceNames=c("OR", "nBeta", "Auto"), choiceValues=c("or_median", "n_beta", "auto"), selected="auto", inline=T),
        radioButtons("markerSizeBy", "MarkerSizeBy", choiceNames=c("N_study", "RCRAS", "None"), choiceValues=c("n_study", "rcras", NA), selected="n_study", inline=T)
      ),
	wellPanel(htmlOutput(outputId="logHtm")),
	wellPanel(htmlOutput(outputId="resultHtm"))
	),
    column(9,
	tabsetPanel(id="tabset", type="tabs",
		tabPanel(value="plot", title=textOutput("plotTabTxt"), plotlyOutput("tigaPlot", height = "500px")),
		tabPanel(value="hits", title=textOutput("hitsTabTxt"), DT::dataTableOutput("hitrows"), br(), downloadButton("hits_file", label="Download Hits")),
		tabPanel(value="detail", title="Provenance", htmlOutput("detail_summary"), DT::dataTableOutput("detail_studies")),
		tabPanel(value="traits", title="Traits (browse)", DT::dataTableOutput("traits")),
		tabPanel(value="genes", title="Genes (browse)", DT::dataTableOutput("genes")),
		tabPanel(value="studies", title="Studies (browse)", DT::dataTableOutput("studies")),
		tabPanel(value="download", title="Download",
			h1("Downloads"),
			p(downloadButton("gt_file", label="Gene-Trait Associations (all)"), textOutput("gtFileInfoTxt")),
			p(downloadButton("traits_file", label="Traits (all)"), textOutput("traitFileInfoTxt")),
			p(downloadButton("genes_file", label="Genes (all)"), textOutput("geneFileInfoTxt")),
			p(downloadButton("provenance_file", label="Provenance (association-to-study)"), textOutput("provFileInfoTxt"))
		),
		tabPanel(value="help", title="Help", htmlOutput("helpHtm"))
	))),
  hr(),
  fluidRow(
    column(12, tags$em(strong(sprintf("%s", APPNAME)), " web app from ", 
        tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
        " and ",
        tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
        " built from ",
        tags$a(href="https://www.ebi.ac.uk/gwas/", target="_blank", span("GWAS Catalog", tags$img(id="gwas_catalog_logo", height="50", valign="bottom", src="GWAS_Catalog_logo.png"))),
        " and ",
        tags$a(href="https://www.ebi.ac.uk/efo/", target="_blank", span("EFO", tags$img(id="efo_logo", height="50", valign="bottom", src="EFO_logo.png")))
        ))),
  #bsTooltip("randTraitQry", "Random query trait or gene", "right"),
  bsTooltip("goReset", "Reset.", "right"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

##########################################################################################
server <- function(input, output, session) {

  #Modal popup for genetrait details.
  observeEvent(input$showGeneTraitProvenance, {
    showModal(modalDialog(easyClose=T, footer=tagList(modalButton("Dismiss")),
	title=HTML("<H2>Gene-trait detail</H2>"),
		htmlOutput("detail_summary"), DT::dataTableOutput("detail_studies")))
  })

  output$helpHtm <- reactive({ paste(sprintf("<H2>%s Help</H2>", APPNAME), HelpHtm()) })

  Sys.sleep(1) #Needed?
  traitQryRand_count <- 0 # initialize once per session
  i_query <- 0 # initialize once per session
  
  # ?trait=EFO_0000341
  # ?trait=EFO_1000654&gene=ENSG00000094914
  httpQstr <- reactive({
    qStr <- getQueryString(session) #named list
    #if (length(qStr)>0)
    #  for (key in names(qStr))
    #    message(sprintf("DEBUG: qStr[[\"%s\"]]=\"%s\"", key, qStr[[key]]))
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
    dqshiny::update_autocomplete_input(session, "traitQry", value="")
    dqshiny::update_autocomplete_input(session, "geneQry", value="")
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
    htm <- "<h3>Provenance and details</h3>\n"
    htm <- paste(htm, sprintf("<b>TRAIT:</b> <tt>%s</tt> &harr; <b>GENE:</b> <tt>%s</tt>\n<br/>", efoId2Name(efoId_this), ensemblId2Symbol(ensemblId_this)), "\n")
    if (nrow(Hits())==0) {
      htm <- paste(htm, "<B>NO ASSOCIATIONS FOUND.</B>")
    } else {
      for (tag in names(Hits()))
        if (tag != "ok2plot")
          htm <- paste(htm, sprintf("<b>%s:</b> <tt>%s</tt>", tag, Hits()[[tag]][1]), sep=" &#8226; ")
    }
    return(htm)
  }

  processUrlParams <- function() {
    #message(sprintf("DEBUG: url: \"%s\"", urlText()))
    qStr <- httpQstr()
    if ("trait" %in% names(qStr))
      dqshiny::update_autocomplete_input(session, "traitQry", value=efoId2Name(qStr[["trait"]]))
    if ("gene" %in% names(qStr))
      dqshiny::update_autocomplete_input(session, "geneQry", value=sprintf("%s:%s", ensemblId2Symbol(qStr[["gene"]]), ensemblId2Name(qStr[["gene"]])))
  }

  # Returns both input fields as a list(gene = ***, trait = ***)
  qryIds <- reactive({
    input$goSubmit #trigger with goSubmit button?
    ids <- list(trait = NA, gene = NA)
    if (i_query==0 & length(names(httpQstr()))>0) { processUrlParams() } #1st query may be via URL http param.
    i_query <<- i_query + 1  # Must assign to up-scoped variable.
    if (!grepl("^\\s*$", input$traitQry)) ids$trait <- input$traitQry
    if (!grepl("^\\s*$", input$geneQry)) ids$gene <- input$geneQry
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
      gt_this <- gt_this[, .(efoId, trait, n_study, pvalue_mlog_median, rcras, traitMeanRank, traitMeanRankScore, study_N_mean, n_snp, n_snpw, traitNgene, or_median, n_beta)]
      setnames(gt_this, old=c("traitMeanRank", "traitMeanRankScore"), new=c("meanRank", "meanRankScore"))
      setorder(gt_this, -meanRankScore)
      gt_this[, ok2plot := as.logical(.I <= input$maxHits)]
    } else if (hitType()=="gene") {
      gt_this <- gt[efoId==qryIds()$trait]
      if (nrow(gt_this)==0) { return(NULL) }
      gt_this$TDL <- factor(gt_this$TDL, levels=c("NA", "Tdark", "Tbio", "Tchem", "Tclin"), ordered=T)
      gt_this <- gt_this[, .(ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, pvalue_mlog_median, rcras, geneMeanRank, geneMeanRankScore, study_N_mean, n_snp, n_snpw, geneNtrait, or_median, n_beta)]
      setnames(gt_this, old=c("geneMeanRank", "geneMeanRankScore"), new=c("meanRank", "meanRankScore"))
      setorder(gt_this, -meanRankScore)
      gt_this[, ok2plot := as.logical(.I <= input$maxHits)]
    } else { #hitType=="genetrait"
      gt_this <- gt[ensemblId==qryIds()$gene & efoId==qryIds()$trait]
      gt_this <- gt_this[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, pvalue_mlog_median, rcras, study_N_mean, n_snp, n_snpw, traitNgene, or_median, n_beta)]
      gt_this[, ok2plot := T]
    }
    if (hitType() %in% c("trait", "gene") & input$yAxis=="auto") { #auto-set yAxis
      n_or <- gt_this[!is.na(or_median), .N]
      n_nbeta <- gt_this[n_beta>0, .N]
      updateRadioButtons(session, "yAxis", selected=ifelse(n_nbeta>n_or, "n_beta", "or_median"))
    }
    gt_this[, or_median := ifelse(is.na(or_median), 0, or_median)] #NA encoded as zero for plotting.
    return(gt_this)
  })

  output$hitCount <- renderText({ as.character(nrow(Hits())) })
  
  output$plotTabTxt <- renderText({ ifelse(is.null(Hits()), "Plot", sprintf("Plot (%d %ss)", Hits()[(ok2plot), .N], hitType())) })
  output$hitsTabTxt <- renderText({ ifelse(!is.null(Hits()), sprintf("Hits (%d %ss)", Hits()[, .N], hitType()), "Hits (0)") })

  output$traitFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(trait_table), ncol(trait_table)) })
  output$geneFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(gene_table), ncol(gene_table)) })
  output$gtFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(gt), ncol(gt)) })
  output$provFileInfoTxt <- renderText({ sprintf("rows: %d; cols: %d", nrow(gt_prov), ncol(gt_prov)) })
  output$detail_summary <- reactive({ DetailSummaryHtm(qryIds()$trait, qryIds()$gene) })

  #Hits table has links to tiga:trait+gene, and to external resources EFO and Pharos.
  HitsWithHtm <- reactive({
    hwh <- data.table(Hits()) #copy
    if (hitType() == "trait") {
      hwh <- hwh[, efoId := sprintf("%s<a href=\"%s?trait=%s&gene=%s\"><i class=\"fa fa-search\"></i></a><a href=\"%s\" target=\"_blank\"><i class=\"fa fa-external-link\"></i></a>", efoId, urlBase(), efoId, qryIds()$gene, sapply(efoId, efoId2Uri))]
    } else if (hitType() == "gene") {
      hwh <- hwh[, geneSymbol := sprintf("%s<a href=\"%s?trait=%s&gene=%s\"><i class=\"fa fa-search\"></i></a><a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\"><i class=\"fa fa-external-link\"></i></a>", geneSymbol, urlBase(), qryIds()$trait, ensemblId, geneSymbol)]
    }
    return(hwh)
  })

  #Here we hide/show/select tabs.
  output$resultHtm <- reactive({
    message(sprintf("TraitQuery: \"%s\" (%s)", traitQryName(), qryIds()$trait))
    message(sprintf("GeneQuery: \"%s\" (%s)", geneQryName(), qryIds()$gene))
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      hideTab("tabset", "plot")
      hideTab("tabset", "hits")
      hideTab("tabset", "detail")
      updateTabsetPanel(session, "tabset", selected="traits")
      return("No query. Search? Browse?")
    }
    htm <- sprintf("<B>Results:</B>")
    if (qryType() %in% c("trait", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", traitQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", efoId2Uri(qryIds()$trait), qryIds()$trait))
    }
    if (qryType() %in% c("gene", "genetrait")) {
      htm <- paste0(htm, sprintf("\"%s\"", geneQryName()))
      htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"https://pharos.nih.gov/targets/%s\">%s</a>)", qryIds()$gene, qryIds()$gene))
    }
    if (qryType() %in% c("gene", "trait")) {
      showTab("tabset", "plot")
      showTab("tabset", "hits")
      hideTab("tabset", "detail")
      updateTabsetPanel(session, "tabset", selected="plot")
    } else if (qryType()=="genetrait") {
      hideTab("tabset", "plot")
      hideTab("tabset", "hits")
      showTab("tabset", "detail")
      updateTabsetPanel(session, "tabset", selected="detail")
    }
    if (!is.null(Hits())) {
      htm <- paste0(htm, sprintf("; N_%s: %d plotted (%d total)", hitType(), Hits()[(ok2plot), .N], Hits()[, .N]))
      if (!is.null(Hits()[["or_median"]]) & !is.null(Hits()[["n_beta"]]))
        htm <- paste0(htm, sprintf("; ORs: %d; N_betas>0: %d", Hits()[or_median>0, .N], Hits()[n_beta>0, .N]))
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
  
  #Errors here?
  markerSize <- reactive({
    if (input$markerSizeBy=="n_study") {
      size <- 5*Hits()[(ok2plot), n_study]
    } else if (input$markerSizeBy=="rcras") {
      size <- 5*Hits()[(ok2plot), rcras]
    } else { #NA
      size <- rep(10, nrow(Hits()[(ok2plot)]))
    }
    size <- pmax(size, rep(10, nrow(Hits()[(ok2plot)]))) #min
    size <- pmin(size, rep(80, nrow(Hits()[(ok2plot)]))) #max
    message(sprintf("DEBUG: markerSizeBy: %s", as.character(input$markerSizeBy)))
    message(sprintf("DEBUG: nrow(Hits()[(ok2plot)]): %d", nrow(Hits()[(ok2plot)])))
    message(sprintf("DEBUG: length(markerSize): %d", length(size)))
    return(size)
  })

  markerTextGenes <- reactive({
    text <- paste0(
    paste0("<b>", Hits()[(ok2plot), geneSymbol], "</b> (", Hits()[(ok2plot), ensemblId], ")"), 
    paste0("<br><b>", Hits()[(ok2plot), geneName], "</b>"),
    paste0("<br>Fam:", Hits()[(ok2plot), geneFamily]),
    paste0(", TDL:", Hits()[(ok2plot), TDL]),
    paste0("; N_trait = ", Hits()[(ok2plot), geneNtrait]),
    ";<br>meanRank = ", round(Hits()[(ok2plot), meanRank], digits=3),
    "; meanRankScore = ", round(Hits()[(ok2plot), meanRankScore], digits=3),
    ";<br>N_study = ", Hits()[(ok2plot), n_study],
    "; study_N = ", Hits()[(ok2plot), study_N_mean], 
    "; N_snp = ", Hits()[(ok2plot), n_snp],
    ";<br>N_snpw = ", Hits()[(ok2plot), n_snpw],
    "; OR = ", round(Hits()[(ok2plot), or_median], digits=2), 
    "; N_beta = ", Hits()[(ok2plot), n_beta],
    "; pVal = ", sprintf("%.2g", 10^(-Hits()[(ok2plot), pvalue_mlog_median])),
    "; RCRAS = ", round(Hits()[(ok2plot), rcras], digits=2),
    ";<br>DEBUG: .I = ", Hits()[(ok2plot), .I]
      )
    message(sprintf("DEBUG: length(markerTextGenes): %d", length(text)))
    return(text)
  })
  markerTextTraits <- reactive({
    text <- paste0(
    paste0("<b>", Hits()[(ok2plot), efoId], "</b>"), 
    paste0("<br><b>", Hits()[(ok2plot), trait], "</b>"),
    paste0("; N_gene = ", Hits()[(ok2plot), traitNgene]),
    ";<br>meanRank = ", round(Hits()[(ok2plot), meanRank], digits=3),
    "; meanRankScore = ", round(Hits()[(ok2plot), meanRankScore], digits=3),
    ";<br>N_study = ", Hits()[(ok2plot), n_study],
    "; study_N = ", Hits()[(ok2plot), study_N_mean], 
    "; N_snp = ", Hits()[(ok2plot), n_snp],
    ";<br>N_snpw = ", Hits()[(ok2plot), n_snpw],
    "; OR = ", round(Hits()[(ok2plot), or_median], digits=2), 
    "; N_beta = ", Hits()[(ok2plot), n_beta],
    "; pVal = ", sprintf("%.2g", 10^(-Hits()[(ok2plot), pvalue_mlog_median])),
    "; RCRAS = ", round(Hits()[(ok2plot), rcras], digits=2))
    return(text)
  })

  output$tigaPlot <- renderPlotly({
    xaxis <- list(title="Evidence (meanRankScore)", type="normal", zeroline=F, showline=F)
    yaxis <- list(title=ifelse(input$yAxis=="n_beta", "N_Beta", "Effect (OddsRatio)"), type="normal")
    axis_none <- list(zeroline=F, showline=F, showgrid=F, showticklabels=F)
    if (is.na(qryIds()$trait) & is.na(qryIds()$gene)) {
      title <- "<I>(No query.)</I>"
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (is.na(qryIds()$trait) & qryIds()$gene %in% filtered$id) {
      title <- sprintf("<I>(No hits.)</I><br>Gene %s:<br>%s<br><I>Filtered by TIGA preprocessing<br>Reason: %s</I>", qryIds()$gene, geneQryName(), filtered[id==qryIds()$gene, reason])
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    } else if (is.null(Hits())) {
      title <- "<I>(No hits.)</I>"
      return(plot_ly(type="scatter", mode="marker") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20)))
    }
    
    #if (qryType()=="trait")
    #  message(sprintf("DEBUG: %s", paste(collapse=",", paste(Hits()[(ok2plot)]$geneSymbol, as.character(markerSize()), sep=":"))))
    #else if (qryType()=="gene")
    #  message(sprintf("DEBUG: %s", paste(collapse=",", paste(Hits()[(ok2plot)]$efoId, as.character(markerSize()), sep=":"))))
    
    if (hitType()=="gene") {
      if (input$yAxis=="n_beta")
        p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok2plot)], 
          #x=~meanRankScore, 
          x = 100*jitter(Hits()[(ok2plot), meanRankScore], 10),
          #y=~n_beta,
          y = jitter(Hits()[(ok2plot), n_beta], 10),
                marker=list(symbol="circle", size=markerSize()), text=markerTextGenes(),
                color=~TDL, colors=c("gray", "black", "red", "green", "blue"))
      else # or_median or auto
        p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok2plot)], 
          #x=~meanRankScore, 
          x = 100*jitter(Hits()[(ok2plot), meanRankScore], 10),
          #y=~or_median,
          y = jitter(Hits()[(ok2plot), or_median], 10),
                marker=list(symbol="circle", size=markerSize()), text=markerTextGenes(),
                color=~TDL, colors=c("gray", "black", "red", "green", "blue"))
      p <- config(p, displayModeBar=F) %>% layout(xaxis=xaxis, yaxis=yaxis, 
          title=paste0(traitQryName(), "<br>", "(", qryType(), ":", qryIds()$trait, ")"),
          margin=list(t=80,r=50,b=60,l=60), showlegend=T,
	  legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
          font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits()[(ok2plot)]), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else if (hitType()=="trait") {
      if (input$yAxis=="n_beta")
        p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok2plot)], x=~meanRankScore, y=~n_beta,
                marker=list(symbol="circle", size=markerSize()), text=markerTextTraits())
      else # or_median or auto
        p <- plot_ly(type='scatter', mode='markers', data=Hits()[(ok2plot)], x=~meanRankScore, y=~or_median,
                marker=list(symbol="circle", size=markerSize()), text=markerTextTraits())
      p <- config(p, displayModeBar=F) %>% layout(xaxis=xaxis, yaxis=yaxis, 
          title=paste0(geneQryName(), "<br>", "(", qryType(), ":", qryIds()$gene, ")"),
          margin=list(t=80,r=50,b=60,l=60), showlegend=F,
	  legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
          font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(Hits()), "; ", nrow(Hits()[(ok2plot)]), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    } else { # "genetrait"
      title <- "<I>(No plot in gene+trait mode.)</I>"
      p <- plot_ly(type="scatter", mode="markers") %>% config(displayModeBar=F) %>% layout(title=title, xaxis=axis_none, yaxis=axis_none, margin=list(t=120,b=20))
    }
    return(p)
  })

  #DT numbers cols from 0.
  #"ensemblId","geneSymbol","geneName","geneFamily","TDL","n_study","pvalue_mlog_median","rcras","meanRank","meanRankScore","study_N_mean","n_snp","n_snpw","geneNtrait","or_median","n_beta"
  output$hitrows <- DT::renderDataTable({
    if (is.null(Hits())) return(NULL)
    if (hitType()=="gene") {
      return(DT::datatable(data=HitsWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap",
          	selection=list(target="row", mode="multiple", selected=NULL),
		colnames=c("ENSG", "GSYMB", "GeneName", "idgFam", "idgTDL", "N_study", "pVal_mlog", "RCRAS", "meanRank", "meanRankScore", "study_N", "N_snp", "N_snpw", "N_trait", "OR", "N_beta", "ok2plot"),
          	options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
          	               list(className='dt-center', targets=c(0, 1, 3:(ncol(HitsWithHtm())-2))),
          	               list(visible=F, targets=c(0, ncol(HitsWithHtm())-1)) #Hide EnsemblId
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_median", "or_median", "rcras", "n_snpw", "meanRank", "meanRankScore"), digits=2)
    %>% DT::formatStyle(c("n_study", "pvalue_mlog_median", "rcras"), backgroundColor="skyblue", fontWeight="bold")
    %>% DT::formatStyle(c("meanRank", "meanRankScore"), color="black", backgroundColor="orange", fontWeight="bold")
    %>% DT::formatStyle("TDL", backgroundColor=styleEqual(c("Tclin", "Tchem", "Tbio", "Tdark"), c("#4444DD", "#11EE11", "#EE1111", "gray"))
  )
        )
  } else if (hitType()=="trait") {
      return(DT::datatable(data=HitsWithHtm(), escape=F, rownames=F, class="cell-border stripe", style="bootstrap",
		selection=list(target="row", mode="multiple", selected=NULL),
		colnames=c("efoId", "trait", "N_study", "pVal_mlog", "RCRAS", "meanRank", "meanRankScore", "study_N", "N_snp", "N_snpw", "N_gene", "OR", "N_beta", "ok2plot"),
		options=list(
			autoWidth=T, dom='tip',
			columnDefs=list(
				list(className='dt-center', targets=c(0, 1, 3:(ncol(HitsWithHtm())-2))),
				list(visible=F, targets=c(ncol(HitsWithHtm())-1))
				)
			)
	) %>% DT::formatRound(columns=c("pvalue_mlog_median", "or_median", "rcras", "n_snpw", "meanRank", "meanRankScore"), digits=2)
    %>% DT::formatStyle(c("n_study", "pvalue_mlog_median", "rcras"), backgroundColor="skyblue", fontWeight="bold")
    %>% DT::formatStyle(c("meanRank", "meanRankScore"), color="black", backgroundColor="orange", fontWeight="bold")
	) 
  }
  }, server=T)

  #All-traits table has tiga-trait links.
  trait_tableHtm <- reactive({
    dt <- data.table(trait_table) #copy
    dt[, idHtm := sprintf("<a href=\"%s?trait=%s\">%s</a>", urlBase(), efoId, efoId)]
    dt[, .(efoId = idHtm, trait, N_study, N_gene)][order(trait)]
  })

  output$traits <- DT::renderDataTable({
    DT::datatable(data=trait_tableHtm()[, .(efoId, trait, N_study, N_gene)], rownames=F, width="100%", options=list(autoWidth=F, dom='tipf'), escape=F)
  }, server=T)

  #All-genes table has tiga-gene links.
  gene_tableHtm <- reactive({
    dt <- data.table(gene_table)[order(geneSymbol)]
    dt[, symbHtm := sprintf("<a href=\"%s?gene=%s\">%s</a>", urlBase(), ensemblId, geneSymbol)]
    dt[, .(ensemblId, geneSymbol = symbHtm, geneName, geneFamily, TDL, N_study, N_trait, filtered)]
  })

  output$genes <- DT::renderDataTable({
    DT::datatable(data=gene_tableHtm()[, .(geneSymbol, geneName, geneFamily, TDL, N_study, N_trait, filtered)], rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)
  
  study_tableHtm <- reactive({
    dt <- data.table(study_table)
    dt[, gcHtm := sprintf("<a href=\"https://www.ebi.ac.uk/gwas/studies/%s\">%s</a><i class=\"fa fa-external-link\">", STUDY_ACCESSION, STUDY_ACCESSION)]
    dt[, pubmedHtm := sprintf("<a href=\"https://pubmed.ncbi.nlm.nih.gov/%s\">%s</a><i class=\"fa fa-external-link\">", PUBMEDID, PUBMEDID)]
    dt[, .(Accession=gcHtm, Study=STUDY, PMID=pubmedHtm, DatePublished=DATE_PUBLISHED, DateAdded=DATE_ADDED_TO_CATALOG)][order(-DatePublished)]
  })

  output$studies <- DT::renderDataTable({
    DT::datatable(data=study_tableHtm()[, .(Accession, Study, PMID, DatePublished, DateAdded)], rownames=F, width="100%", options=list(autoWidth=F, dom='tipf'), escape=F)
  }, server=T)

  output$detail_studies <- DT::renderDataTable({
    DT::datatable(data=detail_studies_tableHtm(), rownames=F, width="100%", options=list(autoWidth=F, dom='tip'), escape=F,
      caption = tags$caption(style='caption-side:top; text-align:center;', 'Table: ', tags$b('Studies with association evidence')))
  }, server=T)

  detail_studies_tableHtm <- reactive({
    dt <- data.table(merge(gt_prov[efoId == qryIds()$trait & ensemblId == qryIds()$gene, .(STUDY_ACCESSION)], study_table[, .(STUDY_ACCESSION, STUDY, PUBMEDID, DATE_PUBLISHED, DATE_ADDED_TO_CATALOG)], by="STUDY_ACCESSION", all.x=T, all.y=F))
    dt <- unique(dt)
    dt[, gcHtm := sprintf("<a href=\"https://www.ebi.ac.uk/gwas/studies/%s\">%s</a>", STUDY_ACCESSION, STUDY_ACCESSION)]
    dt[, pubmedHtm := sprintf("<a href=\"https://pubmed.ncbi.nlm.nih.gov/%s\">%s</a>", PUBMEDID, PUBMEDID)]
    dt[, .(Accession=gcHtm, Study=STUDY, PMID=pubmedHtm, DatePublished=DATE_PUBLISHED, DateAdded=DATE_ADDED_TO_CATALOG)][order(-DatePublished)]
  })

  Hits_export <- reactive({
    if (is.null(Hits())) { return(NULL) }
    hits_out <- data.table(Hits()) #copy
    if (hitType()=="gene") {
      hits_out[["efoId"]] <- qryIds()$trait
      hits_out[["trait"]] <- traitQryName()
      hits_out <- hits_out[, .(efoId, trait, ensemblId, geneSymbol, geneName, geneFamily, TDL, n_study, pvalue_mlog_median, rcras, meanRank, meanRankScore, study_N_mean, n_snp, n_snpw, geneNtrait, or_median, n_beta)]
    } else if (hitType()=="trait") {
      hits_out[["ensemblId"]] <- qryIds()$gene
      hits_out[["geneName"]] <- geneQryName()
      hits_out <- hits_out[, .(ensemblId, geneName, efoId, trait, n_study, pvalue_mlog_median, rcras, meanRank, meanRankScore, study_N_mean, n_snp, n_snpw, traitNgene, or_median, n_beta)]
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
      write_delim(gene_table, file, delim="\t")
    }
  )
  output$gt_file <- downloadHandler(
    filename = "tiga_gene-trait_stats.tsv",
    content = function(file) {
      write_delim(gt, file, delim="\t")
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

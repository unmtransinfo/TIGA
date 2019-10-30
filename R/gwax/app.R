##########################################################################################
### GWAX: GWAS Explorer
### gt = gene-trait data
### Dataset gt_stats.tsv from gwax_gt_stats.R.
### Jeremy Yang
##########################################################################################
library(readr)
library(data.table)
library(shiny, quietly=T)
library(shinyBS, quietly=T)
library(DT, quietly=T)
library(plotly, quietly=T)

##########################################################################################
APPNAME <- "GWAX"
MIN_ASSN <- 5
#
t0 <- proc.time()
DEBUG <- F
if (!file.exists("gwax.Rdata") | DEBUG) {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  gt <- read_delim("gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), 
	n_study=col_integer(), n_snp=col_integer(), n_snpw=col_double(),
	geneNtrait=col_integer(), geneNstudy=col_integer(),
	traitNgene=col_integer(), traitNstudy=col_integer(), pvalue_mlog_median=col_double(), or_median=col_double(), study_N_mean=col_double(), rcras=col_double(),
	muScore=col_double(), muRank=col_integer()))
  setDT(gt)
  setnames(gt, old=c("geneSymbol", "geneIdgTdl"), new=c("gsymb", "TDL"))
  #
  trait_table <- gt[, .(N_gene = uniqueN(ensemblId)),  by=c("trait_uri", "trait", "traitNstudy")]
  setnames(trait_table, old=c("traitNstudy"), new=c("N_study"))
  trait_table <- trait_table[N_gene>=MIN_ASSN & N_study>0]
  trait_table[['id']] <- as.factor(sub("^.*/","", trait_table$trait_uri))
  trait_table <- trait_table[, .(trait_uri, id, trait, N_study, N_gene)]
  
  trait_menu <- sub("^.*/", "", trait_table$trait_uri) #named vector
  names(trait_menu) <- trait_table$trait
  trait_menu <- trait_menu[order(names(trait_menu))]
  #
  gene_table <- gt[, .(gsymb, geneName, geneFamily, TDL, N_study = geneNstudy, N_trait = uniqueN(trait_uri)),  by=c("ensemblId")]
  gene_table <- unique(gene_table)
  #
  save(gt, trait_table, trait_menu, gene_table, file="gwax.Rdata")
} else {
  message(sprintf("Loading gwax.Rdata..."))
  load("gwax.Rdata")
  setDT(gt)
}
#
t_elapsed <- (proc.time()-t0)[3]
#
message(sprintf("Gene count, IDs: %d; symbols: %d", uniqueN(gt$ensemblId), uniqueN(gt$gsymb)))
message(sprintf("Trait count (total): %d", uniqueN(gt$trait)))
#trait_table[['ontology']] <- as.factor(sub("_.*$","", trait_table$id))
#trait_counts <- trait_table[, .N, by="ontology"][order(-N)]
#message(sprintf("traits (%10s): %4d / %4d (%4.1f%%)\n", trait_counts$ontology, trait_counts$N, sum(trait_counts$N), 100*trait_counts$N/sum(trait_counts$N)))
message(sprintf("Trait count (filtered; n_assn<%d): %d", MIN_ASSN, uniqueN(gt$trait)-length(trait_menu)))
message(sprintf("Trait count (menu; n_assn>=%d): %d", MIN_ASSN, length(trait_menu)))
message(sprintf("DEBUG: COUNT or_median: %d", sum(!is.na(gt$or_median))))
message(sprintf("DEBUG: COUNT pvalue_mlog_median: %d", sum(!is.na(gt$pvalue_mlog_median))))
message(sprintf("DEBUG: COUNT rcras: %d", sum(!is.na(gt$rcras))))
#
db_htm <- sprintf("<B>Dataset:</B> genes: %d ; traits: %d; top_traits: %d (t_load: %.1fs)",
	uniqueN(gt$gsymb), uniqueN(gt$trait), length(trait_menu), t_elapsed)
###
TDLS <- c("Tclin", "Tchem", "Tbio", "Tdark")
TDLS_Names <- list(
span("Tclin", icon("circle", lib="font-awesome", class="blue_class")),
span("Tchem", icon("circle", lib="font-awesome", class="green_class")),
span("Tbio", icon("circle", lib="font-awesome", class="red_class")),
span("Tdark", icon("circle", lib="font-awesome", class="black_class")))

idgfams <- c("GPCR", "Kinase", "IC", "NR", "Other")
axes <- c("Effect", "Evidence")
#
qryTraitRand <- sample(trait_menu, 1)
message(sprintf("DEBUG: qryTraitRand: %s", qryTraitRand))
#
#############################################################################
HelpHtm <- function() {(
  sprintf("<P><B>GWAX</B>, GWAS Explorer, is designed to facilitate drug target illumination by 
scoring and ranking of protein-coding genes associated with traits from genome-wide association studies
(GWAS). Rather than a comprehensive analysis of GWAS for all biological implications and insights, this
more focused application provides a rational method by which GWAS findings can be 
aggregated and filtered for applicable, actionable intelligence, 
evidence usable by drug discovery scientists to enrich prioritization of target hypotheses.
<P>
<B>Dataset:</B>
Data from <A HREF=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">The NHGRI-EBI GWAS Catalog</A>.
<UL>
<LI>Well studied traits only available via this app, with minimum associations
(MIN_ASSN=%d, may be multiple for each gene). Traits are mapped to EFO, HPO,
Orphanet, PATO or GO, with 91%% mapped to EFO.
<LI>Mapped genes via Ensembl pipeline as per GWAS Catalog documentation. Reported genes ignored for consistency and accountable confidence assessment in this app and downstream.
<LI>In this version, OR required and BETA ignored, due to problem of harmonizing varying BETA units.
</UL>
<B>Datatypes:</B>
<UL>
  <LI><B>pVal_mLog<SUP>*</SUP></B>: median(-Log(pValue)) supporting trait-gene association.
  <LI><B>OR<SUP>*</SUP></B>: median(odds ratio, inverted if &lt;1) supporting trait-gene association.
  <LI><B>N_snp<SUP>*</SUP></B>: SNPs involved with trait-gene association.
  <LI><B>N_snpw<SUP>*</SUP></B>: N_snp weighted by distance inverse exponential.
  <LI><B>N_study<SUP>*</SUP></B>: studies supporting trait-gene association.
  <LI><B>study_N<SUP>*</SUP></B>: mean(SAMPLE_SIZE) supporting trait-gene association.
  <LI><B>RCRAS<SUP>*</SUP></B>: Relative Citation Ratio (RCR) Aggregated Score (iCite-RCR-based)
  <LI><B>geneNtrait<SUP>**</SUP></B>: total traits associated with gene.
  <LI><B>traitNgene<SUP>**</SUP></B>: total genes associated with trait.
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
<LI>Plot markers may be sized by <B>N_study</B> or <B>RCRAS</B>.
<LI>Note that this app will accept query parameter <B>efo_id</B> via URL, e.g.
<B><TT>?efo_id=EFO_0000341</TT></B>.
</UL>
<BR/>
<B>More documentation:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\">GWAS Catalog</a>
<LI><a href=\"https://www.ebi.ac.uk/gwas/docs/fileheaders\">GWAS Catalog data dictionary</a>
<LI><a href=\"https://www.ebi.ac.uk/efo/\">Experimental Factor Ontology (EFO)</a>
</UL>
<BR/>
<B>Authors:</B> Jeremy Yang<SUP>1</SUP>, Stephen Mathias<SUP>1</SUP>, Cristian
Bologa<SUP>1</SUP>, Lars Juhl Jensen<SUP>2</SUP>, Christophe Lambert<SUP>1</SUP>, David Wild<SUP>3</SUP> and Tudor
Oprea<SUP>1</SUP>.<BR/>
<I><SUP>1</SUP>University of New Mexico, Translational Informatics Division, Dept. of
Internal Medicine; <SUP>2</SUP>Novo Nordisk Center for Protein Research, Copenhagen,
Denmark; <SUP>3</SUP>Indiana University, School of Informatics, Computing and Engineering, Integrative Data Science Lab.</I>
<BR/>
<B>Correspondence</B> from users of this app is welcome, and should be directed to 
<a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
  Built with R-Shiny &amp; Plotly.<BR/>
  This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>",
	MIN_ASSN)
)}

##########################################################################################
ui <- fluidPage(
  tags$style(".green_class {color:#00ff00} .blue_class {color:#0000ff} .red_class {color:#ff0000} .black_class {color:black}"),
  titlePanel(h2("IDG", tags$img(height="50", valign="bottom", src="IDG_logo_only.png"), sprintf("%s: GWAS Explorer", APPNAME),tags$img(height="40", valign="bottom", src="GWAS_Catalog_logo.png"), span(style="font-size:18px", "GWAS Catalog-based drug target illumination (BETA)")),
      windowTitle=APPNAME),
  fluidRow(
    column(3, 
      wellPanel(
	selectInput("traitQry", "Query trait", choices=trait_menu, selectize=T, selected=qryTraitRand),
        sliderInput("maxHits", "Max_hits", 25, 200, 50, step=25),
	checkboxGroupInput("logaxes", "LogAxes", choices=axes, selected=NULL, inline=T),
        radioButtons("markerSizeBy", "MarkerSizeBy", choiceNames=c("N_study", "RCRAS"), choiceValues=c("n_study", "rcras"), selected="n_study", inline=T),
        br(),
	actionButton("goRefresh", "Refresh", style='padding:2px;background-color:#DDDDDD;font-weight:bold'),
	actionButton("showHelp", "Help", style='padding:2px;background-color:#DDDDDD; font-weight:bold'),
	actionButton("randQuery", "RandomQuery", style='padding:2px; background-color:#DDDDDD; font-weight:bold')
      ),
	wellPanel(htmlOutput(outputId="log_htm"))
	),
    column(9,
	tabsetPanel(type="tabs",
		tabPanel("Plot", plotlyOutput("traitGenesPlot", height = "600px")),
		tabPanel("Hits", DT::dataTableOutput("generows"), br(), downloadButton("hits_file", label="Download Hits")),
		tabPanel("Traits (All)", DT::dataTableOutput("traits")),
		tabPanel("Genes (All)", DT::dataTableOutput("genes"))
	))),
  fluidRow(column(12, wellPanel(htmlOutput(outputId="result_htm", height="60px")))),
#  fluidRow(column(12, wellPanel(htmlOutput(outputId="log_htm", height="60px")))),
  fluidRow(
    column(12, tags$em(strong(sprintf("%s", APPNAME)), " web app from ", 
        tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
        " and ",
        tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
        " data from ",
        tags$a(href="https://www.ebi.ac.uk/gwas/", target="_blank", span("GWAS Catalog", tags$img(id="gwas_catalog_logo", height="50", valign="bottom", src="GWAS_Catalog_logo.png"))),
        " with ",
        tags$a(href="https://www.ebi.ac.uk/efo/", target="_blank", span("EFO", tags$img(id="efo_logo", height="50", valign="bottom", src="EFO_logo.png")))
        ))),
  bsTooltip("traitQry", "Select from most studied traits.", "top"),
  bsTooltip("maxHits", "Max hits cutoff, filtered by muScore.", "top"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

id2uri <- function(id) {
  if (is.null(id)) { return(NA) }
  if (grepl("^EFO_", id)) {
    return(sprintf("http://www.ebi.ac.uk/efo/%s", id))
  } else if (grepl("(^HP_|^GO_|^CHEBI_)", id)) {
    return(sprintf("http://purl.obolibrary.org/obo/%s", id))
  } else if (grepl("^Orphanet_", id)) {
    return(sprintf("http://www.orpha.net/ORDO/%s", id))
  } else {
    return(NA)
  }
}

##########################################################################################
server <- function(input, output, session) {
  observeEvent(input$showHelp, {
    showModal(modalDialog(easyClose=T, footer=tagList(modalButton("Dismiss")),
	title=HTML(sprintf("<H2>%s Help</H2>", APPNAME)),
	HTML(HelpHtm())))
  })

  Sys.sleep(1)
  qryTraitRand_previous <- 0 # initialize once per session
  i_query <- 0 # initialize once per session
  
  # ?efo_id=EFO_0000341
  httpQstr <- reactive({
    qStr <- getQueryString(session) #named list
    if (length(qStr)>0) {
      message(sprintf("DEBUG: QueryString: \"%s\"", paste(names(qStr), qStr, sep="=", collapse="&")))
      for (key in names(qStr)) {
        message(sprintf("DEBUG: qStr[[\"%s\"]] = \"%s\"", key, qStr[[key]]))
      }
    }
    return(qStr)
  })

  urlText <- reactive({
    sprintf("%s//%s:%s%s%s",
      session$clientData$url_protocol, session$clientData$url_hostname,
      session$clientData$url_port, session$clientData$url_pathname,
      session$clientData$url_search
    )
  })
 
  query_trait_uri <- reactive({
    input$goRefresh # Re-run this and downstream on action button.
    if (input$randQuery>qryTraitRand_previous) {
      qryTraitRand_previous <<- input$randQuery # Must assign to up-scoped variable.
      qryTraitRand <- sample(trait_menu, 1)
      message(sprintf("DEBUG: qryTraitRand_previous=%d; qryTraitRand=%s", qryTraitRand_previous, qryTraitRand))
      updateTextInput(session, "traitQry", value=as.character(qryTraitRand)) #Better than updateSelectizeInput ??
    }
    qStr <- httpQstr()
    if (i_query==0 & "efo_id" %in% names(qStr)) {
      message(sprintf("DEBUG: qStr[[\"efo_id\"]]=%s", qStr[["efo_id"]]))
      updateTextInput(session, "traitQry", value=as.character(qStr[["efo_id"]]))
    }
    i_query <<- i_query + 1
    if (input$traitQry=="") { return(NULL) }
    return(id2uri(input$traitQry))
  })
  query_trait_name <- reactive({
    if (is.null(query_trait_uri())) { return(NULL) }
    if (!(query_trait_uri() %in% gt$trait_uri)) { return(NULL) }
    return(gt[trait_uri==query_trait_uri(), trait][1])
  })
  query_trait_id <- reactive({
    if (is.null(query_trait_uri())) { return(NULL) }
    return(sub("^.*/","", query_trait_uri()))
  })
  
  gHits <- reactive({
    if (is.null(query_trait_uri())) { return(NULL) }
    gt_this <- gt[trait_uri==query_trait_uri()]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this$TDL <- factor(gt_this$TDL, levels=c("NA", "Tdark", "Tbio", "Tchem", "Tclin"), ordered=T)
    gt_this <- gt_this[, .(ensemblId, gsymb, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, geneNtrait, pvalue_mlog_median, or_median, rcras, muScore, muRank)]
    message(sprintf("DEBUG: gHits() COUNT pvalue_mlog_median: %d", sum(!is.na(gt_this$pvalue_mlog_median))))
    message(sprintf("DEBUG: maxHits: %d", input$maxHits))
    gt_this[, ok := as.logical(muRank<=input$maxHits)]
    setorder(gt_this, muRank)
    message(sprintf("DEBUG: nrow(gt_this): %d; n_ok: %d", nrow(gt_this), sum(gt_this$ok)))
    return(gt_this)
  })
  
  gHits_htm <- reactive({
    hh <- gHits()
    hh <- hh[, gsymb := sprintf("<a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\">%s</a>", gsymb, gsymb)]
    return(hh)
  })
  
  output$result_htm <- reactive({
    htm <- sprintf("<B>Results:</B>")
    htm <- paste0(htm, sprintf("\"%s\"", query_trait_name()))
    htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", id2uri(query_trait_id()), query_trait_id()))
    message(sprintf("Query: \"%s\" (%s)", query_trait_name(), query_trait_id()))
    if (!is.null(gHits())) {
      htm <- paste0(htm, sprintf("; found N_genes: %d", nrow(gHits())))
    } else {
      htm <- paste0(htm, sprintf("; No genes found."))
    }
    return(htm)
  })

  output$log_htm <- reactive({
    htm <- db_htm
    message(sprintf("DEBUG: url = \"%s\"", urlText()))
    return(htm)
  })
  
  markerSize <- reactive({
    if (input$markerSizeBy=="n_study") {
      message(sprintf("DEBUG: markerSizeBy: %s", input$markerSizeBy))
      size <- 10*gHits()[(ok), n_study]
    } else if (input$markerSizeBy=="rcras") {
      message(sprintf("DEBUG: markerSizeBy: %s", input$markerSizeBy))
      size <- 10*gHits()[(ok), rcras]
    } else {
      message(sprintf("DEBUG: markerSizeBy: %s", input$markerSizeBy))
      size <- rep(10, nrow(gHits()[(ok)]))
    }
    message(sprintf("DEBUG: nrow(gHits()): %d", nrow(gHits()[(ok)])))
    message(sprintf("DEBUG: length(markerSize): %d", length(size)))
    size <- pmax(size, rep(10, nrow(gHits()[(ok)]))) #min
    size <- pmin(size, rep(80, nrow(gHits()[(ok)]))) #max
    message(sprintf("DEBUG: length(markerSize): %d", length(size)))
    return(size)
  })

  markerText <- reactive({
    text=paste0(
    "<b>", gHits()[(ok)]$gsymb, "</b> (", gHits()[(ok)]$ensemblId, ")", 
    "<br><b>", gHits()[(ok)]$geneName, "</b>",
    "<br>Fam:", gHits()[(ok)]$geneFamily, ", TDL:", gHits()[(ok)]$TDL,
    ";<br>muScore = ", gHits()[(ok)]$muScore,
    "; muRank = ", gHits()[(ok)]$muRank,
    ";<br>N_study = ", gHits()[(ok)]$n_study,
    "; study_N = ", gHits()[(ok)]$study_N_mean, 
    "; N_trait = ", gHits()[(ok)]$geneNtrait,
    "; N_snp = ", gHits()[(ok)]$n_snp,
    ";<br>N_snpw = ", gHits()[(ok)]$n_snpw,
    "; OR = ", round(gHits()[(ok)]$or_median, digits=2), 
    "; pVal = ", sprintf("%.2g", 10^(-gHits()[(ok)]$pvalue_mlog_median)),
    "; RCRAS = ", round(gHits()[(ok)]$rcras, digits=2))
    #message(sprintf("DEBUG: length(text): %d", length(text)))
    return(text)
  })

  output$traitGenesPlot <- renderPlotly({
    if (is.null(gHits())) { return(NULL) }
    #xrange <- NULL
    #xaxis <- list(title="Evidence (muScore)", range=xrange, type=ifelse("Evidence" %in% input$logaxes, "log", "normal"))
    xaxis <- list(title="Evidence (muScore)", type="normal", zeroline=F, showline=F)
    yaxis <- list(title="Effect (OddsRatio)", type=ifelse("Effect" %in% input$logaxes, "log", "normal"))

    message(sprintf("DEBUG: %s", paste(collapse=",", paste(gHits()[(ok)]$gsymb, as.character(markerSize()), sep=":"))))
    
    p <- plot_ly(type='scatter', mode='markers', data=gHits()[(ok)],
        x=(~muScore + rnorm(nrow(gHits()[(ok)]), sd=(.01*(max(gHits()$muScore))))), #Custom jitter
        y=~or_median,
        color=~TDL, colors=c("gray", "black", "red", "green", "blue"),
        marker=list(symbol="circle", size=markerSize()),
        text=markerText()
      ) %>% config(displayModeBar=F) %>%
      layout(xaxis=xaxis, yaxis=yaxis, 
        title=paste0(query_trait_name(), "<br>", "(", input$traitQry, ")"),
        margin=list(t=80,r=50,b=60,l=60), showlegend=T,
	legend=list(x=1, y=1, traceorder="normal", orientation="h", xanchor="right", yanchor="auto", itemsizing="constant", borderwidth=1, bordercolor="gray"),
        font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N: ", nrow(gHits()), "; ", nrow(gHits()[(ok)]), " shown)"), showarrow=F, x=0, y=1, xref="paper", yref="paper")
    return(p)
  })
  
  #"ensemblId","gsymb","geneName","geneFamily","TDL","n_study","study_N_mean","n_snp","n_snpw","geneNtrait","pvalue_mlog_median","or_median","rcras","muScore","muRank"
  output$generows <- DT::renderDataTable({
    if (is.null(gHits())) { return(NULL) }
    DT::datatable(data=gHits_htm(), escape=F, rownames=F,
	selection=list(target="row", mode="multiple", selected=NULL),
	class="cell-border stripe", style="bootstrap",
	options=list(autoWidth=T, dom='tip', #dom=[lftipr]
		columnDefs = list(
			list(className='dt-center', targets=c(0, 3:(ncol(gHits())-2))), #Numbered from 0.
			list(visible=F, targets=c(0, ncol(gHits())-1)))), #Numbered from 0.
	colnames=c("ENSG","GSYMB","GeneName","idgFam","idgTDL","N_study","study_N","N_snp","N_snpw","N_trait","pVal_mlog","OR","RCRAS",
	           "muScore", "muRank", "ok")
  ) %>% formatRound(digits=c(0,2,2,2,2), columns=c(7,9,11,12,13)) #Numbered from 1.
  }, server=T)

  trait_table_htm <- reactive({
    tt <- trait_table[, id_htm := sprintf("<a href=\"%s\" target=\"_blank\">%s</a>", trait_uri, id)][, .(ID = id_htm, trait, N_study, N_gene)][order(trait)]
    return(tt)
  })
  
  gene_table_htm <- reactive({
    gg <- gene_table[, symb_htm := sprintf("<a href=\"https://pharos.nih.gov/targets/%s\" target=\"_blank\">%s</a>", gsymb, gsymb)][, .(Symbol = symb_htm, ensemblId, Name = geneName, Family = geneFamily, TDL, N_study, N_trait)][order(Symbol)]
    return(gg)
  })
  
  output$traits <- DT::renderDataTable({
    DT::datatable(data=trait_table_htm(), rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)

  output$genes <- DT::renderDataTable({
    DT::datatable(data=gene_table_htm(), rownames=F, options=list(autoWidth=T, dom='tipf'), escape=F)
  }, server=T)
  

  gHits_export <- reactive({
    if (is.null(gHits())) { return(NULL) }
    hits_out <- gHits()
    hits_out[["Trait"]] <- query_trait_name()
    hits_out[["TraitId"]] <- query_trait_id()
    hits_out <- hits_out[, .(Trait, TraitId, ensemblId, gsymb, geneName, geneFamily, TDL, n_study, study_N_mean, n_snp, n_snpw, geneNtrait, pvalue_mlog_median, or_median, rcras, muScore, muRank)]
    return(hits_out)
  })

  output$hits_file <- downloadHandler(
    filename = function() { sprintf("gwax_hits_%s.tsv", query_trait_id()) },
    content = function(file) {
      if (is.null(gHits_export())) { return(NULL) }
      write_delim(gHits_export(), file, delim="\t")
    }
  )
}
###
shinyApp(ui, server)

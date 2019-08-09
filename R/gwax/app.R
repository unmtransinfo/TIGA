##########################################################################################
### GWAX: GWAS Explorer
### gt = gene-trait data
### Dataset gt_stats.tsv from gwascat_gt_stats.R.
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
MIN_ASSN <- 20
#
t0 <- proc.time()
if (file.exists("gwax.Rdata")) {
  message(sprintf("Loading dataset from Rdata..."))
  load("gwax.Rdata")
  setDT(gt)
} else {
  message(sprintf("Loading dataset from files, writing Rdata..."))
#  gt <- read_delim("gt_stats.tsv.gz", '\t', col_types=cols(.default=col_character(), 
#    n_study=col_integer(), n_snp=col_integer(), n_traits_g=col_integer(), n_genes_t=col_integer(),
#    pvalue_mlog_median=col_double(), or_median=col_double(), rcras=col_double()))
  gt <- read_delim("gt_stats_mu.tsv.gz", '\t', col_types=cols(.default=col_character(), 
    n_study=col_integer(), n_snp=col_integer(), n_traits_g=col_integer(), n_genes_t=col_integer(),
    pvalue_mlog_median=col_double(), or_median=col_double(), study_N_median=col_integer(), rcras=col_double(),
    mu_score=col_double(), nAbove=col_integer(), nBelow=col_integer(), mu_rank=col_integer()))
  setDT(gt)
  # Why NAs in or_median?
  gt <- gt[!is.na(or_median)]
  ###
  traits_df <- gt[!is.na(or_median), .(.N),  by=c("trait", "trait_uri")]
  traits_df <- traits_df[N>=MIN_ASSN]
  traits <- sub("^.*/", "", traits_df$trait_uri) #named vector
  names(traits) <- traits_df$trait
  traits <- traits[order(names(traits))]
  save(gt, traits, file="gwax.Rdata")
}
#
t_elapsed <- (proc.time()-t0)[3]
#
message(sprintf("Gene count: %d", uniqueN(gt$gsymb)))
message(sprintf("Trait count (total): %d", uniqueN(gt$trait)))
message(sprintf("Trait count (n_assn>=%d): %d", MIN_ASSN, length(traits)))
message(sprintf("DEBUG: COUNT or_median: %d", sum(!is.na(gt$or_median))))
message(sprintf("DEBUG: COUNT pvalue_mlog_median: %d", sum(!is.na(gt$pvalue_mlog_median))))
message(sprintf("DEBUG: COUNT rcras: %d", sum(!is.na(gt$rcras))))
#
db_htm <- sprintf("<B>Dataset:</B> genes: %d ; traits: %d; top_traits: %d (t_load: %.1fs)",
	uniqueN(gt$gsymb), uniqueN(gt$trait), length(traits), t_elapsed)
###
tdls <- c("Tclin", "Tchem", "Tbio", "Tdark")
idgfams <- c("GPCR", "Kinase", "IC", "NR", "Other")
axes <- c("Effect", "Evidence")
#
qryTraitRand <- sample(traits, 1)
#
#############################################################################
HelpHtm <- function() {(
  sprintf("<P><B>GWAX</B>, GWAS Explorer, facilitates visualization and prioritization
of protein-coding genes associated with traits from genome-wide association studies
(GWAS), using the GWAS Catalog dataset from NIH-NHGRI and EBI.
<P>
<B>Dataset:</B>
<UL>
<LI>Well studied traits only available via this app, with minimum associations
(MIN_ASSN=%d, may be multiple for each gene). Traits are mapped to EFO, HPO,
Orphanet, PATO or GO, with 91%% mapped to EFO.
<LI>In this version, OR required and BETA ignored, due to problem of harmonizing varying BETA units.
</UL>
<B>Datatypes:</B>
<UL>
  <LI>GWAS Catalog:
  <UL>
  <LI><B>PVALUE_MLOG</B>: -LOG(p_value)
  <LI><B>OR</B>: odds ratio, inverted if &lt;1.
  <LI><B>INITIAL_SAMPLE_SIZE</B>: N_subjects
  <LI><B>REPORTED GENE(s)</B>: Gene(s) reported by author.
  <LI><B>MAPPED GENE(s)</B>: Gene(s) mapped to the strongest SNP.
  </UL>
  <LI>GWAX:
  <UL>
  <LI><B>N_gene</B>: total genes associated with trait.
  <LI><B>N_trait</B>: total traits associated with gene.
  <LI><B>N_snp</B>: SNPs involved with trait-gene association.
  <LI><B>N_study</B>: studies supporting trait-gene association.
  <LI><B>study_N</B>: median(INITIAL_SAMPLE_SIZE) supporting trait-gene association.
  <LI><B>RCRAS</B>: Relative Citation Ratio (RCR) Aggregated Score (iCite-RCR-based)
  </UL>
</UL>
Note that this app will accept query parameter <B>efo_id</B> via URL, e.g.
<B><TT>?efo_id=EFO_0000341</TT></B>.
Hits are filtered and ranked based on non-parametric multivariate &mu; scores
(Wittkowski, 2008).
<B>References:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\">GWAS Catalog</a>
<LI><a href=\"https://www.ebi.ac.uk/gwas/docs/fileheaders\">GWAS Catalog data dictionary</a>
<LI><a href=\"https://www.ebi.ac.uk/efo/\">Experimental Factor Ontology (EFO)</a>
</UL>
<B>Authors:</B> Jeremy Yang, Stephen Mathias, Cristian Bologa, Lars Juhl Jensen, Christophe Lambert, and Tudor Oprea.<BR/>
<B>Correspondence</B> from users of this app is welcome, and should be directed to 
<a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
  Data from <A HREF=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">The NHGRI-EBI GWAS Catalog</A>.<BR/>
  Built with R-Shiny &amp; Plotly.<BR/>
  This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>",
	MIN_ASSN)
)}

##########################################################################################
ui <- fluidPage(
  titlePanel(h2(sprintf("%s: GWAS Explorer", APPNAME), tags$em("(BETA)"), span(icon("lightbulb", lib="font-awesome"))),
      windowTitle=APPNAME),
  fluidRow(
    column(3, 
      wellPanel(
	selectInput("traitQry", "Query trait", choices=traits, selectize=T, selected=qryTraitRand),
        sliderInput("maxHits", "Max_hits", 50, 300, 100, step=50),
        sliderInput("minStudy", "Min_study", 1, 5, 1, step=1),
        checkboxGroupInput("tdl_filters", "TDL", choices=tdls, selected=tdls, inline=T),
        checkboxGroupInput("fam_filters", "Gene family", choices=idgfams, selected=idgfams, inline=T),
	checkboxGroupInput("logaxes", "LogAxes", choices=axes, selected=axes, inline=T),
        br(),
	actionButton("randQuery", "Demo", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
	actionButton("goRefresh", "Refresh", style='padding:4px;background-color:#DDDDDD;font-weight:bold'),
	actionButton("showHelp", "Help", style='padding:4px;background-color:#DDDDDD; font-weight:bold')
      )),
    column(9, plotlyOutput("plot", height = "600px"))),
  fluidRow(column(12, DT::dataTableOutput("datarows"))),
  fluidRow(column(12, downloadButton("hits_file", label="Download"))),
  fluidRow(column(12, wellPanel(htmlOutput(outputId="result_htm", height="60px")))),
  fluidRow(column(12, wellPanel(htmlOutput(outputId="log_htm", height="60px")))),
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
  bsTooltip("minStudy", "Minimum studies cutoff.", "top"),
  bsTooltip("maxHits", "Max hits cutoff, filtered by muScore.", "top"),
  bsTooltip("tdl_filters", "Filters: IDG Target Development Levels (TDLs).", "top"),
  bsTooltip("fam_filters", "Filters: IDG protein families.", "top"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

id2uri <- function(id) {
  if (grepl("^EFO", id)) {
    return(sprintf("http://www.ebi.ac.uk/efo/%s", id))
  } else if (grepl("^HP", id)) {
    return(sprintf("http://purl.obolibrary.org/obo/%s", id))
  } else if (grepl("^Orphanet", id)) {
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
      qryTraitRand <- sample(traits, 1)
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
    if (!(query_trait_uri() %in% gt$trait_uri)) { return(NULL) }
    return(gt$trait[gt$trait_uri==query_trait_uri()][1])
  })
  query_trait_id <- reactive({
    if (is.null(query_trait_uri())) { return(NULL) }
    return(sub("^.*/","", query_trait_uri()))
  })
  
  hits <- reactive({
    gt_this <- gt[(trait_uri==query_trait_uri()) & (n_study>=input$minStudy)]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[(fam %in% input$fam_filters) | (("Other" %in% input$fam_filters) & (is.na(fam) | !(fam %in% idgfams)))]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[tdl %in% intersect(input$tdl_filters, tdls)]
    if (nrow(gt_this)==0) { return(NULL) }   
    gt_this$tdl <- factor(gt_this$tdl, levels=c("Tclin", "Tchem", "Tbio", "Tdark", "NA"), ordered=T)
    gt_this <- gt_this[, .(gsymb, name, fam, tdl, n_study, n_snp, n_traits_g, pvalue_mlog_median, or_median, study_N_median, rcras, mu_score, nAbove, nBelow, mu_rank)]
    message(sprintf("DEBUG: hits() COUNT pvalue_mlog_median: %d", sum(!is.na(gt_this$pvalue_mlog_median))))
    gt_this[, ok := as.logical(mu_rank<=input$maxHits)]
    setorder(gt_this, mu_rank)
    return(gt_this)
  })

  output$result_htm <- reactive({
    htm <- sprintf("<B>Results:</B>")
    htm <- paste0(htm, sprintf("\"%s\"", query_trait_name()))
    htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", id2uri(query_trait_id()), query_trait_id()))
    message(sprintf("Query: \"%s\" (%s)", query_trait_name(), query_trait_id()))
    htm <- paste0(htm, "; minStudy = ", input$minStudy)
    if (!is.null(hits())) {
      htm <- paste0(htm, sprintf("; found N_genes: %d", nrow(hits())))
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

  output$plot <- renderPlotly({
    if (is.null(hits())) { return(NULL) }
    #xrange <- c(input$minStudy-.8, max(hits()$n_study, input$minStudy+2)+.5)
    xrange <- NULL
    xaxis <- list(title="Evidence (N_study)", range=xrange, type=ifelse("Evidence" %in% input$logaxes, "log", "normal"))
    yaxis <- list(title="Effect (OddsRatio)", type=ifelse("Effect" %in% input$logaxes, "log", "normal"))

    p <- plot_ly(type='scatter', mode='markers', data=hits()[(ok)],
        x=~n_study  + rnorm(nrow(hits()[(ok)]), sd=(.01*(max(hits()$n_study)-input$minStudy))), #Custom jitter
        y=~or_median,
        color=~tdl, colors=c("blue", "green", "red", "black", "gray"),
        marker=list(symbol="circle", size=pmax(10, 50/hits()[(ok)]$n_traits_g)),
        text=paste0(hits()[(ok)]$gsymb, ": ", hits()[(ok)]$name, "<br>",
        hits()[(ok)]$fam, ", ", hits()[(ok)]$tdl, "<br>",
        "N_study = ", hits()[(ok)]$n_study, "; ", "N_trait = ",
hits()[(ok)]$n_traits_g, "; N_snp = ", hits()[(ok)]$n_snp, "<br>",
        "OR = ", round(hits()[(ok)]$or_median, digits=2), "; ", 
        "study_N = ", hits()[(ok)]$study_N_median, "; ", 
        "pVal = ", sprintf("%.2g", 10^(-hits()[(ok)]$pvalue_mlog_median)), "<br>",
        "RCRAS = ", round(hits()[(ok)]$rcras, digits=2), "; ",
        "muScore = ", hits()[(ok)]$mu_score, "; ",
        "nAbove = ", hits()[(ok)]$nAbove, "; ",
        "nBelow = ", hits()[(ok)]$nBelow, "; ",
        "muRank = ", hits()[(ok)]$mu_rank
        )
      ) %>%
      layout(xaxis=xaxis, yaxis=yaxis, 
        title=paste0(query_trait_name(), "<br>", "(", input$traitQry, ")"),
        margin=list(t=100,r=50,b=60,l=60), legend=list(x=.9, y=.95, marker=list(size=20)), showlegend=T,
        font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N_gene = ", nrow(hits()), "; ", nrow(hits()[(ok)]), " shown)"), showarrow=F, x=.5, y=1, xref="paper", yref="paper")
    return(p)
  })
  
  #"gsymb","name","fam","tdl","n_study","n_snp","n_traits_g","pvalue_mlog_median","or_median","study_N_median", "rcras",
  #"mu_score", "nAbove", "nBelow", "mu_rank"
  output$datarows <- renderDataTable({
    if (is.null(hits())) { return(NULL) }
    DT::datatable(data=hits(), rownames=F,
	selection=list(target="row", mode="multiple", selected=NULL),
	class="cell-border stripe", style="bootstrap",
	options=list(autoWidth=T, dom='tip', #dom=[lftipr]
		columnDefs = list(list(className='dt-center', targets=c(0, 2:(ncol(hits())-2))), list(visible=F, targets=ncol(hits())-1)) #Here numbered from 0.
	),
	colnames=c("GSYMB","GeneName","idgFam","idgTDL","N_study","N_snp","N_trait","pVal_mlog","OR","study_N","RCRAS",
	           "muScore", "nAbove", "nBelow", "muRank", "ok")
  ) %>% formatRound(digits=2, columns=c(8,9,11)) #Here numbered from 1.
  }, server=T)

  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    hits_out <- hits()
    hits_out[["Trait"]] <- query_trait_name()
    hits_out[["TraitURI"]] <- query_trait_uri()
    hits_out <- hits_out[,c(ncol(hits_out)-1,ncol(hits_out),1:(ncol(hits_out)-2))]
    names(hits_out) <- c("Trait","TraitURI","GeneSymbol","GeneName","idgFam","idgTDL","N_study","N_snp","N_trait","pVal_mlog","OR","study_N","RCRAS","muScore","nAbove","nBelow","muRank")
    return(hits_out)
  })

  output$hits_file <- downloadHandler(
    filename = function() { "gwax_hits.tsv" },
    content = function(file) {
      if (is.null(hits_export())) { return(NULL) }
      write_delim(hits_export(), file, delim="\t")
    }
  )
}
###
shinyApp(ui, server)

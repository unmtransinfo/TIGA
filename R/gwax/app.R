##########################################################################################
### GWAX: GWAS Explorer
### gt = gene-trait data
### Dataset gt_stats.csv from gwascat_gt_stats.R.
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
} else {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  gt <- read_delim("gt_stats.tsv.gz", '\t')
  ###
  traits_df <- gt[!is.na(gt$or_median),] %>% group_by(trait, trait_uri) %>% summarise(count=n())
  traits_df <- traits_df[traits_df$count>=MIN_ASSN,]
  traits <- sub("^.*/", "", traits_df$trait_uri) #named vector
  names(traits) <- traits_df$trait
  traits <- traits[order(names(traits))]
  save(gt, traits, file="gwax.Rdata")
}
setDT(gt)
#
t_elapsed <- (proc.time()-t0)[3]
#
message(sprintf("Gene count: %d", uniqueN(gt$gsymb)))
message(sprintf("Trait count (total): %d", uniqueN(gt$trait)))
message(sprintf("Trait count (n_assn>=%d): %d", MIN_ASSN, length(traits)))
#
db_htm <- sprintf("<B>Dataset:</B> genes: %d ; traits: %d; top_traits: %d (t_load: %.1fs)",
	uniqueN(gt$gsymb), uniqueN(gt$trait), length(traits), t_elapsed)
###
tdls <- c("Tclin", "Tchem", "Tbio", "Tdark")
idgfams <- c("GPCR", "Kinase", "IC", "Enzyme", "NR", "Transporter", "TF", "Epigenetic", "TF; Epigenetic", "Other")
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
</UL>
<B>Datatypes:</B>
<UL>
  <LI>GWAS Catalog:
  <UL>
  <LI><B>PVALUE_MLOG</B>: -LOG(p_value)
  <LI><B>OR</B>: odds ratio, inverted if &lt;1.
  <LI><B>REPORTED GENE(s)</B>: Gene(s) reported by author.
  <LI><B>MAPPED GENE(s)</B>: Gene(s) mapped to the strongest SNP.
  </UL>
  <LI>GWAX:
  <UL>
  <LI><B>N_gene</B>: total genes associated with trait.
  <LI><B>N_trait</B>: total traits associated with gene.
  <LI><B>N_snp</B>: SNPs involved with trait-gene association.
  <LI><B>N_study</B>: studies supporting trait-gene association.
  </UL>
</UL>
Note that this app will accept query parameter <B>efo_id</B> via URL, e.g.
<B><TT>?efo_id=EFO_0000341</TT></B>.
<B>References:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\">GWAS Catalog</a>
<LI><a href=\"https://www.ebi.ac.uk/gwas/docs/fileheaders\">GWAS Catalog data dictionary</a>
<LI><a href=\"https://www.ebi.ac.uk/efo/\">Experimental Factor Ontology (EFO)</a>
</UL>
<B>Authors:</B> Jeremy Yang, Stephen Mathias, Cristian Bologa, and Tudor Oprea.<BR/>
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
        #sliderInput("minEffect", "Min_effect", 0, 100, 0, step = 10),
        #sliderInput("minSpec", "Min_specificity", 0, 1, .05, step = .05),
        checkboxGroupInput("tdl_filters", "TDL", choices=tdls, selected=tdls, inline=T),
        checkboxGroupInput("fam_filters", "Gene family", choices=idgfams, selected=idgfams, inline=T),
        br(),
	actionButton("randQuery", "Demo", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
	actionButton("goRefresh", "Refresh", style='padding:4px;background-color:#DDDDDD;font-weight:bold'),
	actionButton("showHelp", "Help", style='padding:4px;background-color:#DDDDDD; font-weight:bold')
      )),
    column(9, plotlyOutput("plot", height = "500px"))),
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
  #bsTooltip("minEffect", "Minimum effect cutoff.", "top"),
  #bsTooltip("minSpec", "Minimum specificity cutoff.", "top"),
  bsTooltip("tdl_filters", "Filters: IDG Target Development Levels (TDLs).", "top"),
  bsTooltip("fam_filters", "Filters: IDG protein families.", "top"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("efo_logo", "Experimental Factor Ontology (EFO)", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

tdl2color <- function(tdl) {
  colors <- rep("#CCCCCC", length(tdl))
  colors[tdl=="Tdark"] <- "#111111"
  colors[tdl=="Tbio"] <- "red"
  colors[tdl=="Tchem"] <- "green"
  colors[tdl=="Tclin"] <- "blue"
  return(colors)
}

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
    
  trait_uri <- reactive({
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
  trait_name <- reactive({
    if (!(trait_uri() %in% gt$trait_uri)) { return(NULL) }
    return(gt$trait[gt$trait_uri==trait_uri()][1])
  })
  trait_id <- reactive({
    if (is.null(trait_uri())) { return(NULL) }
    return(sub("^.*/","", trait_uri()))
  })
  
  hits <- reactive({
    #gt_this <- gt[(trait_uri==trait_uri()) & (or_median>input$minEffect) & (n_traits_g<=(1/input$minSpec))]
    gt_this <- gt[(trait_uri==trait_uri())]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[(fam %in% intersect(input$fam_filters, idgfams)) | (("Other" %in% input$fam_filters) & is.na(fam))]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[tdl %in% intersect(input$tdl_filters, tdls)]
    if (nrow(gt_this)==0) { return(NULL) }   

    gt_this$pvalue_mlog_median <- round(gt_this$pvalue_mlog_median, digits=2)
    gt_this$or_median <- round(gt_this$or_median, digits=2)
    gt_this[["tdl_color"]] <- tdl2color(gt_this$tdl)
    gt_this <- gt_this[,c("gsymb","name","fam","tdl","tdl_color","n_study","n_snp","n_traits_g","pvalue_mlog_median","or_median")]
    gt_this <- gt_this[order(-gt_this$or_median,gt_this$n_traits_g),]
    return(gt_this)
  })

  output$result_htm <- reactive({
    htm <- sprintf("<B>Results:</B>")
    htm <- paste0(htm, sprintf("\"%s\"", trait_name()))
    htm <- paste0(htm, sprintf(" (<a target=\"_blank\" href=\"%s\">%s</a>)", id2uri(trait_id()), trait_id()))
    message(sprintf("Query: \"%s\" (%s)", trait_name(), trait_id()))
    #htm <- paste0(htm, "; minEffect = ", input$minEffect)
    #htm <- paste0(htm, "; minSpecificity = ", input$minSpec)
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
    xaxis <- list(title="Specificity (1/N_trait)")
    yaxis <- list(title="Effect (OddsRatio)")

    p <- plot_ly() %>%
      add_trace(x=(1/hits()$n_traits_g), y=hits()$or_median,  
        type='scatter', mode='markers',
        marker=list(symbol="circle", color=hits()$tdl_color, size=20*log(hits()$n_study+1)),
        text=paste0(hits()$gsymb, ": ", hits()$name, "<br>",
        hits()$fam, ", ", hits()$tdl, "<br>",
        "N_trait = ", hits()$n_traits_g, "; N_snp = ",hits()$n_snp, "; N_study = ", hits()$n_study, "<br>",
        "pVal_MLOG = ", hits()$pvalue_mlog_median, "; ", "OR = ", hits()$or_median)
      ) %>%
      layout(xaxis=xaxis, yaxis=yaxis, 
        title=paste0(trait_name(), "<br>", "(", input$traitQry, ")"),
        margin=list(t=100,r=50,b=60,l=60),
        font=list(family="monospace", size=16)
      ) %>%
      add_annotations(text=paste0("(N_gene = ", nrow(hits()), ")"), showarrow=F, x=.5, y=1, xref="paper", yref="paper")
    return(p)
  })
  
  #"gsymb","name","fam","tdl","tdl_color","n_study","n_snp","n_traits_g","pvalue_mlog_median","or_median"
  
  output$datarows <- renderDataTable({
    if (is.null(hits())) { return(NULL) }
    DT::datatable(data=hits(), rownames=F,
	selection=list(target="row", mode="multiple", selected=NULL),
	class="cell-border stripe", style="bootstrap",
	options=list(dom='tip', #dom=[lftipr]
		autoWidth=T,
		columnDefs = list(
		list(className='dt-center', targets=c(0,2:(ncol(hits())-1))),
		list(visible=F, targets=c(4))
		)
	),
	colnames=c("GSYMB","GeneName","idgFam","idgTDL","tdl_color","N_study","N_snp","N_trait","pVal_mlog","OR")
  ) %>% formatRound(digits=2, columns=c(8,9))
  }, server=T)

  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    hits_out <- hits()
    hits_out$tdl_color <- NULL
    hits_out[["Trait"]] <- trait_name()
    hits_out[["TraitURI"]] <- trait_uri()
    hits_out <- hits_out[,c(ncol(hits_out)-1,ncol(hits_out),1:(ncol(hits_out)-2))]
    names(hits_out) <- c("Trait","TraitURI","GeneSymbol","GeneName","idgFam","idgTDL","N_study","N_snp","N_trait","pVal_mlog","OR")
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

##########################################################################################
### GWAX: GWAS Explorer
### gt = gene-trait data
### Dataset gt_stats.csv from gwascat_gt_stats.R.
### Jeremy Yang
##########################################################################################
library(readr)
library(shiny, quietly = T)
library(shinyBS, quietly=T)
library(DT, quietly=T)
library(dplyr, quietly = T)
library(plotly, quietly = T)

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
t_elapsed <- (proc.time()-t0)[3]
#
db_htm <- sprintf("<B>Dataset:</B> genes: %d ; traits: %d; top_traits: %d (t_load: %.1fs)",
	length(unique(gt$gsymb)), length(unique(gt$trait)), length(traits), t_elapsed)
###
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
(MIN_ASSN=%d).
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
<B>References:</B>
<UL>
<LI><a href=\"https://www.ebi.ac.uk/gwas/\">GWAS Catalog</a>
<LI><a href=\"https://www.ebi.ac.uk/gwas/docs/fileheaders\">GWAS Catalog data dictionary</a>
</UL>
<B>Authors:</B> Jeremy Yang, Stephen Mathias, Cristian Bologa, and Tudor Oprea.<BR/>
<B>Correspondence</B> from users of this app is welcome, and should be directed to 
<a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
  Data from <A HREF=\"https://www.ebi.ac.uk/gwas/\" TARGET=\"_blank\">GTEx, The NHGRI-EBI GWAS Catalog</A>.<BR/>
  Built with R-Shiny &amp; Plotly.<BR/>
  This work was supported by the National Institutes of Health grant U24-CA224370.<BR/>",
	MIN_ASSN)
)}

##########################################################################################
ui <- fluidPage(
  titlePanel(h2(sprintf("%s: GWAS Explorer", APPNAME), em("(BETA)")), windowTitle=APPNAME),
  fluidRow(
    column(3, 
      wellPanel(
	selectInput("trait_qry", "Query trait", choices=traits, selectize=T, selected=sample(traits, 1)),
        sliderInput("eff_min", "Min_effect", 0, 100, 0, step = 10),
        sliderInput("spc_min", "Min_specificity", 0, 1, .05, step = .05),
        checkboxGroupInput("options", "Options", choices=c("idgList"="idglist", "idgPfam"="pfam","Tclin","Tchem","Tbio","Tdark"),
		selected=c("Tclin","Tchem","Tbio","Tdark"), inline=T),
        br(),
        actionButton("randQuery", "Demo", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
        actionButton("goRefresh", "Refresh", style='padding:4px; background-color:#DDDDDD;font-weight:bold'),
        actionButton("showHelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
      )),
    column(9, plotlyOutput("plot", height = "500px"))),
  hr(),
  fluidRow(column(12, DT::dataTableOutput("datarows"))),
  fluidRow(column(12, downloadButton("hits_file", label="Download"))),
  fluidRow(column(12, wellPanel(htmlOutput(outputId="result_htm", height="60px")))),
  fluidRow(column(12, wellPanel(htmlOutput(outputId="log_htm", height="60px")))),
  fluidRow(
    column(12, em(strong(sprintf("%s", APPNAME)), " web app from ", 
        tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
        " and ",
        tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
        " data from ",
        tags$a(href="https://www.ebi.ac.uk/gwas/", target="_blank", span("GWAS Catalog", tags$img(id="gwas_catalog_logo", height="50", valign="bottom", src="GWAS_Catalog_logo.png")))
        ))),
  bsTooltip("trait_qry", "Select from most studied traits.", "right"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("gwas_catalog_logo", "GWAS Catalog, The NHGRI-EBI Catalog of published genome-wide association studies", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

##########################################################################################
server <- function(input, output, session) {
  observeEvent(input$showHelp, {
    showModal(modalDialog(title=HTML(sprintf("<H2>%s Help</H2>", APPNAME)),
      HTML(HelpHtm()),
      easyClose=T, footer=tagList(modalButton("Dismiss"))))
  })
  
  trait_uri <- reactive({
    input$goRefresh # Re-run this and downstream on action button.
    if (is.null(input$trait_qry)) { return(NULL) }
    paste0("http://www.ebi.ac.uk/efo/", input$trait_qry)
  })
  trait <- reactive({
    if (!(trait_uri() %in% gt$trait_uri)) { return(NULL) }
    gt$trait[gt$trait_uri==trait_uri()][1]
  })
  
  hits <- reactive({
    message(sprintf("DEBUG: nrow(gt) = %d", nrow(gt)))
    message(sprintf("DEBUG: %s: %s", trait_uri(), trait()))
    
    gt_this <- gt[gt$trait_uri==trait_uri(),]
    #gt_this <- gt_this[!is.na(gt_this$gsymb),]
    if (nrow(gt_this)==0) { return(NULL) }

    gt_this$trait <- NULL
    gt_this$trait_uri <- NULL
    gt_this$n_genes_t <- NULL
    gt_this <- gt_this[order(-gt_this$or_median,gt_this$n_traits_g),]

    #gt_this <- gt_this[gt_this$pvalue_mlog_median>input$sig_min,]
    gt_this <- gt_this[gt_this$or_median>input$eff_min,]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[gt_this$n_traits_g<=(1/input$spc_min),] 
    if (nrow(gt_this)==0) { return(NULL) }
    if ("pfam" %in% input$options) gt_this <- gt_this[!is.na(gt_this$fam),]
    if (nrow(gt_this)==0) { return(NULL) }
    if ("idglist" %in% input$options) gt_this <- gt_this[gt_this$idg2=="TRUE",]
    if (nrow(gt_this)==0) { return(NULL) }
    gt_this <- gt_this[gt_this$tdl %in% input$options,]
    message("DEBUG: nrow(gt_this) = ", nrow(gt_this))
    if (nrow(gt_this)==0) { return(NULL) }   

    gt_this$pvalue_mlog_median <- round(gt_this$pvalue_mlog_median, digits=1)
    gt_this$or_median <- round(gt_this$or_median, digits=2)

    gt_this[["tdl_color"]] <- NA
    gt_this$tdl_color[gt_this$tdl=="Tdark"] <- "gray"
    gt_this$tdl_color[gt_this$tdl=="Tbio"] <- "red"
    gt_this$tdl_color[gt_this$tdl=="Tchem"] <- "green"
    gt_this$tdl_color[gt_this$tdl=="Tclin"] <- "blue"
    gt_this$tdl_color[is.na(gt_this$tdl_color)] <- "#CCCCCC" #unmapped:
   
    gt_this <- gt_this[,c("gsymb","name","fam","tdl","tdl_color","idg2","n_study","n_snp","n_traits_g","pvalue_mlog_median","or_median")]
    return(gt_this)
  })

  output$result_htm <- reactive({
    htm <- sprintf("<B>Results:</B>")
    htm <- paste0(htm, sprintf("Trait: \"%s\" (%s)", trait(), trait_uri()))
    htm <- paste0(htm, "; eff_min = ", input$eff_min)
    htm <- paste0(htm, "; spc_min = ", input$spc_min)
    if (!is.null(hits())) {
      htm <- paste0(htm, sprintf("; found: %d", nrow(hits())))
    } else {
      htm <- paste0(htm, sprintf("; No gene-trait stats found."))
    }
    return(htm)
  })

  output$log_htm <- reactive({
    htm <- db_htm
    return(htm)
  })
  
  output$plot <- renderPlotly({
    if (is.null(hits())) { return(NULL) }

    #xaxis = list(type=ifelse("logx" %in% input$viewoptions, "log", "linear"), title="Specificity (1/n_trait)")
    #yaxis = list(type=ifelse("logy" %in% input$viewoptions, "log", "linear"), title="Effect (median(OR))")
    xaxis = list(title="Specificity (1/n_trait)")
    yaxis = list(title="Effect (median(OR))")

    p <- plot_ly() %>%
      add_trace(x = 1 / hits()$n_traits_g, y = hits()$or_median,  
        type = 'scatter', mode = 'markers',
        marker = list(symbol = "circle", color = hits()$tdl_color, size = 20*log(hits()$n_study+1)),
        text = paste0(hits()$gsymb, ": ", hits()$name, "<br>",
        hits()$fam, ", ", hits()$tdl, "<br>",
        "traits = ", hits()$n_traits_g, " ; SNPs = ",hits()$n_snp, " ; studies = ", hits()$n_study, "<br>",
        "pVal_MLOG = ", hits()$pvalue_mlog_median, "<br>",
        "OR = ", hits()$or_median)
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
        title = paste0(trait(), "<br>", "(", input$trait_qry, ")"),
        margin = list(t=100,r=50,b=60,l=60),
        font = list(family = "monospace", size = 16)
      ) %>%
      add_annotations(text=paste0("(N_gene = ", nrow(hits()), ")"), showarrow=F, x=0.0, y=1.0, xref="paper", yref="paper")
    return(p)
  })
  
  #"gsymb","name","fam","tdl","tdl_color","idg2","n_study","n_snp","n_traits_g","pvalue_mlog_median","or_median"
  
  output$datarows <- renderDataTable({
    if (is.null(hits())) { return(NULL) }
    DT::datatable(data=hits(), rownames=F,
	selection=list(target="row", mode="multiple", selected=NULL),
	class="cell-border stripe", style="bootstrap",
	options=list(
		autoWidth=T,
		columnDefs = list(
		list(className='dt-center', targets=c(0,2:(ncol(hits())-1))),
		list(visible=F, targets=c(4))
		)
	),
	colnames=c("GSYMB","GeneName","idgFam","idgTDL","tdl_color","idgList","N_study","N_snp","N_traits","pVal_mlog","OR")
  ) %>% formatRound(digits=2, columns=c(9,10))
  }, server=T)

  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    hits_out <- hits()
    hits_out$tdl_color <- NULL
    hits_out[["Trait"]] <- trait()
    hits_out[["TraitURI"]] <- trait_uri()
    hits_out <- hits_out[,c(ncol(hits_out)-1,ncol(hits_out),1:(ncol(hits_out)-2))]
    names(hits_out) <- c("Trait","TraitURI","GeneSymbol","GeneName","idgFam","idgTDL","idgList","N_study","N_snp","N_traits","pVal_mlog","OR")
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

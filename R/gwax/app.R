##########################################################################################
### GWAX: GWAS Explorer
### 
### gt = gene-trait data
### 
### options("shiny.port" = 3839)
###
### Jeremy Yang
##########################################################################################
library(shiny, quietly = T)
library(dplyr, quietly = T)
library(plotly, quietly = T)
#library(readr)


#traits <- c(
# "Alzheimer disease" = "EFO_0000249",
# "asthma" = "EFO_0000270",
# "bipolar disorder" = "EFO_0000289",
# "body height" = "EFO_0004339",
# "body mass index" = "EFO_0004340",
# "breast carcinoma" = "EFO_0000305",
# "colorectal cancer" = "EFO_0005842",
# "coronary heart disease" = "EFO_0001645",
# "Crohns disease" = "EFO_0000384",
# "high density lipoprotein cholesterol" = "EFO_0004612",
# "low density lipoprotein cholesterol" = "EFO_0004611",
# "multiple sclerosis" = "EFO_0003885",
# "Parkinsons disease" = "EFO_0002508",
# "prostate carcinoma" = "EFO_0001663",
# "rheumatoid arthritis" = "EFO_0000685",
# "schizophrenia" = "EFO_0000692",
# "triglyceride measurement" = "EFO_0004530",
# "type II diabetes mellitus" = "EFO_0001360",
# "unipolar depression" = "EFO_0003761"
#)

##########################################################################################
#get_gt_stats <- function() {
#  dbcon <- dbConnect(MySQL(), host="localhost", dbname="gwascatalog")
#  gt_stats <- dbGetQuery(dbcon,"SELECT * FROM gt_stats")
#  dbDisconnect(dbcon)
#  return(gt_stats)
#}

##########################################################################################
#On launch, not at search time.  I.e. not reactive.
#gt_stats <- get_t_stats()

#gt <- read_csv("gt_stats.csv.gz")
fgzin <- gzfile("gt_stats.csv.gz","rt")  #"rt" means read-text 
gt <- read.csv(fgzin, header=T)
close(fgzin)

##########################################################################################
### Select traits with most evidence.  Top 100.  Simple metric: rows of gt.
traits <- function() {
  
  trait2uri <- unique(gt[!is.na(gt$or_median),c("trait","trait_uri")])
  traits_df <- gt[!is.na(gt$or_median),] %>% group_by(trait_uri) %>% summarise(count = n())
  traits_df <- merge(traits_df, trait2uri, by="trait_uri")
  traits_df <- traits_df[order(-traits_df$count),]
  traits_vec <- sub("^.*/", "", traits_df$trait_uri) #named vector
  names(traits_vec) <- traits_df$trait
  traits_vec <- traits_vec[1:100]
  traits_vec[order(names(traits_vec))]
}

##########################################################################################
ui <- fluidPage(
  
  titlePanel(HTML("GWAX: GWAS Explorer <I>(DEMO PROTOTYPE)</I>")),
  
  fluidRow(
    column(3, 
      wellPanel(
	      selectInput("trait_qry", "Query trait", choices = traits(), selectize=T, selected = sample(traits(), 1)),
	      HTML("<I><SMALL>(Select from top 100 most studied traits.)</SMALL></I>")),
  
      wellPanel(
        icon("cog", "fa-2x"),
        sliderInput("eff_min", "Min_effect", 0, 100, 0, step = 10),
        sliderInput("spc_min", "Min_specificity", 0, 1, .05, step = .05),
        checkboxGroupInput("tdls", "TDLS", choices=c("Tclin","Tchem","Tbio","Tdark"), selected=c("Tclin","Tchem","Tbio","Tdark"), inline=T),
        checkboxGroupInput("filters","Filters",choices=c("IDGlist"="idglist", "Pfam"="pfam"), inline=T)
      ),
      wellPanel(
        checkboxGroupInput("viewoptions", "Viewoptions", choices=c("LogX"="logx", "LogY"="logy"), selected = c("logx","logy"), inline=T)
      )
    ),
    column(9, plotlyOutput("plot", height = "600px"))),
  hr(),
  fluidRow(column(12, dataTableOutput("hits"))),
  fluidRow(column(12, downloadButton("hits_file", label="Download"))),
  fluidRow(column(12, HTML("<I>Web app built with R-shiny and plotly, at UNM for the NIH IDG project.</I>")))
)

##########################################################################################
server <- function(input, output, session) {

  trait_uri <- reactive({
    paste0("http://www.ebi.ac.uk/efo/", input$trait_qry)
  })
  trait <- reactive({
    gt$trait[gt$trait_uri == trait_uri()][1]
  })
  
  hits <- reactive({
    
    if (is.na(gt) || nrow(gt)==0) # HANDLE ERROR HOW??   
    {
      output$status <- renderText("No gene-trait stats found.")
      message("DEBUG: No gene-trait stats found.")
      return(NULL)
    }
    message(sprintf("DEBUG: nrow(gt) = %d", nrow(gt)))
    
    message(sprintf("DEBUG: %s: %s", trait_uri(), trait()))
    
    gt_this <- gt[gt$trait_uri==trait_uri(),]
    gt_this <- gt_this[!is.na(gt_this$gsymb),]
    gt_this <- gt_this[!is.na(gt_this$or_median),]
    gt_this <- gt_this[!is.na(gt_this$name),] #Removes many unmapped.  Fix?

    gt_this$trait <- NULL
    gt_this$trait_uri <- NULL
    gt_this$n_genes_t <- NULL
    gt_this <- gt_this[order(-gt_this$or_median,gt_this$n_traits_g),]

    message("DEBUG: eff_min = ", input$eff_min)
    message("DEBUG: spc_min = ", input$spc_min)
    if (nrow(gt_this)>0) { 
      #gt_this <- gt_this[gt_this$pvalue_mlog_median>input$sig_min,]
      gt_this <- gt_this[gt_this$or_median>input$eff_min,]
    }
    if (nrow(gt_this)>0) { 
      gt_this <- gt_this[gt_this$n_traits_g<=1/input$spc_min,] 
    }
   
    message("DEBUG: nrow(gt_this) = ", nrow(gt_this))
    #if (nrow(gt)==0) # HANDLE ERROR HOW??   

    gt_this$pvalue_mlog_median <- round(gt_this$pvalue_mlog_median, digits=1)
    gt_this$or_median <- round(gt_this$or_median, digits=2)

    gt_this$tdl_color <- NA
    gt_this$tdl_color[gt_this$tdl=="Tdark"] <- "gray"
    gt_this$tdl_color[gt_this$tdl=="Tbio"] <- "red"
    gt_this$tdl_color[gt_this$tdl=="Tchem"] <- "green"
    gt_this$tdl_color[gt_this$tdl=="Tclin"] <- "blue"
    gt_this$tdl_color[is.na(gt_this$tdl_color)] <- "#CCCCCC" #unmapped:
   
    if ("pfam" %in% input$filters) gt_this <- gt_this[!is.na(gt_this$fam),]
    if ("idglist" %in% input$filters) gt_this <- gt_this[gt_this$idg2=="TRUE",]

    gt_this <- gt_this[gt_this$tdl %in% input$tdls,] 
    gt_this
  })
  
  output$plot <- renderPlotly({
     
    xaxis = list(type=ifelse("logx" %in% input$viewoptions, "log", "linear"), title="Specificity (1/n_trait)")
    yaxis = list(type=ifelse("logy" %in% input$viewoptions, "log", "linear"), title="Effect (median(OR))")
     
    p <- plot_ly() %>%
      add_trace(x = 1 / hits()$n_traits_g, y = hits()$or_median,  
        type = 'scatter', mode = 'markers',
        marker = list(symbol = "circle", color = hits()$tdl_color, size = 20*log(hits()$n_study+1)),
        text = paste0(hits()$gsymb, ": ", hits()$name, "<br>",
        hits()$fam, ", ", hits()$tdl, "<br>",
        "traits = ", hits()$n_traits_g, " ; snps = ",hits()$n_snp, " ; studies = ", hits()$n_study, "<br>",
        "p_mlog_median = ", hits()$pvalue_mlog_median, "<br>",
        "or_median = ", hits()$or_median)
      ) %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
        title = paste0(trait(), "<br>", "(", input$trait_qry, ")"),
        margin = list(t=100,r=50,b=60,l=60),
        font = list(family = "monospace", size = 16)
      ) %>%
      add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=1.0, xref="paper", yref="paper") %>%
      add_annotations(text=paste0("(N_gene = ", nrow(hits()), ")"), showarrow=F, x=0.0, y=1.0, xref="paper", yref="paper")
    return(p)
  })
  
  output$hits <- renderDataTable(hits(), options = list(pageLength=20, pagingType="simple", lengthChange=F, searching=F, info=T))
  
  output$hits_file <- downloadHandler(
    #filename not working...
    filename = function() {
      paste0("gwax_hits-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(hits(), file, row.names=F)
    }
  )
}
###

shinyApp(ui, server)

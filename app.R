#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ActivePathways)
library(visNetwork)
library(shinyWidgets)
library(shinyscreenshot)
library(shinyalert)
library(tools)

# to-do: visEvent function
source(file = "./enrichmentMap/enrichmentMap.R")
options(shiny.maxRequestSize = 100*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Shiny ActivePathways"),
    helpText("Upload score file and gmt file to analyze enrichment"),
    useShinyalert(),
    sidebarLayout(
        sidebarPanel(
            tags$p(
                "ShinyActivePathways can run ActivePathways method for 
                enrichment analysis and visualize result as enrichment map"
            ),
            tags$p(
                "Upload scoresFile and gmtFile, then click 'Run ActivePathways!'"
            ),
            fileInput(inputId = "scoresFile", 
                      label = "Score file",
                      accept = c(
                          ".tsv",
                          ".csv"
                      )
            ),
            fileInput(inputId = "gmtFile", 
                      label = "GMT file",
                      accept = c(
                          ".gmt"
                      )
            ),
            actionButton(inputId = "runPanel", 
                         label = "Run ActivePathways!"),
            downloadLink('downloadNetwork', 'Download network as .html'),
            # sliderInput(inputId = "size", "Label Size: ", min = 0, max = 30, value = 14),
            shinyscreenshot::screenshotButton(id = "visNet", filename = "downloadNetwork", label = "download"),
            shinyWidgets::sliderTextInput(
                inputId = "size",
                label = "Font Size: ", 
                choices = c(0, 7, 14, 21, 28),
                grid = TRUE,
                selected = 14
            )
        ),
        

        # Show a plot
        mainPanel(
            mainPanel(
                visNetworkOutput(outputId = "visNet", 
                                 width = "130%", 
                                 height = "700px"),
                hr(),
                verbatimTextOutput('geneList')
                
               
            )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    # Uploads score files
    observeEvent(input$runPanel, {
        shinyalert::shinyalert(html = TRUE, text = tagList(
            textInput(inputId = "cutoff", value = 0.1, label = "Cutoff"),
            textInput(inputId = "significant", value = 0.05, label = "Significant"),
            selectInput(inputId = "mergeMethod", choices = c("Brown", "Fisher"), selected = "Brown", label = "Merge Method"), 
            selectInput(inputId = "correctionMethod", choices = c("holm", 
                                                                  "fdr", 
                                                                  "hochberg", 
                                                                  "hommel", 
                                                                  "bonferroni", 
                                                                  "BH", 
                                                                  "BY",
                                                                  "none"), 
                        selected = "holm",
                        label = "correction method"), 
            sliderInput(inputId = "filterRange", 
                        "Geneset Filter:",
                        min = 0, max = 1000,
                        value = c(5, 1000)),
            type = "input"
        ))
    })
    
    scores <- reactive({
        print(input$scoresFile)
        data <- NULL
        if (file_ext(input$scoresFile$datapath) == "tsv") {
            data <- as.matrix(
                read.table(file = input$scoresFile$datapath, header = TRUE, sep = '\t', row.names = 'Gene'))
        } else if (file_ext(input$scoresFile$datapath) == "csv") {
            data <- as.matrix(
                read.table(file = input$scoresFile$datapath, header = TRUE, sep = ',', row.names = 'Gene'))
        }
        
        data[is.na(data)] <- 1
        return(data)
    })
    
    # Uploads gmt file
    gmt <- reactive({
        return(read.GMT(input$gmtFile$datapath))
    })
    enrichmentResult <- reactiveValues(data = NULL)
    enrichNetwork <- reactiveValues(data = NULL)
    
    # Render visNetwork
    output$visNet <- renderVisNetwork({
        if (! is.null(enrichNetwork$data)) {
            enrichNetwork$data %>%
                visInteraction(navigationButtons = TRUE) %>%
                visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes); 
            ;}") # for display genes contained in each gene set
        }
    })
    
    # For display genes contained in each gene set
    selectedPathwayId <- reactiveValues(data = NULL)
    observeEvent(input$current_node_id, {
        print(paste("The selected node ID is:", input$current_node_id))
        gmtData <- getDataFromGMT(gmt())
        selectedPathwayId$data <- gmtData$ids[gmtData$pathwayNames == input$current_node_id]
        print(paste("The pathway ID is:", selectedPathwayId$data))
    })
    
    # Render information of gene list in text
    output$geneList <- renderPrint({
        if (! is.null(selectedPathwayId$data)) {
            gmt()[[selectedPathwayId$data]]
        }})
    
    # Modify network properties
    observe({
        visNetworkProxy("visNet") %>% 
            visNodes(font = list(size = input$size)) # for change font size of network
        #  Download network as html files
        output$downloadNetwork <- downloadHandler(
            filename = function() {
                paste('network-', Sys.Date(), '.html', sep='')
            },
            content = function(con) {
                enrichNetwork$data %>% visSave(con)
            }
        )
    })
    
    
    
    observeEvent(input$shinyalert,{
        enrichmentResult$data <- ActivePathways(score = scores(), 
                                                gmt = gmt(), 
                                                cutoff = as.numeric(input$cutoff), 
                                                significant = as.numeric(input$significant), 
                                                merge.method = input$mergeMethod, 
                                                correction.method = input$correctionMethod,
                                                geneset.filter = input$filterRange)
        
        g <- plotEnrichmentMap(gmt(), 
                               enrichmentResult$data,
                               algorithm = "jaccard", 
                               similarityCutoff = 0.25, 
                               pvalueCutoff = NULL, 
                               k = 0.5)
        visigraph <- visNetwork::visIgraph(g)
        nodes <- visigraph$x$nodes
        edges <- visigraph$x$edges
        enrichNetwork$data <- visNetwork(nodes, edges) %>% 
            visIgraphLayout() %>%
            visExport(type = "png", name = "export-network", 
                      float = "left", label = "Save network", style= "") 
    })
    
    observeEvent(input$download, {
        shinyscreenshot::screenshot(download = TRUE)
        print(input$shinyscreenshot)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

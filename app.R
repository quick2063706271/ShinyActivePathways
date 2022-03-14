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
library(dplyr)

# to-do: visEvent function
source(file = "./enrichmentMap/enrichmentMap.R")
source(file = "./color/color.R")
options(shiny.maxRequestSize = 100*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
    navbarPage(title = "ShinyActivePathways",
               collapsible = TRUE,
               tabPanel("Visualization",
    # Application title
    # titlePanel("Shiny ActivePathways"),
    # helpText("Upload score file and gmt file to analyze enrichment"),
    # useShinyalert(),
    
                    sidebarLayout(
                        sidebarPanel(
                            tabsetPanel(
                                        tabPanel("Input Files", wellPanel(
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
                                            actionButton(inputId = "runEx", 
                                                         label = "Run Example"),
                                            
                                            # sliderInput(inputId = "size", "Label Size: ", min = 0, max = 30, value = 14),
                                            
                                            )),
                                        tabPanel("Adjust", wellPanel(
                                            shinyWidgets::sliderTextInput(
                                                inputId = "size",
                                                label = "Font Size: ", 
                                                seq(from = 0, to = 28),
                                                selected = 14
                                            ),
                                            sliderInput(
                                                inputId = "edgeCutoff",
                                                label = "Edge cutoff (Similarity): ", 
                                                min = 0.25,
                                                max = 1.0,
                                                value = 0.25
                                            ),
                                            numericInput("edgePrecise",
                                                         label = "Value:",
                                                         value = 0.25,
                                                         min = 0.25,
                                                         max = 1.0),
                                            checkboxInput("highlightNeighbors",
                                                          label = "highlight neighbors"),
                                            checkboxInput("navigationButtons",
                                                          label = "hide navigation buttons")
                                            )
                                        ),
                                        tabPanel("Download", wellPanel(
                                            downloadLink('downloadNetwork', 'Download network as .html'),
                                            )
                                        )
                                    
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
                ),
                tabPanel("Tutorial",)
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
            hr(),
            selectInput(inputId = "metric",
                        choices = c("jaccard", "overlap", "combined"),
                        selected = "jaccard",
                        label = "Overlap measures"),
            type = "input"
        ))
    })
    
    scores <- reactive({
        print(input$scoresFile)
        preTable <- NULL
        if (file_ext(input$scoresFile$datapath) == "tsv") {
            preTable <- read.table(file = input$scoresFile$datapath, header = TRUE, sep = '\t')
            
        } else if (file_ext(input$scoresFile$datapath) == "csv") {
            preTable <- read.table(file = input$scoresFile$datapath, header = TRUE, sep = ',')
        }
        row.names(preTable) <- preTable[,1]
        preTable <- preTable[, 2: length(colnames(preTable))]
        data <- as.matrix(preTable)
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
    # observe({
    #     enric
    # })
    observe({
        visNetworkProxy("visNet") %>% 
            visNodes(font = list(size = input$size)) %>% # for change font size of network
            visOptions(highlightNearest = input$highlightNeighbors) %>% # for highlight neighbors
            visInteraction(navigationButtons = !input$navigationButtons)
            
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
    
    
    
    
    fullResult <- reactiveValues(nodes = NULL, edges = NULL)
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
                               algorithm = input$metric, 
                               similarityCutoff = 0.25, 
                               pvalueCutoff = NULL, 
                               k = 0.5)
        visigraph <- visNetwork::visIgraph(g)
        nodes <- visigraph$x$nodes
        nodes$color <- getColors(enrichmentResult$data)
        edges <- visigraph$x$edges
        edges$id <- paste0("e", seq_along(edges$from))
        edges$color <- rgb(137,207,240, max = 255)
        fullResult$nodes <- nodes
        fullResult$edges <- edges
        enrichNetwork$data <- visNetwork(nodes, edges) %>% 
            visIgraphLayout() %>%
            visExport(type = "png", name = "export-network", 
                      float = "left", label = "Save network", style= "") 
    })
    
    

    observe({
        updateSliderInput(
            session = session,
            inputId = "edgeCutoff",
            value = input$edgePrecise
        )
    })
    
    
    observe({
        updateSliderInput(
            session = session,
            inputId = "edgePrecise",
            value = input$edgeCutoff
        )
    })
    
    observe({
        
        filteredResult <- fullResult$edges[fullResult$edges$weight < input$edgeCutoff, ]
        print(filteredResult)
        filteredEdges <- filteredResult$id
        hiddenEdges <- fullResult$edges[! fullResult$edges$id %in% filteredEdges,]
        if (!is.null(filteredResult) & !is.null(filteredEdges)) {
            visRemoveEdges(visNetworkProxy("visNet"), id = filteredEdges)
            visUpdateEdges(visNetworkProxy("visNet"), edges = hiddenEdges) # stop at Feb 6 2022
        }
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

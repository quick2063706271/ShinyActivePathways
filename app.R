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

# to-do: visEvent function
source(file = "./enrichmentMap/enrichmentMap.R")
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Shiny ActivePathways"),

    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "scoresFile", 
                      label = "Score file",
                      accept = c(
                          ".tsv"
                      )
            ),
            fileInput(inputId = "gmtFile", 
                      label = "GMT file",
                      accept = c(
                          ".gmt"
                      )
            ),
            actionButton(inputId = "Calculate", 
                         label = "Run ActivePathways!"),
            downloadLink('downloadNetwork', 'Download network as .html'),
            checkboxInput("showVertexLabel",
                          label = "show vertex label",
                          value = TRUE)
        ),
        

        # Show a plot
        mainPanel(
            mainPanel(
                visNetworkOutput(outputId = "visNet", 
                                 width = "130%", 
                                 height = "400px"),
                hr(),
                verbatimTextOutput('geneList')
                
               
            )
        )
    )
)

# Define server logic
server <- function(input, output) {

    scores <- reactive({
        print(input$scoresFile)
        data <- as.matrix(
            read.table(file = input$scoresFile$datapath, header = TRUE, sep = '\t', row.names = 'Gene'))
        data[is.na(data)] <- 1
        return(data)
    })
    gmt <- reactive({
        return(read.GMT(input$gmtFile$datapath))
    })
    enrichmentResult <- reactiveValues(data = NULL)
    enrichNetwork <- reactiveValues(data = NULL)
    output$visNet <- renderVisNetwork({
        if (! is.null(enrichNetwork$data)) {
            enrichNetwork$data %>%
                visInteraction(navigationButtons = TRUE) %>%
                visOptions(nodesIdSelection = TRUE) %>%
                visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_id', nodes.nodes);
            ;}")
        }
    })
    selectedPathwayId <- reactiveValues(data = NULL)
    observeEvent(input$current_node_id, {
        print(paste("The selected node ID is:", input$current_node_id))
        gmtData <- getDataFromGMT(gmt())
        selectedPathwayId$data <- gmtData$ids[gmtData$pathwayNames == input$current_node_id]
        print(paste("The pathway ID is:", selectedPathwayId$data))
    })
    output$geneList <- renderPrint({
        if (! is.null(selectedPathwayId$data)) {
            gmt()[[selectedPathwayId$data]]
        }})
    # if (! input$showVertexLabel) {
    #     visNetworkProxy("visNet") %>%
    #         visNodes(font = list(size = 0))
    # } else {
    #     visNetworkProxy("visNet") %>%
    #         visNodes(font = list(size = 14))
    # }
    
    
    observeEvent(input$Calculate,{
        enrichmentResult$data <- ActivePathways(scores(), gmt())
        
        g <- plotEnrichmentMap(gmt(), 
                               enrichmentResult$data,
                               algorithm = "jaccard", 
                               similarityCutoff = 0.25, 
                               pvalueCutoff = NULL, 
                               k = 0.5)
        visigraph <- visNetwork::visIgraph(g)
        nodes <- visigraph$x$nodes
        edges <- visigraph$x$edges
        enrichNetwork$data <- visNetwork(nodes, edges) %>% visIgraphLayout()
        output$downloadNetwork <- downloadHandler(
            filename = function() {
                paste('network-', Sys.Date(), '.html', sep='')
            },
            content = function(con) {
                visg %>% visSave(con)
            }
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

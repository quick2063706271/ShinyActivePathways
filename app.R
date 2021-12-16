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

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = "scoresFile", label = "Score file",
                      accept = c(
                          ".tsv"
                      )
            ),
            fileInput(inputId = "gmtFile", label = "GMT file",
                      accept = c(
                          ".gmt"
                      )
            ),
            actionButton(inputId = "Calculate", label = "Run ActivePathways!"),
            downloadLink('downloadNetwork', 'Download network as .html')
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
            mainPanel(
                visNetworkOutput(outputId = "visNet", 
                                 width = "130%", 
                                 height = "700px")
                
            )
        )
    )
)

# Define server logic required to draw a histogram
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
    observeEvent(input$Calculate,{
        enrichmentResult$data <- ActivePathways(scores(), gmt())
        
        g <- plotEnrichmentMap(gmt(), 
                               enrichmentResult$data,
                               algorithm = "jaccard", 
                               similarityCutoff = 0.25, 
                               pvalueCutoff = NULL, 
                               k = 0.5)
        visg <- visNetwork::visIgraph(g)
        output$visNet <- renderVisNetwork({
            visg
        })
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

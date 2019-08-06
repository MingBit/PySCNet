library(shiny)
library(shinythemes)
library(DT)
library(reticulate)
library(shinyjs)

ui = tagList(
  # shinythemes::themeSelector(),
  
  navbarPage(
    theme = shinytheme("superhero"),
    "PySCNet Visualization",
    
    tabPanel(h5("Expression Data"),
     sidebarPanel(
      id = 'sidebar',
      fileInput("file1", "Choose gNetData File:",
                multiple = FALSE,
                accept = ".pk"),
      textInput("txt", "Cell Filter:"),
      textInput("txt", "Gene Filter:"),
      textInput("txt", "Search:"),
      actionButton("action1", "Search", class = "btn-primary"),
      actionButton("clean_button", "Clean", class = 'btn-primary')
    ),
    mainPanel(
      titlePanel('Expression table'),
      DT::dataTableOutput('exp')
    )
    ),
    
    
    tabPanel(h5("Cell Attributes"), 
             fluidPage(
               h3('Cell Attributes'),
               # DT::dataTableOutput('cell_info'),
               hr(),
               h3('Summary'),
               plotOutput('cell_summary')
             )), 
    
    tabPanel(h5("Gene Attributes"), 
            fluidPage(
      h3('Gene Attributes'),
      # DT::dataTableOutput('gene_info'),
      hr(),
      h3('Summary'),
      plotOutput('gene_summary')
    )),
    
    tabPanel(h5("Network Attributes"),   fluidPage(
      h3('Linkage Table'),
      DT::dataTableOutput('links'),
      hr(),
      h3('Network')
      
    ))
  )
)


shinyApp(ui = ui, server = server)
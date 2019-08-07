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
      titlePanel('Expression table')
      # DT::dataTableOutput('exp')
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
    
    tabPanel(h5("Network Attributes"),
             sidebarPanel(
               h3("Create your own network"),
               hr(),
               radioButtons('trans_genes', 'Gene or Transcription Factor', c('Gene', 'TF', 'TF/Gene')),
               selectInput('cluster_id', 'Select Cluster ID', choices = 'clusterid'),
               selectInput('gene_module', 'Select Gene Module', choices = 'genemodule'),
               sliderInput('weight_filter', "Weight cut by", 0.2, min = 0, max = 1, step = .05),
               hr(),
               sliderInput("node_opacity", "Node Opacity", 0.8, min = 0.1,
                           max = 1, step = .1),
               sliderInput("label_size", "Label Fontsize", 25, min = 10,
                           max = 50, step = 5),
               selectInput('node_size', "Node Size encoded by", 
                           choices = c('pageRank', 'degree', 'betweenness', 'closeness'))
             ),
             mainPanel(
             
               h3('Network'),
               simpleNetworkOutput('network'),
              fluidPage(
               hr(),
               h3('Linkage Table'),
               DT::dataTableOutput('links', width = "80%"),
               hr(),
               h3('Gene Module Table'),
               DT::dataTableOutput('gene_module')))
     )
  )
)

shinyApp(ui = ui, server = server)
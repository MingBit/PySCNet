#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggpubr)
library(dplyr)
library(igraph)
library(reticulate)
library(shinythemes)
library(DT)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(tagList(
  # shinythemes::themeSelector(),
  
  navbarPage(
    theme = shinytheme("superhero"),
    "PySCNet Visualization",
    
    # tabPanel(h5("Expression Data"),
    #  sidebarPanel(
    #   id = 'sidebar',
    #   fileInput("file1", "Choose gNetData File:",
    #             multiple = FALSE,
    #             accept = ".pk")
    # 
    # ),
    # mainPanel(
    #   titlePanel('Expression table')
    #   # DT::dataTableOutput('exp')
    # ),
    
    
    tabPanel(h5("Cell Attributes"), 
             sidebarPanel(
               id = 'sidebar',
               fileInput("file1", "Choose gNetData File:",
                         multiple = FALSE,
                         accept = ".pk")
               
             ),
             mainPanel(
               fluidPage(
                 h3('Cell Attributes'),
                 DT::dataTableOutput('cell_info'),
                 hr(),
                 h3('Summary'),
                 plotOutput('cell_summary')
               ))), 
    
    tabPanel(h5("Gene Attributes"), 
             fluidPage(
               h3('Gene Attributes'),
               DT::dataTableOutput('gene_info'),
               hr(),
               h3('Summary'),
               plotOutput('gene_summary')
             )),
    
    tabPanel(h5("Network Attributes"),
             sidebarPanel(
               h3("Create your own network"),
               hr(),
               radioButtons('trans_genes', 'Gene or Transcription Factor', c('Gene', 'TF', 'TF/Gene')),
               uiOutput('clusterid'),
               # uiOutput('genemodule'),
               # sliderInput('weight_filter', "Weight cut by", 0, min = 0, max = 1, step = .05),
               selectInput('top_edges', 'Select Top Edges', choices = c('50', '100', '200', '500'), selected = '200'),
               hr(),
               sliderInput("node_opacity", "Node Opacity", 1, min = 0.1,
                           max = 1, step = .1),
               sliderInput("label_size", "Label Fontsize", 15, min = 10,
                           max = 50, step = 5),
               selectInput('node_size', "Node Size encoded by", 
                           choices = c('degree', 'betweenness', 'closeness', 'pageRank'))
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
))


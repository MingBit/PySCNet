library(ggpubr)
library(dplyr)
library(igraph)

server = function(input, output) {
  source_python('~/MING_V9T/PhD_Pro/PySCNet/PySCNet/Preprocessing/gnetdata.py')
  
  pk_upload <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    pickle_data = load_Gnetdata_object(input$file1$datapath)
    return(pickle_data)
  })
  
  get_links <- reactive({
    pickle_data <- pk_upload()
    links <- data.frame(cell_clusterid = pickle_data$NetAttrs$links['cell_clusterid'],
                        source = pickle_data$NetAttrs$links['source'],
                         target = pickle_data$NetAttrs$links['target'],
                         weight = pickle_data$NetAttrs$links['weight'])
    return(links)
  })
  
  output$exp <- DT::renderDataTable({
    
    pickle_data <- pk_upload()
    datatable(pickle_data$ExpMatrix)%>% 
      formatStyle(names(pickle_data$ExpMatrix), backgroundColor = "grey") %>% 
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')})
  
  output$cell_info <- DT::renderDataTable({
    pickle_data <- pk_upload()
    datatable(pickle_data$CellAttrs) %>% formatStyle(names(pickle_data$CellAttrs), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')})
  
  output$gene_info <- DT::renderDataTable({
    pickle_data <- pk_upload()
    datatable(pickle_data$GeneAttrs) %>% formatStyle(names(pickle_data$GeneAttrs), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')})
  
  output$cell_summary <- renderPlot({
    pickle_data <- pk_upload()
    cell_cluster_df = plyr::count(pickle_data$CellAttrs, vars = c('cluster_id', 'cell_type'))
    cell_count <- column_to_rownames(as.data.frame(aggregate(freq ~ cell_type, data = cell_cluster_df, sum)),'cell_type')
    cell_cluster_df$percent <- paste(format((cell_cluster_df$freq / cell_count[cell_cluster_df$cell_type,]) * 100, digits = 2),"%", sep = "")
    
    ggbarplot(cell_cluster_df, x = 'cluster_id', y = "freq",
                    fill = "cell_type",
                    # color = c('#FB4944', '#9659D9'),
                    palette =  c('#033F63', '#BF4342'),
                    xlab = "Clusters",
                    ylab = "Cell Counts",
                    label = "percent", lab.vjust = 1.2, 
                    lab.size = 6) + theme(text = element_text(size = 20))
  })
  
  output$gene_summary <- renderPlot({
    pickle_data <- pk_upload()
    hist(pickle_data$NetAttrs$links['weight'])
  })
  
  output$clusterid <- renderUI({
    
    links <- get_links()
    selectInput('cluster_id', 'Select Cluster ID', choices = as.list(unique(links$cell_clusterid)), selected = '1')
    
  })

  output$genemodule <- renderUI({
    
    links <- get_links()
    selectInput('gene_module', 'Select Gene Module', choices = as.list(unique(links$cell_clusterid)), selected = '1')
    
  })
  
  output$links <- DT::renderDataTable({
    
    links <- get_links() 
    datatable(links, options = list(lengthMenu = c(5, 30, 50), pageLength = 10), filter = 'top', style = 'jqueryui') %>%
      formatStyle(c('cell_clusterid', 'source', 'target', 'weight'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
  })
  
  output$gene_module <- DT::renderDataTable({
    
    links <- get_links()
    datatable(links, options = list(lengthMenu = c(5, 30, 50), pageLength = 10), style = 'jqueryui') %>%
      formatStyle(c('cell_clusterid','source', 'target', 'weight'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
  })
  
  output$network <- renderForceNetwork({
    
    pickle_data <- pk_upload()
    gene_info <- pickle_data$GeneAttrs %>% mutate(id = seq(0, nrow(.)-1))
    
    links <- get_links() %>% dplyr::filter(cell_clusterid == input$cluster_id) %>% 
      mutate(source_id = gene_info[match(source, gene_info$Gene),'id']) %>%
      mutate(target_id = gene_info[match(target, gene_info$Gene),'id'])
    
    links <- links[order(links$weight, decreasing = T), ][1:as.numeric(input$top_edges),]
    
    
    forceNetwork(Links = links,
                 Source = "source_id", Nodes = gene_info, opacityNoHover = 1,linkWidth = 5,
                 Target = "target_id", Value = 1, NodeID = "Gene", Group = 'TF_Gene',
                 height = 500, width = 1000, fontSize = input$label_size, zoom = TRUE,
                 opacity = input$node_opacity, colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"))
  })
  

  
  observeEvent(input$clean_button, {
    shinyjs::reset("sidebar")
  })
}
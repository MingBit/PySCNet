
library(dplyr)
library(ggpubr)
library(igraph)
library(networkD3)
library(reticulate)
library(tibble)
# Sys.setenv("RSTUDIO_PANDOC"="/usr/lib/rstudio/bin/pandoc")
server = function(input, output) {
  use_python('/home/mwu/miniconda3/bin/python')
  import_from_path('PySCNet', '/home/mwu/MING_V9T/PhD_Pro/PySCNet/')
  source_python('/home/mwu/MING_V9T/PhD_Pro/PySCNet/PySCNet/Preprocessing/gnetdata.py')
  
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
  
  create_network <- reactive({
    
    links <- get_links() %>% dplyr::filter(cell_clusterid == input$cluster_id)
    links <- links[order(links$weight, decreasing = T), ][1:as.numeric(input$top_edges),]
    
    graph <- graph_from_edgelist(as.matrix(links[,c('source', 'target')]), directed = FALSE)
    gene_communities <- cluster_louvain(graph)
    
    gene_module <- data.frame(Gene = as.character(names(membership(gene_communities))),
                              module = as.factor(gene_communities$membership)) %>%
      mutate(id = seq(0, nrow(.)-1)) %>%
      mutate(degree = degree(graph)) %>%
      mutate(betweenness = betweenness(graph)) %>%
      mutate(closeness = closeness(graph)) %>%
      mutate(pageRank = page_rank(graph)$vector)
    
    
    links <- links %>% mutate(source_id = gene_module[match(source, gene_module$Gene),'id']) %>%
      mutate(target_id = gene_module[match(target, gene_module$Gene),'id'])
    
    
    return(list(gene_module, links))
    
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
  
  # output$genemodule <- renderUI({
  #   
  #   links <- get_links()
  #   selectInput('gene_module', 'Select Gene Module', choices = as.list(unique(links$cell_clusterid)), selected = 'All')
  #   
  # })
  
  output$links <- DT::renderDataTable({
    
    links <- get_links() 
    datatable(links, options = list(lengthMenu = c(5, 30, 50), pageLength = 10), filter = 'top', style = 'jqueryui') %>%
      formatStyle(c('cell_clusterid', 'source', 'target', 'weight'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
  })
  
  
  output$gene_module <- DT::renderDataTable({
    
    datatable(as.data.frame(create_network()[1]),
              options = list(lengthMenu = c(5, 30, 50), pageLength = 10), style = 'jqueryui') %>%
      formatStyle(c('module', 'id', 'degree', 'betweenness', 'closeness'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
    
  })
  
  output$network <- renderForceNetwork({
    tmp = create_network()
    # if(input$gene_module == 'All') {Nodes = as.data.frame(tmp[1])} 
    # else{Nodes = as.data.frame(tmp[1]) %>% dplyr::filter(module == input$gene_module)}
    # Links = as.data.frame(tmp[2]) %>% dplyr::filter(source %in% Nodes['Gene'] & target %in% Nodes['Gene'])
    net <- forceNetwork(Links = as.data.frame(tmp[2]),
                        Source = "source_id", Nodes = as.data.frame(tmp[1]), opacityNoHover = 1,linkWidth = 5,
                        Target = "target_id", Value = 'weight', NodeID = "Gene", Group = 'module',
                        height = 500, width = 1000, fontSize = input$label_size, zoom = TRUE, Nodesize = input$node_size,
                        opacity = input$node_opacity, colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"))
    
  })
  
  
  observeEvent(input$clean_button, {
    shinyjs::reset("sidebar")
  })
}


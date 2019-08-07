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
    links <- data.frame(source = pickle_data$NetAttrs$links['source'],
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
  
  output$links <- DT::renderDataTable({
    # pickle_data <- pk_upload()
    # datatable(data.frame(source = pickle_data$NetAttrs$links['source'],
    #                      target = pickle_data$NetAttrs$links['target'],
    #                      weight = pickle_data$NetAttrs$links['weight'])) %>% 
    links <- get_links() 
    datatable(links, options = list(lengthMenu = c(5, 30, 50), pageLength = 5)) %>%
      formatStyle(c('source', 'target', 'weight'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
  })
  
  output$gene_module <- DT::renderDataTable({
    # pickle_data <- pk_upload()
    links <- get_links()
    # datatable(data.frame(source = pickle_data$NetAttrs$links['source'],
    #                      target = pickle_data$NetAttrs$links['target'],
    #                      weight = pickle_data$NetAttrs$links['weight'])) %>% 
    datatable(links, options = list(lengthMenu = c(5, 30, 50), pageLength = 5)) %>%
      formatStyle(c('source', 'target', 'weight'), backgroundColor = "grey") %>%
      formatStyle(0, target = 'row', backgroundColor = 'black', fontWeight = 'bold')
  })
  
  output$cell_summary <- renderPlot({
    pickle_data <- pk_upload()
    hist(pickle_data$NetAttrs$links['weight'])
  })
  
  output$gene_summary <- renderPlot({
    pickle_data <- pk_upload()
    hist(pickle_data$NetAttrs$links['weight'])
  })
  
  output$network <- renderForceNetwork({
    # links <- get_links()
    # gene_info <- pk_upload()$GeneAttrs
    data("MisLinks")
    data("MisNodes")
    forceNetwork(Links = MisLinks, Source = "source", Nodes = MisNodes,
                 Target = "target", Value = "value", NodeID = "name",
                 height = 500, width = 1000, fontSize = input$label_size, zoom = TRUE,
                 Group = "group", opacity = input$node_opacity,  colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"))
  })
  
  observeEvent(input$clean_button, {
    shinyjs::reset("sidebar")
  })
}
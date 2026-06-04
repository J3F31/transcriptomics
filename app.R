# app.R
library(shiny)
library(DT)  
library(EnhancedVolcano)
library(bslib)
library(reactable)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)


deg6 <- readRDS("~/Desktop/Shiny Test/Test/deg_6.rds")
deg36 <- readRDS("~/Desktop/Shiny Test/Test/deg_36.rds")


PlotEnhancedVolcano <- function(deg, yCutoff, pCutOff,Fcutoff, genenames) {
  
  EnhancedVolcano(deg, x = 'log2FoldChange', y = 'pvalue',
                  lab = deg$geneName,
                  selectLab = genenames,
                  subtitle = bquote(italic('')),
                  labSize = 4, pointSize = 1.5, axisLabSize=10, titleLabSize=12,
                  subtitleLabSize=8, captionLabSize=10,
                  legendPosition = "none",
                  drawConnectors = TRUE,
                  ylim = c(0, yCutoff),
                  pCutoff = pCutOff, FCcutoff = Fcutoff)
}

FilterDEGs <- function(deg ,Pvalue , Fvalue) {
  
  S <- filter(deg , pvalue <= Pvalue )
  
  S <- filter(S , abs(log2FoldChange) > Fvalue )
  
  return(S)
}

GoEnrichmentAnalysis <- function(group_genes, pvalue , Ont) {
genes_6 <- bitr(row.names(FilterDEGs(group_genes,0.05,0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
  
  go_results <- enrichGO(gene          = genes_6$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENTREZID', 
                         ont           = Ont,
                         pAdjustMethod = "none",
                         pvalueCutoff  = pvalue,
                         qvalueCutoff  = 1,
                         readable      = TRUE)
  
  #View(as.data.frame(go_results))
  return(go_results)
  
}

KeggEnrichmentAnalysis <- function(group_genes, pvalue) {
  genes_6 <- bitr(row.names(FilterDEGs(group_genes,0.05,0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
  
  kegg_results <- enrichKEGG(gene = genes_6$ENTREZID,
                             organism = "hsa",
                             keyType = "kegg",
                             pvalueCutoff = pvalue,
                             pAdjustMethod = "none",
                             minGSSize = 10,
                             maxGSSize = 500,
                             qvalueCutoff = 1,
                             use_internal_data = FALSE)
  
  #View(as.data.frame(kegg_results))
  return(kegg_results)
}

ReactomeEnrichmentAnalysis <- function(group_genes, pvalue) {
  genes_6 <- bitr(row.names(FilterDEGs(group_genes,0.05,0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
  
  ora_results <- enrichPathway(gene = genes_6$ENTREZID, 
                               organism = "human",
                               pvalueCutoff = pvalue,
                               pAdjustMethod = "none",
                               qvalueCutoff = 1)
  
  #View(as.data.frame(ora_results))
  return(ora_results)
}

ui <- fluidPage(
 
  navset_pill(
    
    nav_panel("DEG 6 Hours",
              
      selectInput( 
        "select", 
        "Select options below:", 
        list("Volcano Plot" = "one","ENRICHMENT ANALYSIS" = "1B", "KEGG ANALYSIS" = "1C","Reactome ANALYSIS" = "1D"
          ) 
        ),
      conditionalPanel(
        condition = "input.select == 'one'",
        fluidRow(  
          column(5, numericInput( 
            "PID6", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
          column(5, numericInput( 
            "FID6", 
            "Log2 Fold Change", 
            value = 1, 
            min = 0, 
            max = 1 
          ))),
        plotOutput("volcanoOut6"),
        dataTableOutput("table6"),
        ),
      
      conditionalPanel(
        condition = "input.select == '1B'",
        fluidRow(  
          column(5, numericInput( 
            "PVE", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
          column(5,  selectInput( 
            "selectOnt", 
            "Select options below:", 
            list("biological process" = "BP","cellular component " = "CC", "molecular function" = "MF"
            ) 
          )
          )),
        plotOutput("Goenrich"),
        dataTableOutput("tableEnrich"),
      ),
      
      conditionalPanel(
        condition = "input.select == '1C'",
        fluidRow(  
          column(5, numericInput( 
            "PVK", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
         ),
        plotOutput("Kegg"),
        dataTableOutput("tableKegg"),
      ),
      
      conditionalPanel(
        condition = "input.select == '1D'",
        fluidRow(  
          column(5, numericInput( 
            "PVR", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
         ),
        plotOutput("Reactome"),
        dataTableOutput("tableReactome"),
      )
      
      
      ), 
    
    nav_panel("DEG 36 Hours",
              
      selectInput( 
        "select", 
        "Select options below:", 
        list("Volcano Plot" = "one","ENRICHMENT ANALYSIS" = "1B", "KEGG ANALYSIS" = "1C","Reactome ANALYSIS" = "1D"
        ) 
      ),
      
      
      conditionalPanel(
        condition = "input.select == 'one'",
        fluidRow(  
          column(5, numericInput(  
            "PID36", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )),
          column(5, numericInput(  
            "FID36", 
            "Log2 Fold Change", 
            value = 1, 
            min = 0, 
            max = 1 
          ))),
        plotOutput("volcanoOut36"),
        dataTableOutput("table36")
      )

    ), 
    
    nav_panel("6 Hours VS 36 Hours",
              
    imageOutput("image") 
  ), 
  ), id = "tab",
)


server <- function(input, output, session) {
  observe(print(input$select))
   
    output$volcanoOut6 <- renderPlot({
      
    req(deg6 )

    tryCatch({
      PlotEnhancedVolcano(deg6, 4 , input$PID6, input$FID6, FilterDEGs(deg6, input$PID6, input$FID6)[input$table6_rows_selected,]$geneName)
    }, error = function(e) {
      output$error_msg <- renderText(e$message)
      NULL
    })
  })
  
  output$volcanoOut36 <- renderPlot({
    req(deg36)
    
    tryCatch({
      PlotEnhancedVolcano(deg36, 10.5 , input$PID36, input$FID36, FilterDEGs(deg36, input$PID36, input$FID36)[input$table36_rows_selected,]$geneName)
    }, error = function(e) {
      output$error_msg <- renderText(e$message)
      NULL
    })
  })
  
  output$table6 <- renderDataTable({
    datatable(FilterDEGs(deg6, input$PID6, input$FID6))
  }) 
  
  output$table36 <- renderDataTable({
    datatable(FilterDEGs(deg36, input$PID36, input$FID36))
  }) 
  
  output$image <- renderImage(
    
    { 
      list(src = "6h_vs_36h_venn_diagram.png", height = "100%") 
    }, 
    deleteFile = FALSE 
  )
  
  output$Goenrich <- renderPlot({
    req(deg6)
    
    tryCatch({
      dotplot(GoEnrichmentAnalysis(deg6,input$PVE, input$selectOnt))
    }, error = function(e) {
      output$error_msg <- renderText(e$message)
      NULL
    })
  })
  
  output$tableEnrich <- renderDataTable({
    datatable(as.data.frame(GoEnrichmentAnalysis(deg6,input$PVE, input$selectOnt)))
  }) 
  
  output$Kegg <- renderPlot({
    req(deg6)
    
    tryCatch({
      dotplot(KeggEnrichmentAnalysis(deg6,input$PVK))
    }, error = function(e) {
      output$error_msg <- renderText(e$message)
      NULL
    })
  })
  
  output$tableKegg <- renderDataTable({
    datatable(as.data.frame(KeggEnrichmentAnalysis(deg6,input$PVK)))
  })
  
  output$Reactome <- renderPlot({
    req(deg6)
    
    tryCatch({
      dotplot(ReactomeEnrichmentAnalysis(deg6,input$PVR))
    }, error = function(e) {
      output$error_msg <- renderText(e$message)
      NULL
    })
  })
  
  output$tableReactome <- renderDataTable({
    datatable(as.data.frame(ReactomeEnrichmentAnalysis(deg6,input$PVR)))
  })
}

shinyApp(ui, server)

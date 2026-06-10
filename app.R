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


deg6 <- readRDS("./deg_6.rds")
deg36 <- readRDS("./deg_36.rds")


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

genes6 <- bitr(row.names(FilterDEGs(deg6,0.05,0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
genes36 <- bitr(row.names(FilterDEGs(deg36,0.05,0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

GoEnrichmentAnalysis <- function(group_genes, pvalue , Ont) {
  go_results <- enrichGO(gene          = group_genes$ENTREZID,
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
kegg_results <- enrichKEGG(gene = group_genes$ENTREZID,
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
ora_results <- enrichPathway(gene = group_genes$ENTREZID, 
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
        "select6", 
        "Select options below:", 
        list("Volcano Plot" = "volcano","ENRICHMENT ANALYSIS" = "go", "KEGG ANALYSIS" = "kegg","Reactome ANALYSIS" = "react"
          ) 
        ),
      conditionalPanel(
        condition = "input.select6 == 'volcano'",
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
        condition = "input.select6 == 'go'",
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
        condition = "input.select6 == 'kegg'",
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
        condition = "input.select6 == 'react'",
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
        "select36", 
        "Select options below:", 
        list("Volcano Plot" = "volcano","ENRICHMENT ANALYSIS" = "go", "KEGG ANALYSIS" = "kegg","Reactome ANALYSIS" = "react"
        ) 
      ),
      conditionalPanel(
        condition = "input.select36 == 'volcano'",
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
      ),
      conditionalPanel(
        condition = "input.select36 == 'go'",
        fluidRow(  
          column(5, numericInput( 
            "PVE36", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
          column(5,  selectInput( 
            "selectOnt36", 
            "Select options below:", 
            list("biological process" = "BP","cellular component " = "CC", "molecular function" = "MF"
            ) 
          )
          )),
        plotOutput("Goenrich36"),
        dataTableOutput("tableEnrich36"),
      ), 
      conditionalPanel(
        condition = "input.select36 == 'kegg'",
        fluidRow(  
          column(5, numericInput( 
            "PVK36", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        plotOutput("Kegg36"),
        dataTableOutput("tableKegg36"),
      ),
      
      conditionalPanel(
        condition = "input.select36 == 'react'",
        fluidRow(  
          column(5, numericInput( 
            "PVR36", 
            "P value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        plotOutput("Reactome36"),
        dataTableOutput("tableReactome36"),
      )

    ), 
    
    nav_panel("6 Hours VS 36 Hours",
              
    imageOutput("image") 
  ), 
  ), id = "tab",
)


server <- function(input, output, session) {
  
  #Outputs for 6 hours
  #Volcano 
  output$volcanoOut6 <- renderPlot({
    PlotEnhancedVolcano(deg6, 4 , input$PID6, input$FID6, FilterDEGs(deg6, input$PID6, input$FID6)[input$table6_rows_selected,]$geneName)
  })
  output$table6 <- renderDataTable({
    datatable(FilterDEGs(deg6, input$PID6, input$FID6))
  }) 
  #Go
  output$Goenrich <- renderPlot({
    dotplot(GoEnrichmentAnalysis(genes6,input$PVE, input$selectOnt))
  })
  output$tableEnrich <- renderDataTable({
    datatable(as.data.frame(GoEnrichmentAnalysis(genes6,input$PVE, input$selectOnt)))
  }) 
  #Kegg
  output$Kegg <- renderPlot({
    dotplot(KeggEnrichmentAnalysis(genes6,input$PVK))
  })
  output$tableKegg <- renderDataTable({
    datatable(as.data.frame(KeggEnrichmentAnalysis(genes6,input$PVK)))
  })
  #Reactome
  output$Reactome <- renderPlot({
    dotplot(ReactomeEnrichmentAnalysis(genes6,input$PVR))
  })
  output$tableReactome <- renderDataTable({
    datatable(as.data.frame(ReactomeEnrichmentAnalysis(genes6,input$PVR)))
  })
  
  #Outputs for 36 hours
  #Volcano
  output$volcanoOut36 <- renderPlot({
    PlotEnhancedVolcano(deg36, 10.5 , input$PID36, input$FID36, FilterDEGs(deg36, input$PID36, input$FID36)[input$table36_rows_selected,]$geneName)
  })
  output$table36 <- renderDataTable({
    datatable(FilterDEGs(deg36, input$PID36, input$FID36))
  })
  #Go
  output$Goenrich36 <- renderPlot({
    dotplot(GoEnrichmentAnalysis(genes36,input$PVE36, input$selectOnt36))
  })
  output$tableEnrich36 <- renderDataTable({
    datatable(as.data.frame(GoEnrichmentAnalysis(genes36,input$PVE36, input$selectOnt36)))
  }) 
  #Kegg
  output$Kegg36 <- renderPlot({
    dotplot(KeggEnrichmentAnalysis(genes36,input$PVK36))
  })
  output$tableKegg36 <- renderDataTable({
    datatable(as.data.frame(KeggEnrichmentAnalysis(genes36,input$PVK36)))
  })
  #Reactome
  output$Reactome36 <- renderPlot({
    dotplot(ReactomeEnrichmentAnalysis(genes36,input$PVR36))
  })
  output$tableReactome36 <- renderDataTable({
    datatable(as.data.frame(ReactomeEnrichmentAnalysis(genes36,input$PVR36)))
  })

  output$image <- renderImage(
    { 
      list(src = "6h_vs_36h_venn_diagram.png", height = "100%") 
    }, 
    deleteFile = FALSE 
  )
}

shinyApp(ui, server)

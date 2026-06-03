# app.R
library(shiny)
library(DT)  
library(EnhancedVolcano)
library(bslib)
library(reactable)
library(dplyr)
library(ggplot2)

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

ui <- fluidPage(
 
  navset_pill( 
    nav_panel("DEG 6 Hours",
              
    navset_pill( 
      nav_panel("ENRICHMENT ANALYSIS", "ONE"),
      nav_panel("KEGG ANALYSIS", "TWO"),
      nav_panel("Reactome ANALYSIS", "THREE")
                ),        
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
    
    nav_panel("DEG 36 Hours",
              
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
    
    nav_panel("6 Hours VS 36 Hours",
              
    imageOutput("image") 
  ), 
  ), 
  id = "tab", 
)


server <- function(input, output, session) {
    
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
      list(src = "vnn_plot.png", height = "100%") 
    }, 
    deleteFile = FALSE 
  )
}

shinyApp(ui, server)

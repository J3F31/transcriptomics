# app.R
library(shiny)
library(DT)  
library(EnhancedVolcano)
library(bslib)
library(reactable)
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(ggVennDiagram)
library(pathview)
library(png)
library(magick)

deg6 <- readRDS("./rds-files/deg_6.rds")

deg36 <- readRDS("./rds-files/deg_36.rds")

goBP6 <- readRDS("./rds-files/go_bp_6.rds")
goCC6 <- readRDS("./rds-files/go_cc_6.rds")
goMF6 <- readRDS("./rds-files/go_mf_6.rds")
go6 <- goBP6
kegg6 <- readRDS("./rds-files/kegg_6.rds")
react6 <- readRDS("./rds-files/react_6.rds")

goBP36 <- readRDS("./rds-files/go_bp_36.rds")
goCC36 <- readRDS("./rds-files/go_cc_36.rds")
goMF36 <- readRDS("./rds-files/go_mf_36.rds")
go36 <- goBP36
kegg36 <- readRDS("./rds-files/kegg_36.rds")
react36 <- readRDS("./rds-files/react_36.rds")

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

genes_6 <- bitr(rownames(FilterDEGs(deg6, 0.05, 0)), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

FilterEnrich <- function(data, pval) {
  return(filter(data, pvalue <= pval))
}

PlotEnrichmentAnalysis <- function(data, pval) {
  #Filter top 10 by padjust value
  t <- filter(data, pvalue <= pval)
  t <- t[order(t$p.adjust, decreasing = F),]
  t <- t[1:10,]
  t <- drop_na(t)
  
  #Parse gene ratio & order entries by count
  t$GeneRatio <- parse_ratio(t$GeneRatio)
  t <- t[order(t$Count),]
  t <- cbind(t, order = 1:nrow(t))
  
  #Plot
  ggplot(t, aes(x = GeneRatio, y = order, size = Count, color = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#C21807", high = "#6495ED") +
    labs(y = "") +
    scale_y_continuous(breaks = 1:nrow(t), labels = t$Description) +
    theme(text = element_text(size = 20))
}

PlotVennDiagram <- function(pval, fval) {
  d6 <- FilterDEGs(deg6, pval, fval)
  d36 <- FilterDEGs(deg36, pval, fval)
  x <- list(`+6h` = rownames(d6), `+36h` = rownames(d36))
  ggVennDiagram(x) +
    scale_fill_gradient(low = "#C21807", high = "#6495ED") +
    theme(legend.position = "none")
}

ui <- fluidPage(
 
  navset_pill(
    
    nav_panel("Sample +6h",
              
      selectInput( 
        "select6", 
        "Select Visualization", 
        list("Volcano Plot" = "volcano","Go Terms" = "go", "Kegg Analysis" = "kegg", "Kegg Pathway" = "path", "Reactome Analysis" = "react")
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
            value = 0, 
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
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
          column(5,  selectInput( 
            "selectOnt", 
            "Select Ontology", 
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
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
         ),
        plotOutput("Kegg"),
        dataTableOutput("tableKegg"),
      ),
      
      conditionalPanel(
        condition = "input.select6 == 'path'",
        fluidRow(  
          column(5, numericInput( 
            "PVpath", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        h3("Select entry to generate pathway"),
        plotOutput("pathview", width = 1626, height = "auto"),
        dataTableOutput("pathviewTable"),
      ),
      
      conditionalPanel(
        condition = "input.select6 == 'react'",
        fluidRow(  
          column(5, numericInput( 
            "PVR", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
         ),
        plotOutput("Reactome"),
        dataTableOutput("tableReactome"),
      )
      
      
      ), 
    
    nav_panel("Sample +36h",
              
      selectInput( 
        "select36", 
        "Select Visualization", 
        list("Volcano Plot" = "volcano","Go Terms" = "go", "Kegg Analysis" = "kegg", "Kegg Pathway" = "path","Reactome Analysis" = "react") 
      ),
      conditionalPanel(
        condition = "input.select36 == 'volcano'",
        fluidRow(  
          column(5, numericInput(  
            "PID36", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )),
          column(5, numericInput(  
            "FID36", 
            "Log2 Fold Change", 
            value = 0, 
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
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
          column(5,  selectInput( 
            "selectOnt36", 
            "Select Ontology", 
            list("Biological Process" = "BP","Cellular Component " = "CC", "Molecular Function" = "MF"
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
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        plotOutput("Kegg36"),
        dataTableOutput("tableKegg36"),
      ),
      
      conditionalPanel(
        condition = "input.select36 == 'path'",
        fluidRow(  
          column(5, numericInput( 
            "PVpath36", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        h3("Select entry to generate pathway"),
        plotOutput("pathview36", width = 1626, height = "auto"),
        dataTableOutput("pathviewTable36"),
      ),
      
      conditionalPanel(
        condition = "input.select36 == 'react'",
        fluidRow(  
          column(5, numericInput( 
            "PVR36", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1 
          )), 
        ),
        plotOutput("Reactome36"),
        dataTableOutput("tableReactome36"),
      )

    ), 
    
    nav_panel("Samples Compared",
      fluidRow(  
        column(
          5, 
          numericInput(  
            "PVvenn", 
            "P Value", 
            value = 0.05, 
            min = 0, 
            max = 1
          )
        ),
        column(
          5, 
          numericInput(  
            "LFCvenn", 
            "Log2 Fold Change", 
            value = 0, 
            min = 0, 
            max = 1 
          )
        )
      ),
      plotOutput("vennDiagram"),
      h3("Overlapping DEGs"),
      dataTableOutput("vennDiagramTable")
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
    PlotEnrichmentAnalysis(go6, input$PVE)
  })
  output$tableEnrich <- renderDataTable({
    datatable(FilterEnrich(go6, input$PVE))
  }) 
  #Kegg
  output$Kegg <- renderPlot({
    PlotEnrichmentAnalysis(kegg6, input$PVK)
  })
  output$tableKegg <- renderDataTable({
    datatable(FilterEnrich(kegg6, input$PVK))
  })
  #Pathview
  my_img <- reactive({
    id <- rownames(FilterEnrich(kegg6, input$PVpath)[input$pathviewTable_rows_selected,])
    if (identical(id, character(0))) return()
    
    pathview(gene.data = genes_6$ENTREZID, pathway.id = id, species = "hsa")
    img <- image_read(paste0("./", id, ".pathview.png"))
    invisible(file.remove(paste0("./", id, ".png")))
    invisible(file.remove(paste0("./", id, ".pathview.png")))
    invisible(file.remove(paste0("./", id, ".xml")))
    w <- image_info(img)$width
    h <- image_info(img)$height
    list(
      raster = as.raster(img),
      w = w,
      h = h
    )
  })
  img_dim_f <- function(parm) {
    function() {
      p <- 0
      i <- my_img()
      if (isTruthy(i)) {
        if (isTruthy(i[[parm]])) {
          p <- i[[parm]]
        }
      }
      p
    }
  }
  observe({
    id <- rownames(FilterEnrich(kegg6, input$PVpath)[input$pathviewTable_rows_selected,])
    if (identical(id, character(0))) return()
    output$pathview <- renderPlot(
      expr = {
        i <- req(my_img())
        r <- req(i$raster)
        plot(r)
      }, 
      width = img_dim_f("w"),
      height = img_dim_f("h")
    )
  })
  
  
  output$pathviewTable <- renderDataTable({
    datatable(FilterEnrich(kegg6, input$PVpath), selection = 'single')
  })
  #Reactome
  output$Reactome <- renderPlot({
    PlotEnrichmentAnalysis(react6, input$PVR)
  })
  output$tableReactome <- renderDataTable({
    datatable(FilterEnrich(react6, input$PVR))
  })
  #Go Ontology
  observe({
    input$selectOnt
    switch(
      input$selectOnt,
      "BP" = {
        go6 <- goBP6
        output$Goenrich <- renderPlot({
          PlotEnrichmentAnalysis(go6, input$PVE)
        })
        output$tableEnrich <- renderDataTable({
          datatable(FilterEnrich(go6, input$PVE))
        }) 
      },
      "CC" = {
        go6 <- goCC6
        output$Goenrich <- renderPlot({
          PlotEnrichmentAnalysis(go6, input$PVE)
        })
        output$tableEnrich <- renderDataTable({
          datatable(FilterEnrich(go6, input$PVE))
        }) 
      },
      "MF" = {
        go6 <- goMF6
        output$Goenrich <- renderPlot({
          PlotEnrichmentAnalysis(go6, input$PVE)
        })
        output$tableEnrich <- renderDataTable({
          datatable(FilterEnrich(go6, input$PVE))
        }) 
      }
    )
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
    PlotEnrichmentAnalysis(go36, input$PVE36)
  })
  output$tableEnrich36 <- renderDataTable({
    datatable(FilterEnrich(go36, input$PVE36))
  }) 
  #Kegg
  output$Kegg36 <- renderPlot({
    PlotEnrichmentAnalysis(kegg36, input$PVK36)
  })
  output$tableKegg36 <- renderDataTable({
    datatable(FilterEnrich(kegg36, input$PVK36))
  })
  #Pathview
  my_img36 <- reactive({
    id <- rownames(FilterEnrich(kegg36, input$PVpath36)[input$pathviewTable36_rows_selected,])
    if (identical(id, character(0))) return()
    
    pathview(gene.data = genes_36$ENTREZID, pathway.id = id, species = "hsa")
    img <- image_read(paste0("./", id, ".pathview.png"))
    invisible(file.remove(paste0("./", id, ".png")))
    invisible(file.remove(paste0("./", id, ".pathview.png")))
    invisible(file.remove(paste0("./", id, ".xml")))
    w <- image_info(img)$width
    h <- image_info(img)$height
    list(
      raster = as.raster(img),
      w = w,
      h = h
    )
  })
  img_dim_f36 <- function(parm) {
    function() {
      p <- 0
      i <- my_img36()
      if (isTruthy(i)) {
        if (isTruthy(i[[parm]])) {
          p <- i[[parm]]
        }
      }
      p
    }
  }
  observe({
    id <- rownames(FilterEnrich(kegg36, input$PVpath36)[input$pathviewTable36_rows_selected,])
    if (identical(id, character(0))) return()
    output$pathview36 <- renderPlot(
      expr = {
        i <- req(my_img36())
        r <- req(i$raster)
        plot(r)
      }, 
      width = img_dim_f36("w"),
      height = img_dim_f36("h")
    )
  })
  
  output$pathviewTable36 <- renderDataTable({
    datatable(FilterEnrich(kegg36, input$PVpath36), selection = 'single')
  })
  #Reactome
  output$Reactome36 <- renderPlot({
    PlotEnrichmentAnalysis(react36, input$PVR36)
  })
  output$tableReactome36 <- renderDataTable({
    datatable(FilterEnrich(react36, input$PVR36))
  })
  #Go Ontology
  observe({
    input$selectOnt36
    switch(
      input$selectOnt36,
      "BP" = {
        go36 <- goBP36
        output$Goenrich36 <- renderPlot({
          PlotEnrichmentAnalysis(go36, input$PVE36)
        })
        output$tableEnrich36 <- renderDataTable({
          datatable(FilterEnrich(go36, input$PVE36))
        })
      },
      "CC" = {
        go36 <- goCC36
        output$Goenrich36 <- renderPlot({
          PlotEnrichmentAnalysis(go36, input$PVE36)
        })
        output$tableEnrich36 <- renderDataTable({
          datatable(FilterEnrich(go36, input$PVE36))
        }) 
      },
      "MF" = {
        go36 <- goMF36
        output$Goenrich36 <- renderPlot({
          PlotEnrichmentAnalysis(go36, input$PVE36)
        })
        output$tableEnrich36 <- renderDataTable({
          datatable(FilterEnrich(go36, input$PVE36))
        }) 
      }
    )
  })

  #Venn Diagram
  output$vennDiagram <- renderPlot({
    PlotVennDiagram(input$PVvenn, input$LFCvenn)
  })
  output$vennDiagramTable <- renderDataTable({
    datatable(deg6[intersect(rownames(FilterDEGs(deg6, input$PVvenn, input$LFCvenn)), rownames(FilterDEGs(deg36, input$PVvenn, input$LFCvenn))),])
  })
}

shinyApp(ui, server)

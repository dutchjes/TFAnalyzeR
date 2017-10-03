

library(shiny)
library(reshape)
options(shiny.maxRequestSize = 100*1024^2)

ui <- fluidPage(
  
  titlePanel("Plotting Calibration Curves"),
  
    sidebarLayout(
      sidebarPanel(
        fileInput(inputId = "fileID", label = "Choose TraceFinder file:"),
        selectInput(inputId = "TargetCompound", label = "Choose target compound:", choices = c()),
        selectInput(inputId = "ISTD", label = "Choose internal standard:", choices = c()),
        selectInput(inputId = "CalLevels", label = "Calibration Levels", choices = c(), multiple = TRUE),
        selectInput(inputId = "CalModel", label = "Calibration Model", choices = c("X", "1/X", "1/X2"))
       
     ),
    
    mainPanel(
      plotOutput(outputId = "CalibrationCurve")
    )
  )
)


server <- function(input, output, session){
  
  TFdata <- eventReactive(input$fileID, {
      inFile <- input$fileID
      newFile <- read.csv(inFile$datapath, header = TRUE)
      subset(newFile, newFile$Peak.Label == "T1")
      
  })
  

  rv <- reactiveValues(targets = c(),
                       istds = c(),
                       cal.levels = c(),
                       target.area = c(),
                       istd.area = c(),
                       response.area = c(),
                       weights = c()
                       )
  
  observeEvent(input$fileID,{
    
    rv$targets <- unique(TFdata()[which(TFdata()$Type == "Target Compound"), 1])
    updateSelectInput(session, "TargetCompound", choices = rv$targets)
    
    rv$istds <- unique(TFdata()[which(TFdata()$Type == "Internal Standard"), 1])
    updateSelectInput(session, "ISTD", choices = rv$istds)
    
    rv$cal.levels <- unique(TFdata()[which(TFdata()$Sample.Type == "Cal Std"), "Level"])
    updateSelectInput(session, "CalLevels", choices = rv$cal.levels)
    
  })
  
  
  observeEvent(c(input$TargetCompound, input$CalLevels),{
    rv$target.area <- #as.numeric(as.character(unlist(
      subset(TFdata(), 
                                                            TFdata()[,1] == input$TargetCompound & 
                                                              TFdata()$Sample.Type == "Cal Std" & 
                                                              TFdata()$Level %in% input$CalLevels#, 
                                                           # select = c(Area)
             )
    #)))
    
    rv$target.area <- as.data.frame(cast(rv$target.area, Level~.,
                                         function(x) mean(na.omit(as.numeric(as.character(x)))),
                                         value = "Area"))[,2]
    
    rv$response.area <- rv$target.area/rv$istd.area
  })
  
  observeEvent(c(input$ISTD,input$CalLevels),{
    rv$istd.area <- #as.numeric(as.character(unlist(
      subset(TFdata(), 
                                                          TFdata()[,1] == input$ISTD &
                                                            TFdata()$Sample.Type == "Cal Std" &
                                                            TFdata()$Level %in% input$CalLevels#, 
                                                         # select = c(Area)
             )
    #)))
    
    rv$istd.area <- as.data.frame(cast(rv$istd.area, Level~.,
                                       function(x) mean(na.omit(as.numeric(as.character(x)))),
                                       value = "Area"))[,2]
    
    rv$response.area <- rv$target.area/rv$istd.area
  })
  
 observeEvent(c(input$CalModel, input$CalLevels),{
    if(input$CalModel == "X"){
      rv$weights <- as.numeric(as.character(input$CalLevels))
    }
    if(input$CalModel == "1/X"){
      rv$weights <- as.numeric(as.character(input$CalLevels))
      rv$weights <- 1/rv$weights
    }
   if(input$CalModel == "1/X2"){
     rv$weights <- as.numeric(as.character(input$CalLevels))
     rv$weights <- 1/(rv$weights^2)

   }

  })

  
  
  output$CalibrationCurve <- renderPlot({
    
    lm <- lm(rv$response.area ~ as.numeric(as.character(input$CalLevels)), weights = rv$weights
             )
    plot(x = input$CalLevels, rv$response.area, xlab = "Calibration Levels", ylab = "Response Area (Target/ISTD)")
    abline(lm)
    legend("topleft", legend = bquote("R"^2 ~ "=" ~ .(round(summary(lm)$r.squared, 2))))
      
  })

}


shinyApp(ui = ui, server = server)

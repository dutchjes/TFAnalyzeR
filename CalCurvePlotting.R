

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
        selectInput(inputId = "CalModel", label = "Calibration Model", choices = c("X", "1/X", "1/X2"), selected = "1/X")
       
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
  
  rv <- reactiveValues(targets = c(), #set based on TFdata, should contain all
                       istds = c(), #set based on TFdata, should contain all
                       cal.levels = c(), #set based on TFdata, should contain all
                       target.area = c(), #subset based on user selected target
                       istd.area = c(), #subset based on user selected istd
                       response.area = c(), #subset based on user selected target and istd
                       weights = c(), #calculated based on selected cal standards and cal model
                       take.level = c(), #subset based on user selected cal standards
                       data = data.frame(c())
                       )

    
## first, define contents of dropdown lists from TFdata contents
  observeEvent(input$fileID,{
    
    rv$targets <- unique(TFdata()[which(TFdata()$Type == "Target Compound"), 1])
    updateSelectInput(session, "TargetCompound", choices = rv$targets)
    
    rv$istds <- unique(TFdata()[which(TFdata()$Type == "Internal Standard"), 1])
    updateSelectInput(session, "ISTD", choices = rv$istds)
    
    rv$cal.levels <- unique(TFdata()[which(TFdata()$Sample.Type == "Cal Std"), "Level"])
    updateSelectInput(session, "CalLevels", choices = rv$cal.levels)
    
  })
 
## second, calculate 
  observeEvent(c(input$TargetCompound, input$CalLevels, input$ISTD),{
    rv$target.area <- #as.numeric(as.character(unlist(
      subset(TFdata(), 
             TFdata()[,1] == input$TargetCompound & 
             TFdata()$Sample.Type == "Cal Std" & 
             TFdata()$Level %in% input$CalLevels, 
             select = c(Area, Filename, Level)
             )

    rv$istd.area <- #as.numeric(as.character(unlist(
      subset(TFdata(), 
             TFdata()[,1] == input$ISTD &
               TFdata()$Sample.Type == "Cal Std" &
               TFdata()$Level %in% input$CalLevels, 
             select = c(Area, Filename, Level)
      )
    
    rv$data <- merge(rv$target.area, rv$istd.area, by = "Filename")
    rv$take.level <- as.numeric(as.character(rv$data$Level.x))
    rv$response.area <- as.numeric(as.character(rv$target.area$Area))/as.numeric(as.character(rv$istd.area$Area))
  })
  

 observeEvent(c(input$CalModel, input$CalLevels),{
    if(input$CalModel == "X"){
      rv$weights <- as.numeric(as.character(rv$take.level))
    }
    if(input$CalModel == "1/X"){
      rv$weights <- as.numeric(as.character(rv$take.level))
      rv$weights <- 1/rv$weights
    }
   if(input$CalModel == "1/X2"){
     rv$weights <- as.numeric(as.character(rv$take.level))
     rv$weights <- 1/(rv$weights^2)

   }

  })

  output$CalibrationCurve <- renderPlot({
    
    lm <- lm(as.numeric(as.character(rv$response.area)) ~ as.numeric(as.character(rv$take.level)), weights = rv$weights
             )
    plot(x = as.numeric(as.character(rv$take.level)), 
         as.numeric(as.character(rv$response.area)), xlab = "Calibration Levels", ylab = "Response Area (Target/ISTD)")
    abline(lm)
    legend("topleft", legend = bquote("R"^2 ~ "=" ~ .(round(summary(lm)$r.squared, 3))))
      
  })

}


shinyApp(ui = ui, server = server)

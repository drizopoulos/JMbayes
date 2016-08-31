# Define UI for miles per gallon application
shinyUI(fluidPage(
    
    # Application title
    headerPanel("Dynamic Predictions using Linear Mixed Effects (LME) models"),
    
    sidebarPanel(
        wellPanel(
            fileInput('RDfile', 'Load the R Workspace with the fitted LME model',
                      accept = NULL),
            
            fileInput('patientFile', 'Load subject data',
                      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

            fluidRow(column(4, radioButtons('sep', 'Separator', c(Comma = ',', Semicolon = ';', Tab = '\t'), ',')),
                     column(4, radioButtons('dec', 'Decimal', c(Dot = '.', Comma = ','), '.')),
                     column(4, radioButtons('quote', 'Quote', c(None = '', 'Double Quote' = '"', 'Single Quote' = "'"), '"')))
        ),
        
        wellPanel(uiOutput("timeVarChoose")),
        
        wellPanel(
            uiOutput("obsChoose"),
            fluidRow(column(4, radioButtons("interval", "Interval Type:", 
                                            c("confidence", "prediction"))),
                     column(4, checkboxInput("marginal", label = "Marginal Prediction", value = FALSE)),
                     column(4, numericInput("M", "Monte Carlo samples:", 500)))
        )
    ),
    
    mainPanel(
        tabsetPanel(
          tabPanel("Data", dataTableOutput('data')),
          tabPanel("Plot", plotOutput('plot'))
        )
    )
))

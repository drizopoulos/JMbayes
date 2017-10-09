# Define UI for miles per gallon application
shinyUI(fluidPage(
    
    # Application title
    headerPanel("Dynamic Predictions using Joint Models"),
    
    sidebarPanel(
        wellPanel(
            fileInput('RDfile', 'Load the R Workspace with the fitted joint model',
                      accept = NULL),
            
            fileInput('patientFile', 'Load subject data',
                      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

            fluidRow(column(4, radioButtons('sep', 'Separator', c(Comma = ',', Semicolon = ';', Tab = '\t'), ',')),
                     column(4, radioButtons('dec', 'Decimal', c(Dot = '.', Comma = ','), '.')),
                     column(4, radioButtons('quote', 'Quote', c(None = '', 'Double Quote' = '"', 'Single Quote' = "'"), '"')))
        ),
        
        wellPanel(
            uiOutput("predictEventTime")
        ),
        
        wellPanel(
            uiOutput("modelChoose"),
            
            uiOutput("outcomeChoose"),
            
            uiOutput("obsChoose"),
                        
            fluidRow(column(4, radioButtons('TypePlot', 'Type of Plot', 
                                            c("Survival" = 'surv', "Cumulative Incidence" = "cumInc",
                                              "Smiley Faces" = "smFace",
                                              "100 Clones" = 'stickMan', "Longitudinal" = 'longitudinal'), 
                                            "surv")),
                     column(4, numericInput("windowTime", "Target window time:", NULL)),
                     column(4,  numericInput("time", "Target horizon time:", NULL))),
            
            uiOutput("lastTime"),
            
            numericInput("M", "Monte Carlo samples:", 200),
                   
            fluidRow(column(6, downloadButton('downloadData', 'Download Event-free Probabilities')),
                     column(6, downloadButton('downloadPlot', 'Download Plot')))
            
        )
    ),
    
    mainPanel(
        tabsetPanel(
            tabPanel("Data", tableOutput('contents'), uiOutput("message")),
            tabPanel("Event-free Probabilities", tableOutput('survprobs'), uiOutput("message2")),
            tabPanel("Plot", plotOutput('plot')),
            tabPanel("Help", htmlOutput('help'))
        )
    )
))

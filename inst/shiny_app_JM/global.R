library("shiny")
library("JMbayes")
library("splines")
options(shiny.maxRequestSize = 1000*1024^2)

actionButton <- function (inputId, label, style = "" , additionalClass = "") {
    if (style %in% c("primary","info","success","warning","danger","inverse","link")) {
        class.style <- paste("btn",style,sep="-")
    } else class.style = ""
    
    tags$button(id = inputId, type = "button", 
                class = paste("btn action-button", class.style, additionalClass), label)
}

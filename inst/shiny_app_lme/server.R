shinyServer(function(input, output) {
    loadModel <- reactive({
        if (!is.null(input$RDfile)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- ls()
            ss <- sapply(objs, function (o) class(get(o))) == "lme"
            if (all(!ss)) {
                stop("\nIt seems that there is no model fitted by lme() ", 
                     "in the workspace that you loaded ...")
            }
            if (sum(ss) > 1) {
                warning("\nIt seems that in the loaded workspace more that one mixed ",
                        "models are present. The following calculations are based on ",
                        "model ", objs[ss][1L])
            }
            get(objs[ss][1L])
        }
    })
    
    loadPatient <- reactive({
        if (!is.null(input$patientFile)) {
            inFile <- input$patientFile
            read.csv(inFile$datapath, sep = input$sep, quote = input$quote, 
                     dec = input$dec)
        }
    })
    
    output$obsChoose <- renderUI({
        if (!is.null(input$patientFile) && input$timeVar != "<select>") {
            nd <- loadPatient()
            nr <- nrow(nd)
            if (!is.null(nr) && nr > 1) {
                sliderInput("obs", "Number of observations to use in prediction:", 
                            min = 1, max = nr, value = 1, step = 1,
                            animate = animationOptions(loop = TRUE))
            }
        }
    })
    
    output$timeVarChoose <- renderUI({
        if (!is.null(input$RDfile)) {
            lmeObject <- loadModel()
            selectInput("timeVar", "Select Time variable:", 
                        choices = c("<select>", all.vars(formula(lmeObject))))
        }
    })
    
    output$data <- renderDataTable({
        loadPatient()
    })
    
    calculate_preds <- reactive({
        if (!is.null(input$patientFile) && !is.null(input$timeVar) 
            && input$timeVar != "<select>") {
            nd <- loadPatient()
            lmeObject <- loadModel()
            for (colname in colnames(nd)) {
                if (is.factor(nd[[colname]])) {
                    nd[[colname]] <- factor(nd[[colname]], 
                                            levels = levels(lmeObject$data[[colname]]))
                }
            }
            nr <- nrow(nd)
            DFs <- vector("list", nr)
            for (i in seq_len(nr)) {
                DFs[[i]] <- IndvPred_lme(lmeObject, newdata = nd[1:i,], timeVar = input$timeVar,
                                   M = input$M, interval = input$interval, return_data = TRUE)
            }
            nd <- rbind(nd, nd[nrow(nd), ])
            nd[nrow(nd), input$timeVar] <- max(lmeObject$data[[input$timeVar]])
            nd$marginal <- predict(lmeObject, newdata = nd, level = 0)
            c(DFs, list(nd = nd), respVar = as.character(formula(lmeObject))[2L])
        }
    })
    
    output$plot <- renderPlot({
        if (!is.null(input$patientFile) && !is.null(input$timeVar) 
            && input$timeVar != "<select>") {
            DFs <- calculate_preds()
            times <- DFs[[input$obs]][[input$timeVar]]
            na_ind <- is.na(DFs[[input$obs]][["low"]])
            last_time <- times[!na_ind][1L]
            active_DF <- DFs[[input$obs]]
            nd <- DFs[["nd"]]
            nd_active <- nd[seq_len(input$obs), ]
            matplot(active_DF[[input$timeVar]], active_DF[c("pred", "low", "upp")], 
                    type = "l", col = c(2, 1, 1), lty = c(1, 2, 2), lwd = 2,
                    xlab = input$timeVar, ylab = DFs$respVar)
            abline(v = last_time, lty = 2)
            form <- as.formula(paste(DFs$respVar, "~", input$timeVar))
            points(form, data = nd_active, pch = 8, col = 4)
            if (input$marginal) {
                form <- as.formula(paste("marginal ~", input$timeVar))
                lines(form, data = nd, col = 3, lty = 2, lwd = 2)
            }
        }
    })
})
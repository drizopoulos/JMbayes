shinyServer(function(input, output) {
    
    findJMbayesObjs <- reactive({
        if (!is.null(input$RDfile)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- ls()
            ss <- sapply(objs, function (o) class(get(o))) %in% c("JMbayes", "mvJMbayes")
            if (all(!ss)) {
                stop("\nIt seems that there is no joint model fitted by jointModelBayes() ", 
                     "in the workspace that you loaded ...")
            }
            out <- objs[ss]
            names(out) <- out
            out
        }
    })
    
    output$modelChoose <- renderUI({
        if (!is.null(input$RDfile)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- findJMbayesObjs()
            selectInput("model", "Choose joint model:", objs, objs[1L])
        }
    })
    
    loadObject <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- findJMbayesObjs()
            choice <- input$model
            get(objs[input$model])
        }
    })
    
    output$outcomeChoose <- renderUI({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            object <- loadObject()
            if (inherits(object, "mvJMbayes")) {
                components <- object$model_info$mvglmer_components
                outcomes <- components[grep("respVar", names(components), fixed = TRUE)]
                names(outcomes) <- unlist(outcomes, use.names = FALSE)
                outcomes[] <- lapply(seq_along(outcomes), function (i) i)
                selectInput("outcome", "Choose outcome:", outcomes, outcomes[1L])
            }
        }
    })
        
    dataObject <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            object <- loadObject()
            if (inherits(object, "JMbayes")) {
                forms <- object$Forms[c("formYx", "formYz", "formT")]
                forms$formT <- reformulate(attr(delete.response(object$Terms$termsT), "term.labels"))
                nams <- unique(unlist(lapply(forms, all.vars), use.names = FALSE))
                d <- object$Data$data[1:3, c(nams)]
            } else {
                mvglmer_components <- object$model_info$mvglmer_components
                Terms_comp <- mvglmer_components[grep("Term", names(mvglmer_components))]
                mvglmer_vars <- unlist(lapply(Terms_comp, all.vars), use.names = FALSE)
                cox_vars <- all.vars(delete.response(object$model_info$coxph_components$Terms))
                nams <- unique(c(mvglmer_vars, cox_vars))
                d <- object$model_info$mvglmer_components$data[1:3, c(nams)]
            }
            d[] <- lapply(d, function (x) {x[] <- NA; x})
            d
        }
    })
    
    ND <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            if (is.null(input$patientFile)) {
                dataObject()
            } else {
                d <- dataObject()
                indF <- sapply(d, is.factor)
                namsFactors <- names(d[indF])
                levelsF <- lapply(d[indF], levels) 
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, 
                                             quote=input$quote, dec = input$dec), TRUE)
                if (!inherits(tt, "try-error")) {
                    f <- function (f, l) factor(f, levels = l)
                    inData[namsFactors] <- mapply(f, inData[namsFactors], 
                                                  levelsF[namsFactors], SIMPLIFY = FALSE)
                    inData <- inData[names(inData) %in% names(d)]
                    inData$id <- rep(1, nrow(inData))
                    inData
                }
            }
        }
    })
    
    output$obsChoose <- renderUI({
        if (!is.null(input$patientFile)) {
            nr <- nrow(ND())
            if (!is.null(nr) && nr > 1) {
                sliderInput("obs", "Number of observations to use in prediction:", 
                            min = 1, max = nr, value = 1, step = 1,
                            animate = animationOptions(loop = TRUE))
            }
        }
    })
    
    output$lastTime <- renderUI({
        if (!is.null(input$patientFile)) {
            object <- loadObject()
            nd <- ND()
            times <- if (inherits(object, "JMbayes")) nd[[object$timeVar]] else nd[[object$model_info$timeVar]]
            if (!is.null(times))
                numericInput("lasttime", "Last time point without event:", 
                             round(max(times, na.rm = TRUE), 2))
        }
    })
    
    output$message <- renderPrint({
        if (is.null(input$RDfile)) {
            cat("<br /><h3> <span style='color:black'> Welcome to the web interface",
                  "for producing dynamic predictions from joint models using package",
                  "<a href='http://cran.r-project.org/package=JMbayes'",
                  "style='text-decoration:none;' target='_blank'><b>JMbayes</b></a></h3>",
                "<br /><br /><h4>Use the menus on the left to load the R workspace containing the fitted joint model",
                "to continue ... (for further details & instructions check the 'Help' tab above)</h4>")
        } else {
            if (is.null(input$patientFile)) {
                dataset <- ND()
                classes <- sapply(dataset, class)
                classes[classes == "numeric" | classes == "integer"] <- "a numeric (continuous) variable."
                ind <- classes == "factor"
                classes[ind] <- paste("a factor (categorical) variable, with levels", 
                    sapply(dataset[ind], function (x) paste0("<i>", levels(x), "</i>", collapse = ", ")))
                msg1 <- "The data of the new subject should be stored in a CSV file with columns:"
                msg2 <- "<ul style='list-style-type:circle'>"
                msg3 <- paste("<li><b>", names(dataset), "</b>:", classes, "</li>", collapse = " ")
                msg4 <- "</ul>"
                msg5 <- "<br />You can use as a template the table above (e.g., copy-paste it in Excel and save it as CSV)." 
                msg6 <- "<u>Note:</u> R is case sensitive"
                cat("<br /><br /><p>", paste(msg1, msg2, msg3, msg4, msg5, msg6), "</p>")
            } else {
                d <- dataObject()
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, 
                                             colClasses = sapply(d, class), dec = input$dec), TRUE)
                if (inherits(tt, "try-error"))
                    cat("<br /><span style='color:red'>Something went wrong when loading the data of the new subject;</span>",
                        "<br>try tweaking the 'Separator' and 'Decimal' options on the left panel ...")
            }
        } 
    })
    
    output$message2 <- renderPrint({
        if (is.null(input$RDfile) || is.null(input$patientFile))
            cat("<br />No subject data have been loaded yet ...")
    })
    
    output$help <- renderPrint({
        cat("<br /><h3>Dynamic Predictions in a Nutshell</h3>")
        cat("<br />The dynamic predictions calculated by this app are applicable in the ", 
            "context of follow-up studies where the aim is to utilize longitudinally ", 
            "measured outcomes to predict future values of these outcomes or future ", 
            "events for the sample units. For example, a cohort of patients is followed-up ",
            "in time and it is of interest to predict events, such as relapse of disease ",
            "or death using longitudinally measured biomarkers.", "These predictions ", 
            "are derived under the framework of joint models for longitudinal and ", 
            "time-to-event data, and more technical details can be found in the ", 
            "References.")
        #############
        cat("<br /><br /><h3>How this App Works</h3>")
        cat("<br />On the left hand side there are (or the will appear) a number of ", 
            "control arguments.", "In the main panel there are three tabs, named, <span style='color:blue'>'Data'</span>,
            <span style='color:blue'>'Event-free Probabilities'</span> and <span style='color:blue'>'Plots'</span>, in which the results appear.")
        p1 <- paste("Load the R workspace that contains the joint model(s) based on which", 
                    "you would like to calculated the predictions.")
        p21 <- paste("Load the comma separated file (.csv) containing the available data", 
                     "including both the longitudinally measured outcomes and potentially",
                     "other baseline covariates, of the sample unit (e.g., patient) for", 
                     "which predictions are to be calculated.")
        p22 <- paste("The names of the variables/columns in this file and the coding of", 
                     "the categorical variables must be the same as in the database used", 
                     "to fit the joint model. More exact information on how this file", 
                     "should look like appears in the <span style='color:blue'>'Data'</span> tab after you load the R", 
                     "workspace in 'Step I'.")
        p23 <- paste("Depending on the locale of your machine it may be required to tweak", 
                     "the options <span style='color:red'>'Separator'</span>,", 
                     "<span style='color:red'>'Decimal'</span> and <span style='color:red'>'Quote'</span>.")
        p31 <- paste("After the R workspace has been loaded a select dialog box appears", 
                     "on the left hand side that allows the user to <span style='color:red'>choose the joint", 
                     "model</span> contained in the workspace based on which the predictions", 
                     "will be calculated. In addition, in the main panel, under the", 
                     "<span style='color:blue'>'Data'</span> tab, an example of the .csv file that the user should load", 
                     "appears, that gives the names of the columns that should be used", 
                     "and the coding for categorical variables.")
        p32 <- paste("After the user has loaded the data of the 'new' sample unit its", 
                     "data are depicted under the <span style='color:blue'>'Data'</span> tab.", 
                     "The calculation of the", 
                     "prediction is initiated after the user selects either the", 
                     "<span style='color:blue'>'Event-free Probabilities'</span> or", 
                     "<span style='color:blue'>'Plot'</span> tab. In addition, in the left", 
                     "hand side a slider appears, named", 
                     "<span style='color:red'>'Number of observations to use in prediction'</span>,", 
                     "which selects the number of measurements", 
                     "used in the prediction. By pressing the small play button at the", 
                     "end of this slider the dynamic predictions are depicted in the", 
                     "'Plot' tab.")
        p33 <- paste("The <span style='color:red'>'Type of Plot'</span> options on the left control the type of appearing",
                     "in the 'Plot' tab. The default option <span style='color:#FE2EF7'>'Survival'</span> depicts survival", 
                     "probabilities (aka event-free probabilities) for the event of", 
                     "interest based on the avalable longitudinal measurements.", 
                     "Similarly, option <span style='color:#FE2EF7'>'Cumulative Incidence'</span> depicts cumulative", 
                     "incidence probabilities for the event of interest. Options", 
                     "<span style='color:#FE2EF7'>'Smiley Faces'</span> and", 
                     "<span style='color:#FE2EF7'>'100 Clones'</span> are used in conjuction with the", 
                     "<span style='color:red'>'Target Window Time'</span> option below and depict, in a simple manner, the risk that the event",
                     "will occur within the time interval <i>(t, t + Dt]</i>, where <i>t</i> denotes",
                     "the time of the last available longitudinal measurement (specified", 
                     "according to the slider 'Number of observations to use in prediction')", 
                     "and <i>Dt</i> is given in the 'Target Window Time'. The idea is that if at",
                     "the particular time point t we would create 100 clones of the subject at hand", 
                     "the 'Smiley Faces' and '100 Clones' plots show how many of them would die within the", 
                     "time interval from <i>t</i> up to <i>t + Dt</i>.",
                     "Finally, option <span style='color:#FE2EF7'>'Longitudinal'</span>", 
                     "can be selected to depict predictions for the longitudinal outcome based", 
                     "on the recored information.")
        p34 <- paste("The field <span style='color:red'>'Target horizon time'</span> specifies a horizon time point at", 
                     "which predictions are of interest. For example, 10-year survival probability.")
        p35 <- paste("The field <span style='color:red'>'Monte Carlo samples'</span> specifies the number of Monte Carlo", 
                     "samples used in calculating the predictions-more details can be found", 
                     "in the References.")
        cat("<dl><dt>Step I</dt>", "<dd>- ", p1, "</dd>", "<dd>  </dd>",
            "<dt>Step II</dt>", "<dd>- ", p21, "</dd>", "<dd>- ", p22, "</dd>", "<dd>- ", p23, "</dd>", "<dd>  </dd>",
            "<dt>Options & Widgets</dt>", "<dd>- ", p31, "</dd>", "<dd>- ", p32, "</dd>", 
            "<dd>- ", p33, "</dd>", "<dd>- ", p34, "</dd>", "<dd>- ", p35, "</dd>", "</dl>")
        #############
        cat("<br /><h3>References</h3>")
        p1 <- paste0("<li>", "Yu M, Taylor J, Sandler H. (2008). ",
                     "Individualized prediction in prostate cancer studies using a joint longitudinal-survival-cure model. ",
                     "<a href='http://dx.doi.org/10.1198/016214507000000400' ", 
                     "style='text-decoration:none;' target='_blank'>", 
                     "<i>Journal of the American Statistical Association</i> <b> 103</b>, 178-187.</a>",
                     "</li>")
        p2 <- paste0("<li>", "Rizopoulos D. (2011). Dynamic predictions and prospective ",
                     "accuracy in joint models for longitudinal and time-to-event data. ",
                     "<a href='http://dx.doi.org/10.1111/j.1541-0420.2010.01546.x' ", 
                     "style='text-decoration:none;' target='_blank'>", 
                     "<i>Biometrics</i> <b> 67</b>, 819-829.</a>",
                     "</li>")
        p3 <- paste0("<li>", "Rizopoulos D. (2012). ",
                     "<a href='http://www.crcpress.com/product/isbn/9781439872864' ", 
                     "style='text-decoration:none;' target='_blank'>", 
                     "<i>Joint Models for Longitudinal and Time-to-Event Data, with Applications in R</i></a>.",
                     " Boca Raton: Chapman & Hall/CRC.",
                     "</li>")
        p4 <- paste0("<li>", "Taylor J, Park Y, Ankerst D, Proust-Lima C, ", 
                     "Williams S, Kestin L, Bae K, Pickles T, Sandler H. (2013). ",
                     "Real-time individual predictions of prostate cancer recurrence", 
                     "using joint models. ", "<a href='http://dx.doi.org/10.1111/j.1541-0420.2012.01823.x' ", 
                     "style='text-decoration:none;' target='_blank'>", 
                     "<i>Biometrics</i> <b> 69</b>, 206-213.</a>",
                     "</li>")
        cat("<ol>", p1, p2, p3, p4, "</ol> ")
    })
    
    sfits <- reactive({
        if (!is.null(input$patientFile) && input$TypePlot != 'longitudinal') {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            lastTimeData <- if (inherits(object, "JMbayes")) max(nd[[object$timeVar]]) 
            else max(nd[[object$model_info$timeVar]])
            sfits <- vector("list", n)
            for (i in 1:n) {
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], M = input$M, 
                                        last.time = lt)
            }
            sfits
        }
    })
    
    lfits <- reactive({
        if (!is.null(input$patientFile) && input$TypePlot == 'longitudinal') {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            lastTimeData <- if (inherits(object, "JMbayes")) max(nd[[object$timeVar]]) 
            else max(nd[[object$model_info$timeVar]])
            lfits <- vector("list", n)
            for (i in 1:n) {
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                lfits[[i]] <- predict(object, newdata = nd[1:i, ], M = input$M, 
                                      returnData = TRUE, type = "Subject", 
                                      interval = "prediction", last.time = lt)
            }
            lfits
        }
    })
    
    sfits2 <- reactive({
        if (!is.na(input$time) || !is.na(input$windowTime)) {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            sfits <- vector("list", n)
            for (i in 1:n) {
                timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                lastTimeData <- max(nd[1:i, timeVar])
                target.time <- if (!is.na(input$windowTime)) lastTimeData + input$windowTime else input$time
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData && !is.na(input$windowTime))
                    target.time <- lastTimeUser + input$windowTime
                if (target.time > lastTimeData)
                    sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], 
                                            survTimes = target.time, M = input$M,
                                            last.time = lt)
            }
            sfits
        }
    })
    
    output$contents <- renderTable({
        if (!is.null(input$RDfile)) {
            nd <- ND()
            if (!is.null(input$patientFile)) {
                tt <- try(nd[, head(1:ncol(nd), -1)], TRUE)
                if (!inherits(tt, "try-error"))
                    nd <- tt
            }
            nd
        }
    })
    
    sprobs <- reactive({
        if (!is.null(input$patientFile)) {
            if (is.na(input$time) && is.na(input$windowTime)) {
                sfits. <- sfits()
            } else {
                sfits. <- sfits2()
            }
            nn <- if(is.na(input$obs)) length(sfits.) else input$obs
            ss <- sfits.[[nn]]
            if (!is.null(ss)) {
                sf <- sfits.[[nn]]
                f <- function (d, t) {
                    dd <- d[1, , drop = FALSE]
                    dd[1, ] <- c(as.vector(t), rep(1, ncol(dd) - 1))
                    round(rbind(dd, d), 4)
                }
                d <- mapply(f, sf$summaries, sf$last.time, SIMPLIFY = FALSE)
                d <- d[[1L]][, c(1,2,4,5), drop = FALSE]
                colnames(d) <- c("time", "Surv", "95% low", "95% upp")
                nr <- nrow(d)
                out <- if (nr > 10) d[round(seq(1, nr, length.out = 10)), ] else d
                rownames(out) <- NULL
                out
            } else NULL
        }
    })
    
    output$survprobs <- renderTable({
        sprobs()
    })
    
    output$plot <- renderPlot({
        if (!is.null(input$patientFile)) {
            if (input$TypePlot != 'longitudinal') {
                if (!is.na(input$windowTime) && (input$TypePlot == "stickMan" 
                                                          || input$TypePlot == "smFace")) {
                    stickMan <- input$TypePlot == "stickMan"
                    sfits. <- sfits2()
                    nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                    ss <- round(100 * (1 - sfits.[[nn]][["summaries"]][[1]][1, 2]))
                    draw.stick <- JMbayes:::draw.stick
                    smilyface <- JMbayes:::smilyface
                    if (stickMan) {
                        xx <- seq(0.2, 1.1, 0.1)
                        yy <- seq(0.0, 1.8, 0.2)
                    } else {
                        xx <- seq(0.3, 0.8, len = 10)
                        yy <- seq(0.45, 1.8, len = 10)
                    }
                    cords <- data.matrix(expand.grid(xx, rev(yy)))
                    cols <- rep("blue", 100) 
                    cols[seq_len(ss)] <- "red"
                    op <- par(mar = c(0, 0, 2.1, 0))
                    plot(c(.25, 1.25), c(0, 2), type = "n", xaxt = 'n', yaxt = 'n', ann = FALSE)
                    title(paste0("Risk = ", ss, "%"), cex.main = 1.7)
                    for (cc in 1:100) {
                        if (stickMan) {
                            draw.stick(cords[cc, 1], cords[cc, 2], linecol = cols[cc], scale = 0.2)
                        } else {
                            m <- if (cols[cc] == "blue") "happy" else "sad"
                            smilyface(cords[cc, ], 0.01, mood = m)
                        }
                    }
                    par(op)
                } else {
                    sfits. <- sfits()
                    nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                    if (input$TypePlot == "surv") {
                        if (inherits(sfits.[[nn]], "survfit.mvJMbayes")) {
                            plot(sfits.[[nn]], which_outcomes = as.numeric(input$outcome))
                        } else {
                            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, 
                                 fill.area = TRUE, include.y = TRUE, lwd = 2, ask = FALSE, 
                                 cex = 2, main = "")                
                        }
                    }
                    if (input$TypePlot == "cumInc") {
                        if (inherits(sfits.[[nn]], "survfit.mvJMbayes")) {
                            plot(sfits.[[nn]], which_outcomes = as.numeric(input$outcome),
                                 fun = function (s) 1 - s)
                        } else {
                            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, 
                                 fill.area = TRUE, include.y = TRUE, lwd = 2, ask = FALSE,
                                 cex = 2, main = "", fun = function (s) 1 - s, 
                                 ylab = "Cumulative Incidence")
                        }
                    }
                    if (!is.na(input$windowTime) || !is.na(input$time)) {
                        object <- loadObject()
                        nd <- ND()
                        nr <- nrow(nd)
                        lastTimeUser <- input$lasttime
                        timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                        lastTimeData <- max(nd[1:nn, timeVar])
                        target.time <- if (nn == nr && !is.na(lastTimeUser) && 
                                               lastTimeUser > lastTimeData) {
                            if (!is.na(input$windowTime)) lastTimeUser + input$windowTime else input$time
                        } else {
                            if (!is.na(input$windowTime)) lastTimeData + input$windowTime else input$time
                        }
                        abline(v = target.time, lty = 2, col = 2, lwd = 2)
                    }
                }
            } else {
                object <- loadObject()
                lfits. <- lfits()
                require("lattice")
                nn <- if(is.na(input$obs)) length(lfits.) else input$obs
                timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                resp <- paste(object$Forms$formYx)[2L]
                tv <- lfits.[[nn]][[timeVar]]
                lastTimeData <- with(lfits.[[nn]], tv[!is.na(low)][1])
                lastTimeUser <- input$lasttime
                lastTime <- if (nn == length(lfits.)) 
                    max(lastTimeData, lastTimeUser, na.rm = TRUE) else lastTimeData
                form <- as.formula(paste("pred + low + upp +", resp, "~", timeVar))
                yl <- range(c(data.matrix(do.call(rbind, lfits.)[c("pred", "low", "upp")])),
                            na.rm = TRUE)
                yl <- yl + c(-0.05, 0.05) * yl
                target.time <- if (!is.na(input$windowTime) || !is.na(input$time)) {
                    object <- loadObject()
                    nd <- ND()
                    nr <- nrow(nd)
                    lastTimeUser <- input$lasttime
                    timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                    lastTimeData <- max(nd[1:nn, timeVar])
                    if (nn == nr && !is.na(lastTimeUser) && 
                                           lastTimeUser > lastTimeData) {
                        if (!is.na(input$windowTime)) lastTimeUser + input$windowTime else input$time
                    } else {
                        if (!is.na(input$windowTime)) lastTimeData + input$windowTime else input$time
                    }
                } else NA
                print(xyplot(form, data = lfits.[[nn]], last.time = lastTime, target.time = target.time,
                             panel = function (..., last.time, target.time) {
                                 xx <- ..1; yy <- ..2; ind <- ..4
                                 xx <- do.call(cbind, split(xx, ind))
                                 yy <- do.call(cbind, split(yy, ind))
                                 na.ind <- is.na(yy[, 2])
                                 lx <- xx[!na.ind, 2]; ll <- yy[!na.ind, 2]; lu <- yy[!na.ind, 3]
                                 lpolygon(c(lx, rev(lx)), c(ll, rev(lu)), border = "transparent", 
                                          col = "lightgrey")
                                 panel.xyplot(xx[, 1], yy[, 1], type = "l", lty = 1, col = 2, lwd = 2)
                                 panel.xyplot(xx[, 2], yy[, 2], type = "l", lty = 2, col = 1, lwd = 2)
                                 panel.xyplot(xx[, 3], yy[, 3], type = "l", lty = 2, col = 1, lwd = 2)
                                 panel.xyplot(xx[na.ind, 4], yy[na.ind, 4], type = "p", pch = 8, 
                                              cex = 1.1, col = 1)
                                 panel.abline(v = last.time, lty = 3)
                                 if (!is.na(target.time))
                                     panel.abline(v = target.time, lty = 2, col = 2, lwd = 2)
                             }, xlab = "Time", ylim = yl,
                             ylab = paste("Predicted", resp)))
            }
        }
    })
        
    output$downloadData <- downloadHandler(
        filename = function () { paste(input$dataset, '.csv', sep = '') },
        content = function (file) {
            write.table(sprobs(), file, sep = input$sep, dec = input$dec,
                        row.names = FALSE)
        }
    )
    
    output$downloadPlot <- downloadHandler(
        filename = function () { paste(input$dataset, '.png', sep = '') },
        content = function (file) {
            png(file)
            #sfits. <- sfits()
            #nn <- if(is.na(input$obs)) length(sfits.) else input$obs
            #plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
            #     include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "")
            if (!is.null(input$patientFile)) {
                if (input$TypePlot != 'longitudinal') {
                    if (!is.na(input$windowTime) && (input$TypePlot == "stickMan" 
                                                     || input$TypePlot == "smFace")) {
                        stickMan <- input$TypePlot == "stickMan"
                        sfits. <- sfits2()
                        nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                        ss <- round(100 * (1 - sfits.[[nn]][["summaries"]][[1]][1, 2]))
                        draw.stick <- JMbayes:::draw.stick
                        smilyface <- JMbayes:::smilyface
                        if (stickMan) {
                            xx <- seq(0.2, 1.1, 0.1)
                            yy <- seq(0.0, 1.8, 0.2)
                        } else {
                            xx <- seq(0.3, 0.8, len = 10)
                            yy <- seq(0.45, 1.8, len = 10)
                        }
                        cords <- data.matrix(expand.grid(xx, rev(yy)))
                        cols <- rep("blue", 100) 
                        cols[seq_len(ss)] <- "red"
                        op <- par(mar = c(0, 0, 2.1, 0))
                        plot(c(.25, 1.25), c(0, 2), type = "n", xaxt = 'n', yaxt = 'n', ann = FALSE)
                        title(paste0("Risk = ", ss, "%"), cex.main = 1.7)
                        for (cc in 1:100) {
                            if (stickMan) {
                                draw.stick(cords[cc, 1], cords[cc, 2], linecol = cols[cc], scale = 0.2)
                            } else {
                                m <- if (cols[cc] == "blue") "happy" else "sad"
                                smilyface(cords[cc, ], 0.01, mood = m)
                            }
                        }
                        par(op)
                    } else {
                        sfits. <- sfits()
                        nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                        if (input$TypePlot == "surv") {
                            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                                 include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "")                
                        }
                        if (input$TypePlot == "cumInc") {
                            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                                 include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "",
                                 fun = function (s) 1 - s, ylab = "Cumulative Incidence")                
                        }
                        if (!is.na(input$windowTime) || !is.na(input$time)) {
                            object <- loadObject()
                            nd <- ND()
                            nr <- nrow(nd)
                            lastTimeUser <- input$lasttime
                            timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                            lastTimeData <- max(nd[1:nn, timeVar])
                            target.time <- if (nn == nr && !is.na(lastTimeUser) && 
                                                   lastTimeUser > lastTimeData) {
                                if (!is.na(input$windowTime)) lastTimeUser + input$windowTime else input$time
                            } else {
                                if (!is.na(input$windowTime)) lastTimeData + input$windowTime else input$time
                            }
                            abline(v = target.time, lty = 2, col = 2, lwd = 2)
                        }
                    }
                } else {
                    object <- loadObject()
                    lfits. <- lfits()
                    require("lattice")
                    nn <- if(is.na(input$obs)) length(lfits.) else input$obs
                    timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                    resp <- paste(object$Forms$formYx)[2L]
                    tv <- lfits.[[nn]][[timeVar]]
                    lastTimeData <- with(lfits.[[nn]], tv[!is.na(low)][1])
                    lastTimeUser <- input$lasttime
                    lastTime <- if (nn == length(lfits.)) 
                        max(lastTimeData, lastTimeUser, na.rm = TRUE) else lastTimeData
                    form <- as.formula(paste("pred + low + upp +", resp, "~", timeVar))
                    yl <- range(c(data.matrix(do.call(rbind, lfits.)[c("pred", "low", "upp")])),
                                na.rm = TRUE)
                    yl <- yl + c(-0.05, 0.05) * yl
                    target.time <- if (!is.na(input$windowTime) || !is.na(input$time)) {
                        object <- loadObject()
                        nd <- ND()
                        nr <- nrow(nd)
                        lastTimeUser <- input$lasttime
                        timeVar <- if (inherits(object, "JMbayes")) object$timeVar else object$model_info$timeVar
                        lastTimeData <- max(nd[1:nn, timeVar])
                        if (nn == nr && !is.na(lastTimeUser) && 
                                lastTimeUser > lastTimeData) {
                            if (!is.na(input$windowTime)) lastTimeUser + input$windowTime else input$time
                        } else {
                            if (!is.na(input$windowTime)) lastTimeData + input$windowTime else input$time
                        }
                    } else NA
                    print(xyplot(form, data = lfits.[[nn]], last.time = lastTime, target.time = target.time,
                                 panel = function (..., last.time, target.time) {
                                     xx <- ..1; yy <- ..2; ind <- ..4
                                     xx <- do.call(cbind, split(xx, ind))
                                     yy <- do.call(cbind, split(yy, ind))
                                     na.ind <- is.na(yy[, 2])
                                     lx <- xx[!na.ind, 2]; ll <- yy[!na.ind, 2]; lu <- yy[!na.ind, 3]
                                     lpolygon(c(lx, rev(lx)), c(ll, rev(lu)), border = "transparent", 
                                              col = "lightgrey")
                                     panel.xyplot(xx[, 1], yy[, 1], type = "l", lty = 1, col = 2, lwd = 2)
                                     panel.xyplot(xx[, 2], yy[, 2], type = "l", lty = 2, col = 1, lwd = 2)
                                     panel.xyplot(xx[, 3], yy[, 3], type = "l", lty = 2, col = 1, lwd = 2)
                                     panel.xyplot(xx[na.ind, 4], yy[na.ind, 4], type = "p", pch = 8, 
                                                  cex = 1.1, col = 1)
                                     panel.abline(v = last.time, lty = 3)
                                     if (!is.na(target.time))
                                         panel.abline(v = target.time, lty = 2, col = 2, lwd = 2)
                                 }, xlab = "Time", ylim = yl,
                                 ylab = paste("Predicted", resp)))
                }
            }
            dev.off()
        }
    )
})
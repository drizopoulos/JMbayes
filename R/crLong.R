crLong <-
function (data, statusVar, censLevel, 
                    nameStrata = "strata", nameStatus = "status2") {
    n <- nrow(data)
    status <- data[[statusVar]]
    unqLevs <- unique(status)
    unqLevs <- unqLevs[unqLevs != censLevel]
    ncr <- length(unqLevs)
    dataOut <- data[rep(seq_len(n), each = ncr), ]
    dataOut[[nameStrata]] <- rep(unqLevs, n)
    dataOut[[nameStatus]] <- as.numeric(dataOut[[statusVar]] == dataOut[[nameStrata]])
    dataOut[[nameStrata]] <- factor(dataOut[[nameStrata]])
    dataOut
}

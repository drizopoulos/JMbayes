fixNames <-
function (model, names.betas, names.D, names.gammas, names.Bs.gammas, names.alphas, names.Dalphas, 
                      names.shapes, names.id) {
    colnames(model$mcmc$betas) <- names(model$postMeans$betas) <- names(model$postModes$betas) <-
        names(model$StErr$betas) <- names(model$EffectiveSize$betas) <- names(model$StDev$betas) <-
        colnames(model$CIs$betas) <- names.betas
    #######
    if (ss <- !is.null(model$mcmc$sigma)) {
        colnames(model$mcmc$sigma) <- names(model$postMeans$sigma) <- names(model$postModes$sigma) <-
            names(model$StErr$sigma) <- names(model$EffectiveSize$sigma) <- names(model$StDev$sigma) <-
            colnames(model$CIs$sigma) <- "sigma"
    }
    #######
    if (!is.null(model$mcmc[['b']]))
        dimnames(model$mcmc[['b']]) <- list(names.id, names.D)
    dimnames(model$postMeans[['b']]) <- list(names.id, names.D)
    dimnames(model$postVarsRE) <- list(names.D, names.D, names.id)
    #######
    dimnames(model$postMeans$D) <- dimnames(model$postModes$D) <- list(names.D, names.D)
    colnames(model$mcmc$D) <- names(model$StErr$D) <- names(model$EffectiveSize$D) <- names(model$StDev$D) <-
        colnames(model$CIs$D) <- paste0("D[", row(model$postMeans$D), ", ", col(model$postMeans$D), "]")
    #######
    if (!is.null(names.gammas)) {
        colnames(model$mcmc$gammas) <- names(model$postMeans$gammas) <- names(model$postModes$gammas) <-
            names(model$StErr$gammas) <- names(model$EffectiveSize$gammas) <- names(model$StDev$gammas) <-
            colnames(model$CIs$gammas) <- names.gammas        
    }
    #######
    colnames(model$mcmc$Bs.gammas) <- names(model$postMeans$Bs.gammas) <- names(model$postModes$Bs.gammas) <-
        names(model$StErr$Bs.gammas) <- names(model$EffectiveSize$Bs.gammas) <- names(model$StDev$Bs.gammas) <-
        colnames(model$CIs$Bs.gammas) <- names.Bs.gammas
    #######
    if (ss2 <- !is.null(model$mcmc$tauBs)) {
        colnames(model$mcmc$tauBs) <- names(model$postMeans$tauBs) <- names(model$postModes$tauBs) <-
            names(model$StErr$tauBs) <- names(model$EffectiveSize$tauBs) <- names(model$StDev$tauBs) <-
            colnames(model$CIs$tauBs) <- "tauBs"
    }
    #######    
    if (!is.null(names.alphas)) {
        colnames(model$mcmc$alphas) <- names(model$postMeans$alphas) <- names(model$postModes$alphas) <-
            names(model$StErr$alphas) <- names(model$EffectiveSize$alphas) <- names(model$StDev$alphas) <-
            colnames(model$CIs$alphas) <- names.alphas        
    }
    #######    
    if (!is.null(names.Dalphas)) {
        colnames(model$mcmc$Dalphas) <- names(model$postMeans$Dalphas) <- names(model$postModes$Dalphas) <-
            names(model$StErr$Dalphas) <- names(model$EffectiveSize$Dalphas) <- names(model$StDev$Dalphas) <-
            colnames(model$CIs$Dalphas) <- names.Dalphas        
    }
    #######
    if (!is.null(names.shapes)) {
        colnames(model$mcmc$shapes) <- names(model$postMeans$shapes) <- names(model$postModes$shapes) <-
            names(model$StErr$shapes) <- names(model$EffectiveSize$shapes) <- names(model$StDev$shapes) <-
            colnames(model$CIs$shapes) <- names.shapes       
    }
    #######
    ind <- lower.tri(model$postMeans$D, TRUE)
    nn <- c(names.betas, if (ss) "sigma", colnames(model$mcmc$D)[ind], names.gammas, 
            names.Bs.gammas, if (ss2) "tauBs", names.alphas, names.Dalphas, names.shapes)
    dimnames(model$vcov) <- list(nn, nn)
    #######
    model
}

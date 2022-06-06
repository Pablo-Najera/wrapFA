coef_lambda <- function(coef, namesX, namesF){
  nX <- length(namesX)
  nF <- length(namesF)
  lambda <- matrix(NA, nrow = nX, ncol = nF, dimnames = list(namesX, namesF))
  for(f in namesF){
    for(x in namesX){
      lambda[x, f] <- coef$est[grepl(paste0(f, "."), coef$paramHeader) & coef$param == x]
    }
  }
  return(lambda)
}

refine.H <- function(lambda){
  B <- lambda
  C <- t(apply(B, 1, function(j) j / sum(j)))
  v <- apply(C, 2, function(x) mean(x^2))
  s <- apply(C, 2, function(x) sd(x^2))
  t <- v + s/4
  H <- C^2 >= t
  return(H)
}

r2.loads <- function(fit, phi.grid, selector = "BIC", estimator = NULL, mimic.mplus = NULL, mplus.path = NULL, max.iter = 10000, rm.files = FALSE, suppressMessages = TRUE){

  #--------------------
  # 1. Check arguments
  #--------------------

  if(class(fit) != "wrapFA"){stop("fit must be of class 'wrapFA'.")}
  type <- fit$specifications$type
  if(!type %in% c("EFA", "BSEM")){stop("r2.loads can be only applied to EFA or BSEM models.")}
  if(!is.numeric(phi.grid) | any(phi.grid < 0 | phi.grid > 1)){stop("phi.grid must contains numbers between 0 and 1")}
  if(!selector %in% c("ChiSq/df", "CFI", "TLI", "AIC", "BIC", "RMSEA", "SRMR")){stop("selector must be 'ChiSq/df', 'CFI', 'TLI', 'AIC', 'BIC', 'RMSEA', or 'SRMR'.")}
  software <- fit$specifications$software
  if(is.null(mimic.mplus)){
    if(type == "EFA"){
      mimic.mplus <- fit$specifications$mimic.mplus
    } else {
      mimic.mplus <- TRUE
    }
  }
  if(is.null(mplus.path)){
    mplus.path <- "wrapFA_ecfa"
  } else {
    if(!is.character(mplus.path)){stop("mplus.path must be of type character.")}
  }
  categorical <- fit$specifications$categorical
  if(is.null(estimator)){
    if(type == "EFA"){
      estimator <- fit$specifications$estimator
    } else {
      if(is.null(categorical)){
        estimator <- "ML"
      } else {
        estimator <- "WLSMV"
      }
    }
  }
  if(!is.numeric(max.iter)){stop("max.iter must be numeric.")}

  #---------------------------
  # 2. Apply R-squared method
  #---------------------------

  lambda <- fit$lambda
  phi <- fit$phi
  if(software == "lavaan"){
    nobs <- lavaan::lavInspect(fit$out, what = "nobs")
    data <- lavaan::lavInspect(fit$out, what = "data")
  } else if(software == "mplus"){
    nobs <- nrow(fit$out$rdata)
    data <- fit$out$rdata
  }
  J <- nrow(lambda)
  K <- ncol(lambda)
  L <- GDINA::attributepattern(K)[-1,]
  nL <- nrow(L)
  S <- lambda %*% phi

  R2 <- matrix(NA, nrow = J, ncol = (2^K - 1), dimnames = list(paste0("x", 1:J), apply(L, 1, paste, collapse = "")))
  for(l in 1:nL){
    tmpF <- as.numeric(which(L[l,] == 1))
    phi_star <- phi[tmpF, tmpF, drop = FALSE]
    S_star <- S[, tmpF , drop = FALSE]
    R2[,l] <- diag(S_star %*% solve(phi_star) %*% t(S_star))
  }
  PVAF <- t(apply(R2, 1, function(j) j / j[nL]))

  cand.PVAF <- matrix(0, nrow = J, ncol = K + 1, dimnames = list(rownames(PVAF), paste0("K", 0:K)))
  cand.PVAF[,K + 1] <- 1
  cand.q <- matrix(paste(rep("0", K), collapse = ""), nrow = J, ncol = K + 1, dimnames = list(rownames(PVAF), paste0("K", 0:K)))
  cand.q[,K + 1] <- paste(rep("1", K), collapse = "")
  for(j in 1:J){
    for(k in 1:(K - 1)){
      cand.PVAF[j, k + 1] <- max(PVAF[j, rowSums(L) == k])
      cand.q[j, k + 1] <- names(which.max(PVAF[j, rowSums(L) == k]))
    }
  }

  sug.Q <- models <- fits <- list()
  fit.indices <- matrix(NA, nrow = length(phi.grid), ncol = 9, dimnames = list(paste0("phi = ", phi.grid), c("ChiSq", "df", "ChiSq/df", "CFI", "TLI", "AIC", "BIC", "RMSEA", "SRMR")))
  for(p in 1:length(phi.grid)){
    sel.q <- as.character(sapply(1:J, function(j1) cand.q[j1, which((cand.PVAF >= phi.grid[p])[j1,])][1]))
    sel.q.K <- stringr::str_count(sel.q, "1")
    sug.Q.p <- matrix(as.numeric(unlist(strsplit(sel.q, ""))), nrow = J, ncol = K, byrow = TRUE, dimnames = list(rownames(lambda), colnames(lambda)))
    if(p > 1){
      if(any(sapply(1:length(sug.Q), function(l) identical(sug.Q.p, sug.Q[[l]])))){
        sel.p <- which(sapply(1:length(sug.Q), function(l) identical(sug.Q.p, sug.Q[[l]])))[1]
        sug.Q[[p]] <- sug.Q.p; names(sug.Q)[p] <- paste0("phi.grid = ", phi.grid[p])
        models[[p]] <- models[[sel.p]]; names(models)[p] <- paste0("phi.grid = ", phi.grid[p])
        fits[[p]] <- fits[[sel.p]]; names(fits)[p] <- paste0("phi.grid = ", phi.grid[p])
        fit.indices[p,] <- fit.indices[sel.p,]
      } else {
        sug.Q[[p]] <- sug.Q.p; names(sug.Q)[p] <- paste0("phi.grid = ", phi.grid[p])
        models.p <- modFA(sug.Q.p, software = software)$CFA
        models[[p]] <- models.p; names(models)[p] <- paste0("phi.grid = ", phi.grid[p])
        fits.p <- list()
        try({
          fits.p <- CFA(models.p, data, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, mplus.path = paste0(mplus.path, "_", phi.grid[p]), max.iter = max.iter, rm.files = rm.files, suppressMessages = suppressMessages)
        })
        fits[[p]] <- fits.p; names(fits)[p] <- paste0("phi.grid = ", phi.grid[p])
        fit.indices.p <- as.vector(unlist(c(fits.p$model.fit["ChiSq"], fits.p$model.fit["ChiSq_df"], fits.p$model.fit["ChiSq"]/fits.p$model.fit["ChiSq_df"], fits.p$model.fit["CFI"],
                                            fits.p$model.fit["TLI"], fits.p$model.fit["AIC"], fits.p$model.fit["BIC"], fits.p$model.fit["RMSEA"], fits.p$model.fit["SRMR"])))
        if(length(fit.indices.p) > 0){fit.indices[p,] <- fit.indices.p}
      }
    } else {
      sug.Q[[p]] <- sug.Q.p; names(sug.Q)[p] <- paste0("phi.grid = ", phi.grid[p])
      models.p <- modFA(sug.Q.p, software = software)$CFA
      models[[p]] <- models.p; names(models)[p] <- paste0("phi.grid = ", phi.grid[p])
      fits.p <- list()
      try({
        fits.p <- CFA(models.p, data, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, mplus.path = paste0(mplus.path, "_", phi.grid[p]), max.iter = max.iter, rm.files = rm.files, suppressMessages = suppressMessages)
      })
      fits[[p]] <- fits.p; names(fits)[p] <- paste0("phi.grid = ", phi.grid[p])
      fit.indices.p <- as.vector(unlist(c(fits.p$model.fit["ChiSq"], fits.p$model.fit["ChiSq_df"], fits.p$model.fit["ChiSq"]/fits.p$model.fit["ChiSq_df"], fits.p$model.fit["CFI"],
                                          fits.p$model.fit["TLI"], fits.p$model.fit["AIC"], fits.p$model.fit["BIC"], fits.p$model.fit["RMSEA"], fits.p$model.fit["SRMR"])))
      if(length(fit.indices.p) > 0){fit.indices[p,] <- fit.indices.p}
    }
  }

  if(all(apply(fit.indices, 1, function(x) all(is.na(x))))){stop("No convergence for any phi.grid.")}

  sug.phi.grid <- c("ChiSq/df" = phi.grid[as.numeric(which.min(fit.indices[,"ChiSq/df"]))],
                    "CFI" = phi.grid[as.numeric(which.max(fit.indices[,"CFI"]))],
                    "TLI" = phi.grid[as.numeric(which.max(fit.indices[,"TLI"]))],
                    "AIC" = phi.grid[as.numeric(which.min(fit.indices[,"AIC"]))],
                    "BIC" = phi.grid[as.numeric(which.min(fit.indices[,"BIC"]))],
                    "RMSEA" = phi.grid[as.numeric(which.min(fit.indices[,"RMSEA"]))],
                    "SRMR" = phi.grid[as.numeric(which.min(fit.indices[,"SRMR"]))])
  sel.fit <- fits[[paste0("phi.grid = ", as.numeric(sug.phi.grid[selector]))]]

  #-------------------
  # 3. Export results
  #-------------------

  specifications <- list(type = "ECFA", method = "r2", phi.grid = phi.grid, selector = selector, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, mplus.path = mplus.path, max.iter = max.iter)
  res <- list(model.fit = sel.fit$model.fit, lambda = sel.fit$lambda, phi = sel.fit$phi, psi = sel.fit$psi, lambda.p = sel.fit$lambda.p, phi.p = sel.fit$phi.p, psi.p = sel.fit$psi.p, MI = sel.fit$MI, mods = sel.fit$mods, out = sel.fit$out, phi.fit = fit.indices, phi.sug = sug.phi.grid, phi.outs = fits, phi.models = models, phi.Q = sug.Q, specifications = specifications)
  class(res) <- "wrapFA"
  return(res)

}

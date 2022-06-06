#' @include BSEM.R CFA.R ECFA.R EFA.R
#' @export
print.bootFA <- function(x, ...){
  cat(paste0("================================================================", "\n",
             "Congruent coefficient for factor loadings:", "\n", "\n"))
  print(x$lambdaCC)

  cat(paste0("\n", "\n",
             "Congruent coefficient for factor correlations:", "\n", "\n"))
  print(x$phiCC)
  cat("================================================================", "\n")
}

#' @export
print.wrapFA <- function(wrapFA){
  type <- wrapFA$specifications$type
  software <- wrapFA$specifications$software
  if(software == "lavaan"){nobs <- lavaan::lavInspect(wrapFA$out, what = "nobs")}
  if(software == "mplus"){nobs <- nrow(wrapFA$out$rdata)}
  if(software == "mplus"){
    software <- "Mplus (via MplusAutomation)"
  }
  estimator <- wrapFA$specifications$estimator
  npar <- wrapFA$model.fit["npar"]
  nmods <- nrow(wrapFA$mods)
  rotation <- wrapFA$specifications$rotation
  if(is.null(wrapFA$specifications$method)){wrapFA$specifications$method <- "none"}
  if(type %in% c("CFA", "EFA", "ECFA")){
    cat(paste0(
      "================================================================", "\n",
      type, " fitted with ", software, "\n", "\n",
      "Number of observations = ", nobs, "\n",
      "Number of parameters = ", npar, "\n",
      "Estimator = ", estimator, "\n",
      if(!is.null(rotation)){paste0("Rotation procedure = ", rotation, "\n")},
      if(!is.null(nmods)){paste0("Number of modifications introduced in the model = ", nmods, "\n")},
      if(!is.null(wrapFA$phi.sug)){paste0("ECFA method = R-squared; selector = ", wrapFA$specifications$r2.selector, "; phi cutoff = ", as.numeric(wrapFA$phi.sug[wrapFA$specifications$r2.selector]), "\n")},
      if(type == "ECFA" & wrapFA$specifications$method != "r2"){paste0("ECFA method = p-value (", wrapFA$specifications$method, ")", "\n")},
      "\n",
      "CFI = ", wrapFA$model.fit["CFI"], "; TLI = ", wrapFA$model.fit["CFI"], "; SRMR = ", wrapFA$model.fit["SRMR"], "\n",
      "RMSEA = ", wrapFA$model.fit["RMSEA"], "; 90% CI = (", wrapFA$model.fit["RMSEA_low"], ", ", wrapFA$model.fit["RMSEA_up"], "); p-value (< 0.05) = ", wrapFA$model.fit["RMSEA_p"], "\n",
      "AIC = ", wrapFA$model.fit["AIC"], "; BIC = ", wrapFA$model.fit["BIC"], "\n",
      "================================================================"
    ))
  } else if(type == "BSEM"){
    cat(paste0(
      "============================================", "\n",
      type, " fitted with ", software, "\n", "\n",
      "Number of observations = ", nobs, "\n",
      "Number of parameters = ", npar, "\n",
      "Estimator = ", estimator, "\n", "\n",
      "Posterior predictive p-value (PPP) = ", wrapFA$model.fit["PostPredictive_p"], "\n",
      "BIC = ", wrapFA$model.fit["BIC"], "; DIC = ", wrapFA$model.fit["DIC"], "\n",
      "============================================"
    ))
  }
}

#' @export
print.detFA <- function(detFA){
  est.r <- round(as.data.frame(detFA$est.r), 3)
  est.sigma <- round(as.data.frame(detFA$est.sigma), 3)
  if(is.null(detFA$true.r)){
    cat(paste0(
      "==================================================================", "\n",
      "Empirical determinacy of factor score estimates", "\n", "\n",
      "Using the observed correlation matrix:", "\n",
      "Correlation = ", paste(est.r$p, collapse = " "), " (Mean = ", round(mean(est.r$p), 3), ")", "\n",
      "Squared-correlation = ", paste(est.r$p2, collapse = " "), " (Mean = ", round(mean(est.r$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(est.r$minp, collapse = " "), " (Mean = ", round(mean(est.r$minp), 3), ")", "\n", "\n",
      "Using the implied correlation matrix:", "\n",
      "Correlation = ", paste(est.sigma$p, collapse = " "), " (Mean = ", round(mean(est.sigma$p), 3), ")", "\n",
      "Squared-correlation = ", paste(est.sigma$p2, collapse = " "), " (Mean = ", round(mean(est.sigma$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(est.sigma$minp, collapse = " "), " (Mean = ", round(mean(est.sigma$minp), 3), ")", "\n",
      "=================================================================="
    ))
  } else {
    true.r <- round(as.data.frame(detFA$true.r), 3)
    true.sigma <- round(as.data.frame(detFA$true.sigma), 3)
    cat(paste0(
      "==================================================================", "\n",
      "True determinacy of factor score estimates", "\n", "\n",
      "Using the observed correlation matrix:", "\n",
      "Correlation = ", paste(true.r$p, collapse = " "), " (Mean = ", round(mean(true.r$p), 3), ")", "\n",
      "Squared-correlation = ", paste(true.r$p2, collapse = " "), " (Mean = ", round(mean(true.r$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(true.r$minp, collapse = " "), " (Mean = ", round(mean(true.r$minp), 3), ")", "\n", "\n",
      "Using the implied correlation matrix:", "\n",
      "Correlation = ", paste(true.sigma$p, collapse = " "), " (Mean = ", round(mean(true.sigma$p), 3), ")", "\n",
      "Squared-correlation = ", paste(true.sigma$p2, collapse = " "), " (Mean = ", round(mean(true.sigma$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(true.sigma$minp, collapse = " "), " (Mean = ", round(mean(true.sigma$minp), 3), ")", "\n", "\n", "\n",
      "Empirical determinacy of factor score estimates", "\n", "\n",
      "Using the observed correlation matrix:", "\n",
      "Correlation = ", paste(est.r$p, collapse = " "), " (Mean = ", round(mean(est.r$p), 3), ")", "\n",
      "Squared-correlation = ", paste(est.r$p2, collapse = " "), " (Mean = ", round(mean(est.r$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(est.r$minp, collapse = " "), " (Mean = ", round(mean(est.r$minp), 3), ")", "\n", "\n",
      "Using the implied correlation matrix:", "\n",
      "Correlation = ", paste(est.sigma$p, collapse = " "), " (Mean = ", round(mean(est.sigma$p), 3), ")", "\n",
      "Squared-correlation = ", paste(est.sigma$p2, collapse = " "), " (Mean = ", round(mean(est.sigma$p2), 3), ")", "\n",
      "Minimum correlation = ", paste(est.sigma$minp, collapse = " "), " (Mean = ", round(mean(est.sigma$minp), 3), ")", "\n",
      "=================================================================="
    ))
  }
}

#' @export
print.modFA <- function(modFA){
  software <- modFA$specifications$software
  switch(software,
         "lavaan" = cat(paste0("===================================", "\n",
                               "CFA and EFA models in lavaan syntax", "\n", "\n",
                               "CFA model:", "\n",
                               modFA$CFA, "\n", "\n",
                               "EFA model (mechanical rotation):", "\n",
                               modFA$EFA, "\n",
                               "===================================")),
         "mplus" = cat(paste0("=========================================", "\n",
                              "CFA, EFA, and BSEM models in Mplus syntax", "\n", "\n",
                              "CFA model:", "\n",
                              modFA$CFA, "\n", "\n",
                              "EFA model (mechanical rotation):", "\n",
                              modFA$EFA, "\n", "\n",
                              "EFA model (target rotation):", "\n",
                              modFA$EFAt, "\n", "\n",
                              "BSEM model:", "\n",
                              modFA$BSEM$model, "\n",
                              modFA$BSEM$priors, "\n",
                              "========================================="))
  )
}

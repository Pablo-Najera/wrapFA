#' Nonparametric bootstrap to evaluate the stability of the internal structure
#'
#' @description Evaluates the stability of a factor analysis solution (e.g., CFA, EFA, ECFA, BSEM) by conducting nonparametric
#' bootstrapping, resampling with a replacement from the original dataset. This is used in Nájera et al. (2022), based on the
#' bootstrap EGA (Christensen & Golino, 2021).
#'
#' @param fit A \code{wrapFA} object or a \code{list} of \code{wrapFA} objects.
#' @param R A \code{double} indicating the number of bootstrapping replications. Default is 100.
#' @param n.cores A \code{double} indicating the number of CPU processors to be used for the bootstrapping. Default is 2.
#' @param plots A \code{logical} indicating whether result plots should be generated. Default is \code{TRUE}.
#' @param rm.files A \code{logical} indicating whether the stored files should be deleted. The default is \code{TRUE}.
#' @param seed A \code{double} containing a random seed. Default is \code{NULL}.
#' @param verbose A \code{logical} indicating whether a progress bar should be printed. Default is \code{TRUE}.
#' @param digits A \code{double} indicating the number of decimal digits in the outcome. Default is 3.
#'
#' @return \code{bootFA} returns an object of class \code{bootFA}.
#' \describe{
#' \item{\code{lambdaCC}}{Average congruent coefficient of factor loadings between pairs of bootraspping solutions (\code{matrix}).}
#' \item{\code{lambdaMAD}}{Average mean absolute deviation of factor loadings  between pairs of bootraspping solutions (\code{matrix}).}
#' \item{\code{phiCC}}{Average congruent coefficient of factor correlations between pairs of bootraspping solutions (\code{matrix}).}
#' \item{\code{phiMAD}}{Average mean absolute deviation of factor correlations between pairs of bootraspping solutions (\code{matrix}).}
#' \item{\code{lambdaCC.r}}{Congruent coefficient of factor loadings between pairs of bootraspping solutions (\code{list}).}
#' \item{\code{lambdaMAD.r}}{Mean absolute deviation of factor loadings  between pairs of bootraspping solutions (\code{list}).}
#' \item{\code{phiCC.r}}{Congruent coefficient of factor correlations between pairs of bootraspping solutions (\code{list}).}
#' \item{\code{phiMAD.r}}{Mean absolute deviation of factor correlations between pairs of bootraspping solutions (\code{list}).}
#' \item{\code{plots}}{Plots summarizing the results. Only printed if \code{plots = TRUE} (\code{list}).}
#' \item{\code{specifications}}{Specifications used to run the function (\code{list}).}
#' }
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021). Estimating the stability of the number of factors via bootstrap exploratory graph analysis: A tutorial. \emph{Psych}, \emph{3}(3), 479-500. https://doi.org/10.3390/psych3030032
#' Nájera, P., Abad, F. J., & Sorrel, M. A. (2022). Is EFA always to be preferred? A systematic comparison of factor analytic techniques throughout the confirmatory-exploratory continuum. \emph{Manuscript submitted for publication}.
#'
#' @import foreach
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(MBESS)
#' data(HS)
#' HS <- HS[, -c(1:8)]
#' colnames(HS) <- paste0("t", 1:26)
#' modHS <- data.frame(spatial = c(rep(1, 4), rep(0, 20), rep(1, 2)),
#'                     verbal = c(rep(0, 4), rep(1, 5), rep(0, 17)),
#'                     speed = c(rep(0, 9), rep(1, 4), rep(0, 13)),
#'                     memory = c(rep(0, 13), rep(1, 6), rep(0, 7)),
#'                     maths = c(rep(0, 19), rep(1, 5), rep(0, 2)), row.names = colnames(HS))
#' modHS.lav <- modFA(lambda = modHS, software = "lavaan", keep.names = TRUE)
#' cfa.lav <- CFA(model = modHS.lav$CFA, data = HS)
#' boot.cfa <- bootFA(fit = cfa.lav, n.cores = 4)
#' }
bootFA <- function(fit, R = 100, n.cores = 2, plots = TRUE, rm.files = TRUE, seed = NULL, verbose = TRUE, digits = 3){

  if(!is.null(seed)){set.seed(seed)}
  check.class <- c()

  if(!inherits(fit, "list")){fit <- list(fit)}
  for(f in fit){
    check.class <- c(inherits(f, "wrapFA"))
  }
  if(!all(check.class)){stop("All objects in fit must be of class wrapFA.")}

  cl <- parallel::makeCluster(n.cores, type = "SOCK")
  doSNOW::registerDoSNOW(cl)

  if(verbose){
    cat("Bootstrapping Progress:", "\n")
    pb <- utils::txtProgressBar(max = R, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }

  fit.R <- foreach::foreach(r = 1:R,
                            .options.snow = opts,
                            .packages = c("lavaan", "MplusAutomation", "fungible", "psych", "MBESS"),
                            .export = c("CFA", "EFA", "ECFA", "BSEM"),
                            .combine = list, .inorder = TRUE) %dopar% {
                              try({
                                out.r <- list()
                                for(f in 1:length(fit)){
                                  s <- fit[[f]]$specifications
                                  mplus.path <- paste0(s$mplus.path, "_R", r)
                                  data.r <- s$data[sample(c(1:nrow(s$data)), size = nrow(s$data), replace = TRUE),]
                                  fit.r <- switch(s$type,
                                                  CFA = CFA(s$model, data.r, s$categorical, s$estimator, s$software, s$mimic.mplus, mplus.path, rm.files = rm.files, s$max.iter, s$mod, s$mod.args),
                                                  EFA = EFA(s$model, data.r, s$n.factors, s$categorical, s$estimator, s$rotation, s$software, s$mimic.mplus, mplus.path,  rm.files = rm.files, s$max.iter, s$mod, s$mod.args),
                                                  ECFA = ECFA(s$fit, s$model, data.r, s$method, s$r2.args, s$categorical, s$estimator, s$software, s$mimic.mplus, mplus.path, rm.files = rm.files, s$max.iter),
                                                  BSEM = BSEM(s$model, data.r, s$categorical, mplus.path, rm.files = rm.files, s$max.iter, s$n.cores, s$n.chains)
                                  )
                                  if(rm.files & s$software == "mplus"){file.remove(mplus.path)}
                                  fit.r$boot.r <- r
                                  out.r[[f]] <- fit.r
                                }
                                return(out.r)
                              })
                            }
  parallel::stopCluster(cl)

  back.fit.R <- fit.R
  class.f <- c()
  for(f in fit.R){class.f <- c(class.f, class(f))}
  out <- fit.R[which(class.f == "wrapFA")]
  if(length(which(class.f == "wrapFA")) > 0){
    fit.R <- fit.R[-which(class.f == "wrapFA")]
  }
  continue <- ifelse(any(class.f == "list"), TRUE, FALSE)
  while(continue){
    fit.R <- unlist(fit.R, recursive = FALSE)
    class.f <- c()
    for(f in fit.R){class.f <- c(class.f, class(f))}
    out <- append(out, fit.R[which(class.f == "wrapFA")])
    if(length(which(class.f == "wrapFA")) > 0){
      fit.R <- fit.R[-which(class.f == "wrapFA")]
    }
    continue <- ifelse(any(class.f == "list"), TRUE, FALSE)
  }

  config <- data.frame(type = sapply(out, function(x) x$specifications$type), r = sapply(out, function(x) x$boot.r))
  out <- out[order(config$type, config$r)]
  config <- config[order(config$type, config$r),]
  conv <- 100 * table(config$type) / R

  lambdaCC.r <- lambdaMAD.r <- phiCC.r <- phiMAD.r <- list()
  for(t in unique(config$type)){
    base <- fit[[which(sapply(fit, function(x) x$specifications$type) == t)]]
    lambdaCC.r[[t]] <- lambdaMAD.r[[t]] <- phiCC.r[[t]] <- phiMAD.r[[t]] <- matrix(NA, nrow = 0, ncol = ncol(base$lambda), dimnames = list(NULL, colnames(base$lambda)))
    for(r1 in 1:(max(config$r) - 1)){
      fit1 <- out[[which(config$type == t & config$r == r1)]]
      align1 <- try(fungible::faAlign(F1 = base$lambda, F2 = fit1$lambda, Phi2 = fit1$phi))
      for(r2 in (r1 + 1):max(config$r)){
        fit2 <- out[[which(config$type == t & config$r == r2)]]
        align2 <- try(fungible::faAlign(F1 = base$lambda, F2 = fit2$lambda, Phi2 = fit2$phi))
        lambdaCC.r[[t]] <- rbind(lambdaCC.r[[t]], diag(psych::factor.congruence(align1$F2, align2$F2)))
        lambdaMAD.r[[t]] <- rbind(lambdaMAD.r[[t]], colMeans(abs(align1$F2 - align2$F2)))
        phiCC.r[[t]] <- rbind(phiCC.r[[t]], diag(psych::factor.congruence(align1$Phi2, align2$Phi2)))
        phiMAD.r[[t]] <- rbind(phiMAD.r[[t]], colMeans(abs(align1$Phi2 - align2$Phi2)))
        rownames(lambdaCC.r[[t]])[nrow(lambdaCC.r[[t]])] <- paste0(r1, "-", r2)
        rownames(lambdaMAD.r[[t]])[nrow(lambdaMAD.r[[t]])] <- paste0(r1, "-", r2)
        rownames(phiCC.r[[t]])[nrow(phiCC.r[[t]])] <- paste0(r1, "-", r2)
        rownames(phiMAD.r[[t]])[nrow(phiMAD.r[[t]])] <- paste0(r1, "-", r2)
      }
    }
  }

  lambdaCC <- round(t(sapply(lambdaCC.r, colMeans)), digits)
  lambdaMAD <- round(t(sapply(lambdaMAD.r, colMeans)), digits)
  phiCC <- round(t(sapply(phiCC.r, colMeans)), digits)
  phiMAD <- round(t(sapply(phiMAD.r, colMeans)), digits)

  if(plots){
    out.plot <- list()
    for(t in unique(config$type)){
      lambdaCC.p <- data.frame(Dimension = rep(colnames(lambdaCC.r[[t]]), each = nrow(lambdaCC.r[[t]])), CC = as.numeric(lambdaCC.r[[t]]))
      lambdaCC.p$Dimension <- factor(lambdaCC.p$Dimension, levels = unique(lambdaCC.p$Dimension))
      lambdaCC.p <- ggplot2::ggplot(lambdaCC.p, ggplot2::aes(x = Dimension, y = CC)) +
        ggplot2::geom_violin(color = "transparent", fill = "grey70") +
        # ggplot2::geom_jitter(pch = 21, fill = "white", height = 0, width = 0.2) +
        ggplot2::stat_summary(geom = "point", fun = mean, pch = 4, size = 3) +
        ggplot2::scale_y_continuous("Loadings CC", limits = c(0, 1)) +
        ggplot2::theme_bw(base_size = 12)
      lambdaMAD.p <- data.frame(Dimension = rep(colnames(lambdaMAD.r[[t]]), each = nrow(lambdaMAD.r[[t]])), MAD = as.numeric(lambdaMAD.r[[t]]))
      lambdaMAD.p$Dimension <- factor(lambdaMAD.p$Dimension, levels = unique(lambdaMAD.p$Dimension))
      lambdaMAD.p <- ggplot2::ggplot(lambdaMAD.p, ggplot2::aes(x = Dimension, y = MAD)) +
        ggplot2::geom_violin(color = "transparent", fill = "grey70") +
        # ggplot2::geom_jitter(pch = 21, fill = "white", height = 0, width = 0.2) +
        ggplot2::stat_summary(geom = "point", fun = mean, pch = 4, size = 3) +
        ggplot2::scale_y_continuous("Loadings MAD", limits = c(0, max(lambdaMAD.p$MAD))) +
        ggplot2::theme_bw(base_size = 12)
      phiCC.p <- data.frame(Dimension = rep(colnames(phiCC.r[[t]]), each = nrow(phiCC.r[[t]])), CC = as.numeric(phiCC.r[[t]]))
      phiCC.p$Dimension <- factor(phiCC.p$Dimension, levels = unique(phiCC.p$Dimension))
      phiCC.p <- ggplot2::ggplot(phiCC.p, ggplot2::aes(x = Dimension, y = CC)) +
        ggplot2::geom_violin(color = "transparent", fill = "grey70") +
        # ggplot2::geom_jitter(pch = 21, fill = "white", height = 0, width = 0.2) +
        ggplot2::stat_summary(geom = "point", fun = mean, pch = 4, size = 3) +
        ggplot2::scale_y_continuous("Correlations CC", limits = c(0, 1)) +
        ggplot2::theme_bw(base_size = 12)
      phiMAD.p <- data.frame(Dimension = rep(colnames(phiMAD.r[[t]]), each = nrow(phiMAD.r[[t]])), MAD = as.numeric(phiMAD.r[[t]]))
      phiMAD.p$Dimension <- factor(phiMAD.p$Dimension, levels = unique(phiMAD.p$Dimension))
      phiMAD.p <- ggplot2::ggplot(phiMAD.p, ggplot2::aes(x = Dimension, y = MAD)) +
        ggplot2::geom_violin(color = "transparent", fill = "grey70") +
        # ggplot2::geom_jitter(pch = 21, fill = "white", height = 0, width = 0.2) +
        ggplot2::stat_summary(geom = "point", fun = mean, pch = 4, size = 3) +
        ggplot2::scale_y_continuous("Correlations MAD", limits = c(0, max(phiMAD.p$MAD))) +
        ggplot2::theme_bw(base_size = 12)
      out.plot[[t]] <- list(lambdaCC = lambdaCC.p, lambdaMAD = lambdaMAD.p, phiCC = phiCC.p, phiMAD = phiMAD.p)
    }
  }

  specs <- list(fit = fit, R = R, n.cores = n.cores, plots = plots, rm.files = rm.files, seed = seed, verbose = verbose, digits = digits)
  res <- list(lambdaCC = lambdaCC, lambdaMAD = lambdaMAD, phiCC = phiCC, phiMAD = phiMAD,
              lambdaCC.r = lambdaCC.r, lambdaMAD.r = lambdaMAD.r, phiCC.r = phiCC.r, phiMAD.r = phiMAD.r,
              plots = out.plot, specifications = specs)
  class(res) <- "bootFA"
  return(res)
}

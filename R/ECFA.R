#' EFA-based confirmatory factor analysis (ECFA)
#'
#' @description Fit an ECFA model (Nájera et al., 2022) via \emph{lavaan} or \emph{MplusAutomation}.
#'
#' @param fit A \code{wrapFA} object from \code{EFA()} or \code{BSEM()} functions. Default is \code{NULL}, which means that a user-specified model will be used.
#' @param model A \code{character} string containing the EFA model using either the \emph{Mplus} or the \emph{lavaan} syntax. Default is \code{NULL}, which means that the argument \code{fit} will be used.
#' @param data A \emph{N} individuals x \emph{J} observed variables \code{matrix} or \code{data.frame}. Missing values need to be coded as \code{NA}.
#' @param method A \code{character} indicating the method to detect the relevant loadings. It can be \code{"r2"} for the R-squared procedure (Citation, Year), \code{"nominal"} for statistically significant loadings with \emph{alpha} = 0.05, or \code{"bonferroni"} for statistically significant loadings using the Bonferroni correction. Default is \code{"r2"}.
#' @param r2.args A \code{list} of arguments about the R-squared procedure:
#' \describe{
#' \item{\code{phi.grid}}{A \code{numeric} vector indicating the phi parameter space to search. Default is from 0.70 to 0.99, in steps of 0.01.}
#' \item{\code{selector}}{A \code{character} indicating fit index to use to select the final model. It can be \code{"ChiSq/df"}, \code{"CFI"}, \code{"TLI"}, \code{"AIC"}, \code{"BIC"}, \code{"aBIC"}, \code{"RMSEA"}, \code{"SRMR"}, and \code{"AICC"}. Default is \code{"BIC"}.}
#' }
#' @param categorical A \code{vector} containing the number or name of the observed variables that need to be treated as categorical. Default is \code{NULL}, which treates all observed variables as continuous.
#' @param estimator A \code{character} indicating the estimator to be used to fit the model. It can be \code{"ML"}, \code{"MLM"}, \code{"MLMV"}, \code{"MLMVS"}, \code{"MLF"}, \code{"MLR"}, \code{"GLS"}, \code{"DWLS"}, \code{"WLS"}, \code{"WLSM"}, \code{"WLSMV"}, \code{"ULS"},\code{"ULSM"}, and \code{"ULSMV"}. Check \code{?lavaan::lavOptions} or the Mplus User's Guide for more information. Default is \code{NULL}, which uses \code{"ML"} if all observed variables are continuous and \code{"WLSMV"} if any observed variable is categorical.
#' @param software A \code{character} indicating the software to use for fitting the CFA model. It can be \code{"lavaan"} or \code{"mplus"}. Only required if design matrices are used in the \code{model} argument; the software is automatically detected if a syntax is used. Default is \code{NULL}.
#' @param mimic.mplus A \code{logical} indicating whether \emph{lavaan} should mimic \emph{Mplus} output. Only applicable to \emph{lavaan}. Default is \code{TRUE}.
#' @param mplus.path A \code{character} indicating the path where the the files generated by Mplus will be stored, including the name of the files. Default is \code{NULL}, which stores the files in the current working directory with the name \emph{wrapECFA}.
#' @param rm.files A \code{logical} indicating whether the stored files should be deleted. The default is \code{FALSE}.
#' @param max.iter A \code{double} indicating the maximum number of iterations. Default is 10000.
#' @param suppressMessages A \code{logical} indicating whether the messages from MplusAutomation should be silent. Default is \code{TRUE}.
#'
#' @return \code{ECFA} returns an object of class \code{wrapFA}.
#' \describe{
#' \item{\code{model.fit}}{Model fit indices (\code{vector}).}
#' \item{\code{lambda}}{The standardized estimated factor loading matrix (\code{matrix}).}
#' \item{\code{phi}}{The standardized estimated factor correlation matrix (\code{matrix}).}
#' \item{\code{psi}}{The standardized estimated error correlation matrix (\code{matrix}).}
#' \item{\code{lambda.p}}{The p-value associated to each free parameter in the estimated factor loading matrix (\code{matrix}).}
#' \item{\code{phi.p}}{The p-value associated to each free parameter in the estimated factor correlation matrix (\code{matrix}).}
#' \item{\code{psi.p}}{The p-value associated to each free parameter in the estimated error correlation matrix (\code{matrix}).}
#' \item{\code{MI}}{Information about the MI and SEPC of the fixed parameters (\code{matrix}).}
#' \item{\code{out}}{The complete output provided by either \emph{lavaan} or \emph{Mplus} (\code{list}).}
#' \item{\code{specifications}}{Specifications used to run the function (\code{list}).}
#' }
#'
#' @references
#' Nájera, P., Abad, F. J., & Sorrel, M. A. (2022). Is EFA always to be preferred? A systematic comparison of factor analytic techniques throughout the confirmatory-exploratory continuum. \emph{Manuscript submitted for publication}.
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
#' efa_geomin <- EFA(model = modHS.lav$EFA, data = HS)
#' ecfa_r2 <- ECFA(fit = efa_geomin, data = HS, method = "r2", r2.args = list(phi.grid = seq(0.85, 0.95, 0.01), selector = "BIC"))
#' }
ECFA <- function(fit = NULL, model = NULL, data, method = "r2", r2.args = list(phi.grid = seq(0.70, 0.99, 0.01), selector = "BIC"), categorical = NULL, estimator = NULL, software = NULL, mimic.mplus = TRUE, mplus.path = NULL, rm.files = FALSE, max.iter = 10000, suppressMessages = TRUE){

  #--------------------
  # 1. Check arguments
  #--------------------

  if(is.null(fit) & is.null(model)){stop("Either fit or model should be specified.")}
  if(!is.null(fit) & !is.null(model)){
    warning("model will be ignored because fit has been specified.")
    model <- NULL
  }
  if(is.null(model) & is.null(software)){software <- fit$specifications$software}
  if(!is.null(fit)){if(class(fit) != "wrapFA"){stop("fit only supports wrapFA objects.")}}
  if(!is.null(fit)){
    type <- names(fit)[1]
    if(is.null(software)){software <- fit$software}
  }
  if(!is.null(model)){
    if(is.character(model)){
      if(grepl("BY", model)){
        software <- "mplus"
        if(any(grepl("(\\*1)", model))){
          if(grepl("~0", model)){
            type <- "EFAt"
          } else {
            type <- "EFAg"
          }
        } else if(length(model) == 2){
          type <- "BSEM"
        } else {
          type <- "CFA"
        }
      } else if(grepl("=~", model)){
        software <- "lavaan"
        if(grepl("block", model)){
          type <- "EFA"
        } else {
          type <- "CFA"
        }
      } else {
        stop("model must be specified either by an Mplus or lavaan syntax.")
      }
    }
  }
  if(type == "CFA"){stop("CFA model has been identified. ECFA() works either with EFA or BSEM models.")}
  if(is.matrix(data)){
    data <- as.data.frame(data)
    if(all(colnames(data) == paste0("V", 1:ncol(data)))){colnames(data) <- paste0("x", 1:ncol(data))}
  } else if(!is.data.frame(data)){
    stop("data must be either a matrix or data.frame.")
  }
  if(!method %in% c("r2", "nominal", "bonferroni")){stop("method must be 'r2', 'nominal', or 'bonferroni'.")}
  if(method == "r2"){
    if(is.null(r2.args$phi.grid)){r2.args$phi.grid <- seq(0.70, 0.99, 0.01)}
    if(is.null(r2.args$selector)){
      if(is.null(categorical)){
        r2.args$selector <- "BIC"
      } else {
        r2.args$selectr <- "RMSEA"
      }
    } else {
      if(is.null(categorical)){
        if(!r2.args$selector %in% c("ChiSq/df", "CFI", "TLI", "AIC", "BIC", "RMSEA", "SRMR")){stop("r2.args$selector must be 'ChiSq/df', 'CFI', 'TLI', 'AIC', 'BIC', 'RMSEA', or 'SRMR'.")}
      } else {
        if(!r2.args$selector %in% c("ChiSq/df", "CFI", "TLI", "RMSEA")){stop("r2.args$selector must be 'ChiSq/df', 'CFI', 'TLI', or 'RMSEA' if categorical variables are involved in the analysis.")}
      }
    }
  }
  if(is.null(categorical)){
    if(is.null(estimator)){estimator <- "ML"}
  } else {
    if(is.null(estimator)){estimator <- "WLSMV"}
    if(is.numeric(categorical)){
      categorical <- colnames(data)[categorical]
    } else if (is.character(categorical)){
      if(!all(categorical %in% colnames(data))){stop("categorical contains some observed variables that are not detected in data.")}
    } else {
      stop("categorical must be either a numerical or character vector.")
    }
  }
  if(!estimator %in% c("ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR", "GLS", "DWLS", "WLS", "WLSM", "WLSMV", "ULS","ULSM", "ULSMV")){stop("estimator is not recognized.")}
  if(!is.numeric(max.iter)){stop("max.iter must be numeric.")}
  if(is.null(mplus.path)){
    mplus.path <- "wrapECFA"
  } else {
    if(!is.character(mplus.path)){stop("mplus.path must be of type character.")}
  }
  if(mimic.mplus){
    mimic <- "Mplus"
  } else {
    mimic <- "lavaan"
  }

  #----------------
  # 2. lavaan ECFA
  #----------------

  if(software == "lavaan"){

    # 2.1. Extract or fit model
    if(!is.null(fit)){
      FIT <- fit
    } else {
      FIT <- EFA(model, data, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, max.iter = max.iter)
    }

    # 2.2. Identify relevant loadings
    if(method == "r2"){
      ECFA <- r2.loads(fit = FIT, phi.grid = r2.args$phi.grid, selector = r2.args$selector, estimator = estimator, mimic.mplus = mimic.mplus, mplus.path = mplus.path, max.iter = max.iter, rm.files = rm.files, suppressMessages = suppressMessages)
    } else if(method %in% c("nominal", "bonferroni")){
      if(method == "nominal"){
        alfa <- 0.05
      } else {
        alfa <- 0.05 / (nrow(FIT$lambda) * ncol(FIT$lambda))
      }
      p.lambda <- matrix(as.numeric(FIT$lambda.p < alfa), ncol = ncol(FIT$lambda))
      rownames(p.lambda) <- rownames(FIT$lambda.p)
      colnames(p.lambda) <- colnames(FIT$lambda.p)
      for(nl in which(rowSums(p.lambda) == 0)){p.lambda[nl, which.min(FIT$lambda.p[nl,])] <- 1}
      model <- modFA(p.lambda, software = software)
      ECFA <- CFA(model$CFA, data, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, max.iter = max.iter, rm.files = rm.files, suppressMessages = suppressMessages)
    }

  }

  #---------------
  # 3. mplus ECFA
  #---------------

  if(software == "mplus"){

    # 3.1. Extract or fit model
    if(!is.null(fit)){
      FIT <- fit
    } else {
      if(type == "EFAg"){FIT <- EFA(model, data, categorical = categorical, estimator = estimator, max.iter = max.iter, software = software, mimic.mplus = mimic.mplus)}
      if(type == "EFAt"){FIT <- EFA(model, data, categorical = categorical, estimator = estimator, rotation = "TARGET", max.iter = max.iter, software = software, mimic.mplus = mimic.mplus)}
      if(type == "BSEM"){FIT <- BSEM(model, data, categorical = categorical, max.iter = max.iter)}
    }

    # 3.2. Identify relevant loadings
    if(method == "r2"){
      ECFA <- r2.loads(fit = FIT, phi.grid = r2.args$phi.grid, selector = r2.args$selector, estimator = estimator, mimic.mplus = mimic.mplus, mplus.path = mplus.path, max.iter = max.iter, rm.files = rm.files, suppressMessages = suppressMessages)
    } else if(method %in% c("nominal", "bonferroni")){
      if(method == "nominal"){
        alfa <- 0.05
      } else {
        alfa <- 0.05 / (nrow(FIT$lambda) * ncol(FIT$lambda))
      }
      p.lambda <- matrix(as.numeric(FIT$lambda.p < alfa), ncol = ncol(FIT$lambda))
      rownames(p.lambda) <- rownames(FIT$lambda.p)
      colnames(p.lambda) <- colnames(FIT$lambda.p)
      for(nl in which(rowSums(p.lambda) == 0)){p.lambda[nl, which.min(FIT$lambda.p[nl,])] <- 1}
      model <- modFA(p.lambda, software = software)
      ECFA <- CFA(model$CFA, data, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, max.iter = max.iter, mplus.path = mplus.path, rm.files = rm.files, suppressMessages = suppressMessages)
    }

  }

  #-------------------
  # 4. Export results
  #-------------------

  if(rm.files & software == "mplus"){file.remove(paste0(mplus.path, c(".dat", ".inp", ".out")))}
  specifications <- list(type = "ECFA", fit = fit, model = model, data = data, method = method, r2.args = r2.args, categorical = categorical, estimator = estimator, software = software, mimic.mplus = mimic.mplus, mplus.path = mplus.path, rm.files = rm.files, max.iter = max.iter)
  res <- ECFA
  res$specifications <- specifications
  class(res) <- "wrapFA"
  return(res)

}

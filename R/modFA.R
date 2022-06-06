#' Translate design matrices into syntax for lavaan or Mplus models
#'
#' @description Translates design matrices (factor loadings, factor correlations, residual correlations) into a syntax that is compatible with either
#' \emph{lavaan} and \emph{Mplus}. Namely, CFA and EFA are supported for \emph{lavaan}, and CFA, EFA, EFA with target rotation, and BSEM are supported for \emph{Mplus}.
#'
#' @param lambda Factor loading \code{matrix} (or \code{data.frame}).
#' @param phi Factor correlation \code{matrix} (or \code{data.frame}). Default is \code{NULL}, which implies that all factor correlations are included in the model.
#' @param psi Residual correlation \code{matrix} (or \code{data.frame}). Default is \code{NULL}, which implies that no residual correlations are included in the model.
#' @param software A \code{character} indicating the software to use for translating the design matrices into syntax. It can be \code{"lavaan"}, \code{"mplus"}, or \code{"lslx"}.
#' @param keep.names A \code{logical} indicating whether the names of the observed and latent variables (taken from \code{lambda}) should be retained or not. Default is \code{TRUE}. If \code{FALSE}, \code{x1, x2, ...} is used for observed variables and \code{F1, F2, ...} for latent variables.
#'
#' @return \code{modFA} returns an object of class \code{modFA}.
#' \describe{
#' \item{\code{CFA}}{The syntax for fitting a CFA (\code{character}).}
#' \item{\code{EFA}}{The syntax for fitting an EFA with mechanical/analytical rotation (\code{character}).}
#' \item{\code{EFAt}}{The syntax for fitting an EFA with target rotation. Only available for \emph{Mplus} (\code{character}).}
#' \item{\code{BSEM}}{The syntax for fitting a BSEM. Only available for \emph{Mplus} (\code{character}).}
#' \item{\code{specifications}}{Specifications used to implement the \code{modFA} function (\code{list}).}
#' }
#'
#' @export
#'
#' @examples
#' bfi <- psych::bfi
#' bfi <- bfi[,1:25]
#' modelbfi <- data.frame(agree = c(rep(1, 5), rep(0, 20)),
#'                        consc = c(rep(0, 5), rep(1, 5), rep(0, 15)),
#'                        extra = c(rep(0, 10), rep(1, 5), rep(0, 10)),
#'                        neuro = c(rep(0, 15), rep(1, 5), rep(0, 5)),
#'                        openn = c(rep(0, 20), rep(1, 5)), row.names = colnames(bfi))
#' modelbfi.lav <- modFA(lambda = modelbfi, software = "lavaan", keep.names = TRUE)
#' modelbfi.mpl <- modFA(lambda = modelbfi, software = "mplus", keep.names = TRUE)
modFA <- function(lambda, phi = NULL, psi = NULL, software, keep.names = TRUE){

  #--------------------
  # 1. Check arguments
  #--------------------

  if(!is.matrix(lambda) & !is.data.frame(lambda)){stop("lambda must be a matrix or data.frame.")}
  if(!is.null(phi) & !is.matrix(phi) & !is.data.frame(phi)){stop("phi must be a matrix or data.frame.")}
  if(!is.null(psi) & !is.matrix(psi) & !is.data.frame(psi)){stop("psi must be a matrix or data.frame.")}
  if(!software %in% c("mplus", "lavaan", "lslx")){stop("sofware must be 'mplus', 'lavaan', or 'lslx'.")}
  if(!is.logical(keep.names)){stop("keep.names must be logical.")}

  #-----------------------
  # 2. Gather information
  #-----------------------

  J <- nrow(lambda)
  K <- ncol(lambda)
  if(is.null(phi)){phi <- matrix(1, K, K)}
  if(is.null(psi)){psi <- diag(J)}
  if(!keep.names | is.null(rownames(lambda)) | is.null(colnames(lambda))){
    colnames(lambda) <- rownames(phi) <- colnames(phi) <- paste0("F", 1:K)
    rownames(lambda) <- rownames(psi) <- colnames(psi) <- paste0("x", 1:J)
  } else {
    colnames(phi) <- rownames(phi) <- colnames(lambda)
    colnames(psi) <- rownames(psi) <- rownames(lambda)
  }
  loads <- apply(lambda, 2, function(j) names(which(j == 1)))
  if(is.matrix(loads)){
    tmp <- list()
    for(k in 1:K){
      tmp[[k]] <- loads[,k]
    }
    names(tmp) <- colnames(loads)
    loads <- tmp
  }
  cors <- unique(t(apply(which(phi == 0, arr.ind = TRUE), 1, sort)))
  resi <- unique(t(apply(which(psi == 1, arr.ind = TRUE), 1, sort)))
  resi <- resi[-which(resi[,1] == resi[,2]),,drop = FALSE]
  specifications <- list(software = software, obs.names = rownames(lambda), lat.names = colnames(lambda))

  #-------------------------
  # 3. Create lavaan models
  #-------------------------

  if(software == "lavaan"){
    CFA.lambda <- EFA.lambda <- c()
    for(k in 1:K){
      CFA.lambda <- c(CFA.lambda, paste(names(loads)[[k]], "=~", paste(loads[[k]], collapse = " + ")))
      EFA.lambda <- c(EFA.lambda, paste0("efa('block1')*", names(loads)[[k]], " =~ ", paste(rownames(lambda), collapse = " + ")))
    }

    CFA.lambda <- paste(CFA.lambda, collapse = "\n")
    EFA.lambda <- paste(EFA.lambda, collapse = "\n")

    m.cors <- paste(as.vector(apply(cors, 1, function(p) paste0(colnames(phi)[p[1]], " ~~ 0*", colnames(phi)[p[2]]))), collapse = "\n")
    m.resi <- paste(as.vector(apply(resi, 1, function(p) paste0(colnames(psi)[p[1]], " ~~ ", colnames(psi)[p[2]]))), collapse = "\n")

    if(nrow(cors) > 0){
      CFA.lavaan <- paste(CFA.lambda, m.cors, sep = "\n")
      EFA.lavaan <- paste(EFA.lambda, m.cors, sep = "\n")
    }
    if(nrow(resi) > 0){
      CFA.lavaan <- paste(CFA.lambda, m.resi, sep = "\n")
      EFA.lavaan <- paste(EFA.lambda, m.resi, sep = "\n")
    }
    if(nrow(cors) == 0 & nrow(resi) == 0){
      CFA.lavaan <- CFA.lambda
      EFA.lavaan <- EFA.lambda
    }

    res <- list(CFA = CFA.lavaan, EFA = EFA.lavaan, specifications = specifications)
    class(res) <- "modFA"
    return(res)
  }

  #------------------------
  # 4. Create mplus models
  #------------------------

  if(software == "mplus"){
    CFA.lambda <- EFAt.lambda <- BSEM.lambda <- BSEM.priors <- c()
    for(k in 1:K){
      CFA.lambda <- c(CFA.lambda, paste(names(loads)[[k]], "BY", paste(loads[[k]], collapse = " ")))
      EFAt.lambda <- c(EFAt.lambda, paste(names(loads)[[k]], "BY", paste(loads[[k]], collapse = " ")))
      EFAt.lambda[k] <- paste(paste(EFAt.lambda[k], paste(paste0(rownames(lambda)[as.numeric(which(lambda[,k] != 1))], "~0"), collapse = " ")), "(*1)")
      BSEM.lambda <- c(BSEM.lambda, paste(names(loads)[[k]], "BY", paste(loads[[k]], collapse = " ")))
      BSEM.lambda <- c(BSEM.lambda, paste(names(loads)[[k]], "BY", paste0(rownames(lambda)[as.numeric(which(lambda[,k] != 1))], collapse = " "), paste0("(", paste0(names(loads)[[k]], "x", c(1, length(as.numeric(which(lambda[,k] != 1)))), collapse = "-"), ")")))
      BSEM.priors <- c(BSEM.priors, paste0(paste0(names(loads)[[k]], "x", c(1, length(as.numeric(which(lambda[,k] != 1)))), collapse = "-"), "~N(0,0.01)"))
    }
    EFA.lambda <- paste(paste(colnames(lambda), collapse = " "), "BY", paste(rownames(lambda), collapse = " "), "(*1)")

    if(any(ceiling((nchar(CFA.lambda) + 1) / 90) > 1)){
      CFA.tmp <- c()
      for(m in 1:length(CFA.lambda)){
        if(ceiling((nchar(CFA.lambda[m]) + 1) / 90) <= 1){
          CFA.tmp <- c(CFA.tmp, CFA.lambda[m])
        } else {
          cut <- 1
          CFA.lambda.trim <- c()
          for(s in 1:(ceiling((nchar(CFA.lambda[m]) + 1) / 90))){
            pos <- which.min(abs(gregexpr(" ", CFA.lambda[m])[[1]] - round(nchar(CFA.lambda[m]) / ceiling((nchar(CFA.lambda[m]) + 1) / 90)) * s))
            cut <- c(cut, gregexpr(" ", CFA.lambda[m])[[1]][pos])
            if(s == 1){
              CFA.lambda.trim <- c(CFA.lambda.trim, paste0(substr(CFA.lambda[m], cut[s], cut[s + 1] - 1), "\n"))
            } else if(s > 1 & s < ceiling((nchar(CFA.lambda[m]) + 1) / 90)) {
              CFA.lambda.trim <- c(CFA.lambda.trim, paste0(substr(CFA.lambda[m], cut[s] + 1, cut[s + 1] - 1), "\n"))
            } else {
              CFA.lambda.trim <- c(CFA.lambda.trim, paste0(substr(CFA.lambda[m], cut[s] + 1, nchar(CFA.lambda[m]))))
            }
          }
          CFA.lambda.trim <- paste(CFA.lambda.trim, collapse = " ")
          CFA.tmp <- c(CFA.tmp, CFA.lambda.trim)
        }
      }
      CFA.lambda <- CFA.tmp
    }

    if(ceiling((nchar(EFA.lambda) + 1) / 90) > 1){
      cut <- 1
      EFA.lambda.trim <- c()
      for(s in 1:(ceiling((nchar(EFA.lambda) + 1) / 90))){
        pos <- which.min(abs(gregexpr(" ", EFA.lambda)[[1]] - round(nchar(EFA.lambda) / ceiling((nchar(EFA.lambda) + 1) / 90)) * s))
        cut <- c(cut, gregexpr(" ", EFA.lambda)[[1]][pos])
        if(s == 1){
          EFA.lambda.trim <- c(EFA.lambda.trim, paste0(substr(EFA.lambda, cut[s], cut[s + 1] - 1), "\n"))
        } else if(s > 1 & s < ceiling((nchar(EFA.lambda) + 1) / 90)) {
          EFA.lambda.trim <- c(EFA.lambda.trim, paste0(substr(EFA.lambda, cut[s] + 1, cut[s + 1] - 1), "\n"))
        } else {
          EFA.lambda.trim <- c(EFA.lambda.trim, paste0(substr(EFA.lambda, cut[s] + 1, nchar(EFA.lambda))))
        }
      }
      EFA.lambda <- paste(EFA.lambda.trim, collapse = " ")
    }

    if(any(ceiling((nchar(EFAt.lambda) + 1) / 90) > 1)){
      EFAt.tmp <- c()
      for(m in 1:length(EFAt.lambda)){
        if(ceiling((nchar(EFAt.lambda[m]) + 1) / 90) <= 1){
          EFAt.tmp <- c(EFAt.tmp, EFAt.lambda[m])
        } else {
          cut <- 1
          EFAt.lambda.trim <- c()
          for(s in 1:(ceiling((nchar(EFAt.lambda[m]) + 1) / 90))){
            pos <- which.min(abs(gregexpr(" ", EFAt.lambda[m])[[1]] - round(nchar(EFAt.lambda[m]) / ceiling((nchar(EFAt.lambda[m]) + 1) / 90)) * s))
            cut <- c(cut, gregexpr(" ", EFAt.lambda[m])[[1]][pos])
            if(s == 1){
              EFAt.lambda.trim <- c(EFAt.lambda.trim, paste0(substr(EFAt.lambda[m], cut[s], cut[s + 1] - 1), "\n"))
            } else if(s > 1 & s < ceiling((nchar(EFAt.lambda[m]) + 1) / 90)) {
              EFAt.lambda.trim <- c(EFAt.lambda.trim, paste0(substr(EFAt.lambda[m], cut[s] + 1, cut[s + 1] - 1), "\n"))
            } else {
              EFAt.lambda.trim <- c(EFAt.lambda.trim, paste0(substr(EFAt.lambda[m], cut[s] + 1, nchar(EFAt.lambda[m]))))
            }
          }
          EFAt.lambda.trim <- paste(EFAt.lambda.trim, collapse = " ")
          EFAt.tmp <- c(EFAt.tmp, EFAt.lambda.trim)
        }
      }
      EFAt.lambda <- EFAt.tmp
    }

    if(any(ceiling((nchar(BSEM.lambda) + 1) / 60) > 1)){
      BSEM.tmp <- c()
      for(m in 1:length(BSEM.lambda)){
        if(ceiling((nchar(BSEM.lambda[m]) + 1) / 60) <= 1){
          BSEM.tmp <- c(BSEM.tmp, BSEM.lambda[m])
        } else {
          cut <- 1
          BSEM.lambda.trim <- nX <- c()
          for(s in 1:(ceiling((nchar(BSEM.lambda[m]) + 1) / 60))){
            pos <- which.min(abs(gregexpr(" ", BSEM.lambda[m])[[1]] - round(nchar(BSEM.lambda[m]) / ceiling((nchar(BSEM.lambda[m]) + 1) / 60)) * s))
            cut <- c(cut, gregexpr(" ", BSEM.lambda[m])[[1]][pos])
            if(grepl("-", BSEM.lambda[m])){
              if(s == 1){
                nameF <- strsplit(substr(BSEM.lambda[m], cut[s], cut[s + 1] - 1), " BY ")[[1]][1]
                nX <- c(nX, length(strsplit(strsplit(substr(BSEM.lambda[m], cut[s], cut[s + 1] - 1), " BY ")[[1]][2], " ")[[1]]))
                BSEM.lambda.trim <- c(BSEM.lambda.trim, paste0(substr(BSEM.lambda[m], cut[s], cut[s + 1] - 1), paste0(" (", paste0(nameF, "x", c(1, nX[1]), collapse = "-"), ")"), ";\n"))
              } else if(s > 1 & s < ceiling((nchar(BSEM.lambda[m]) + 1) / 60)) {
                nX <- c(nX, length(strsplit(substr(BSEM.lambda[m], cut[s] + 1, cut[s + 1] - 1), " ")[[1]]))
                BSEM.lambda.trim <- c(BSEM.lambda.trim, paste0(nameF, " BY ", substr(BSEM.lambda[m], cut[s] + 1, cut[s + 1] - 1), paste0(" (", paste0(nameF, "x", c(sum(nX[(1:s - 1)]) + 1, sum(nX)), collapse = "-"), ")"), ";\n"))
              } else {
                nX <- c(nX, length(strsplit(substr(BSEM.lambda[m], cut[s] + 1, nchar(BSEM.lambda[m])), " ")[[1]]) - 1)
                namesX <- paste(strsplit(substr(BSEM.lambda[m], cut[s] + 1, nchar(BSEM.lambda[m])), " ")[[1]][-(nX + 1)], collapse = " ")
                BSEM.lambda.trim <- c(BSEM.lambda.trim, paste0(nameF, " BY ", namesX, paste0(" (", paste0(nameF, "x", c(sum(nX[(1:s - 1)]) + 1, sum(nX)), collapse = "-"), ")")))
              }
            } else {
              if(s < ceiling((nchar(BSEM.lambda[m]) + 1) / 60)){
                BSEM.lambda.trim <- c(BSEM.lambda.trim, paste0(substr(BSEM.lambda[m], cut[s], cut[s + 1] - 1), "\n"))
              } else {
                BSEM.lambda.trim <- c(BSEM.lambda.trim, paste0(substr(BSEM.lambda[m], cut[s], nchar(BSEM.lambda[m]))))
              }
            }
          }
          BSEM.lambda.trim <- paste(BSEM.lambda.trim, collapse = " ")
          BSEM.tmp <- c(BSEM.tmp, BSEM.lambda.trim)
        }
      }
      BSEM.lambda <- BSEM.tmp
    }

    CFA.lambda <- paste0(paste(CFA.lambda, collapse = ";\n"), ";")
    EFA.lambda <- paste0(paste(EFA.lambda, collapse = ";\n"), ";")
    EFAt.lambda <- paste0(paste(EFAt.lambda, collapse = ";\n"), ";")
    BSEM.lambda <- paste0(paste(BSEM.lambda, collapse = ";\n"), ";")
    BSEM.priors <- paste0(paste(BSEM.priors, collapse = ";\n"), ";")

    m.cors <- ifelse(nrow(cors) > 1, paste0(paste(as.vector(apply(cors, 1, function(p) paste0(colnames(phi)[p[1]], " WITH ", colnames(phi)[p[2]], "@0"))), collapse = ";\n"), ";"), "")
    m.resi <- ifelse(nrow(resi) > 1, paste0(paste(as.vector(apply(resi, 1, function(p) paste0(colnames(psi)[p[1]], " WITH ", colnames(psi)[p[2]]))), collapse = ";\n"), ";"), "")

    if(nrow(cors) > 0){
      CFA.mplus <- paste(CFA.lambda, m.cors, sep = "\n")
      EFA.mplus <- paste(EFA.lambda, m.cors, sep = "\n")
      EFAt.mplus <- paste(EFAt.lambda, m.cors, sep = "\n")
      BSEM.mplus <- paste(BSEM.lambda, m.cors, sep = "\n")
      BSEM.mplus <- list(model = BSEM.mplus, priors = BSEM.priors)
    }
    if(nrow(resi) > 0){
      CFA.mplus <- paste(CFA.lambda, m.resi, sep = "\n")
      EFA.mplus <- paste(EFA.lambda, m.resi, sep = "\n")
      EFAt.mplus <- paste(EFAt.lambda, m.resi, sep = "\n")
      BSEM.mplus <- paste(BSEM.lambda, m.resi, sep = "\n")
      BSEM.mplus <- list(model = BSEM.mplus, priors = BSEM.priors)
    }
    if(nrow(cors) == 0 & nrow(resi) == 0){
      CFA.mplus <- CFA.lambda
      EFA.mplus <- EFA.lambda
      EFAt.mplus <- EFAt.lambda
      BSEM.mplus <- BSEM.lambda
      BSEM.mplus <- list(model = BSEM.mplus, priors = BSEM.priors)
    }

    CFA.mplus <- toupper(CFA.mplus)
    EFA.mplus <- toupper(EFA.mplus)
    EFAt.mplus <- toupper(EFAt.mplus)
    BSEM.mplus$model <- toupper(BSEM.mplus$model)
    BSEM.mplus$priors <- toupper(BSEM.mplus$priors)

    res <- list(CFA = CFA.mplus, EFA = EFA.mplus, EFAt = EFAt.mplus, BSEM = BSEM.mplus, specifications = specifications)
    class(res) <- "modFA"
    return(res)
  }
}

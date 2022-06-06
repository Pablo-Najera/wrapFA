#' Determinacy of factor score estimates
#'
#' @description Computes the determinacy of factor score estimates. Apart from the estimated determinacy (using the implied matrices), the
#' true determinacy can be computed if the generating matrices are provided (e.g., in simulation studies). Both the observed correlation matrix
#' and the implied correlation matrix (Beauducel, 2011) are used to calculate the determinacy.
#'
#' @param data A \emph{N} individuals x \emph{J} observed variables \code{matrix} or \code{data.frame}. Missing values need to be coded as \code{NA}. It can be \code{NULL} if the correlation matrix (\code{r}) is provided. Default is \code{NULL}.
#' @param r A \emph{J} observed variables x \emph{J} observed variables correlation \code{matrix} or \code{data.frame}. It can be \code{NULL} if the \code{data} are provided. Default is \code{NULL}.
#' @param est.lambda A \emph{J} observed variables x \emph{K} latent variables estimated factor loading \code{matrix} or \code{data.frame}.
#' @param est.phi A \emph{K} latent variables x \emph{K} latent variables estimated factor correlation \code{matrix} or \code{data.frame}.
#' @param est.psi A \emph{J} observed variables x \emph{J} observed variables estimated residual correlation \code{matrix} or \code{data.frame}.
#' @param true.lambda A \emph{J} observed variables x \emph{K} latent variables generating factor loading \code{matrix} or \code{data.frame}. Required to compute the true determinacy of factor score estimates. Useful for simulation studies. Default is \code{NULL}.
#' @param true.phi A \emph{K} latent variables x \emph{K} latent variables generating factor correlation \code{matrix} or \code{data.frame}. Required to compute the true determinacy of factor score estimates. Useful for simulation studies. Default is \code{NULL}.
#' @param align Align factor scores? It can be either \code{'auto'} or a \code{numeric vector}. If \code{'auto'}, factor scores will be aligned based on \code{true.lambda} (not aligned if \code{true.lambda = NULL}). Otherwise, the \code{numeric vector} must identify the order for the factor scores. Default is \code{'auto'}.
#'
#' @return \code{detFA} returns an object of class \code{detFA}.
#' \describe{
#' \item{\code{est.r}}{The estimated determinacy based on the observed correlation matrix (\code{matrix}).}
#' \item{\code{est.sigma}}{The estimated determinacy based on the implied correlation matrix (\code{matrix}).}
#' \item{\code{true.r}}{The true determinacy based on the observed correlation matrix. Only if \code{true.lambda} and \code{true.phi}  have been provided (\code{matrix}).}
#' \item{\code{true.r}}{The true determinacy based on the implied correlation matrix. Only if \code{true.lambda} and \code{true.phi} have been provided (\code{matrix}).}
#' }
#'
#' @references
#' Beauducel, A. (2011). Indeterminacy of factor score estimates in slightly misspecified confirmatory factor models. \emph{Journal of Modern Applied Statistical Methods}, \emph{10}, 583-598. https://doi.org/10.22237/jmasm/1320120900
#'
#' @export
#'
#' @examples
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
#' fd_efa <- detFA(HS, est.lambda = efa_geomin$lambda, est.phi = efa_geomin$phi, est.psi = efa_geomin$psi)
detFA <- function(data = NULL, r = NULL, est.lambda, est.phi, est.psi, true.lambda = NULL, true.phi = NULL, align = "auto"){

  #--------------------
  # 1. Check arguments
  #--------------------

  if(is.null(data) & is.null(r)){stop("Either {data} or {r} must be provided.")}
  if(!is.null(data)){if(!is.matrix(data) & !is.data.frame(data)){stop("{data} must be a matrix or data.frame.")}}
  if(!is.null(r)){if(!is.matrix(r) & !is.data.frame(r)){stop("{r} must be a matrix or data.frame.")}}
  if(!is.matrix(est.lambda) & !is.data.frame(est.lambda)){stop("{est.lambda} must be a matrix or data.frame.")}
  if(!is.matrix(est.phi) & !is.data.frame(est.phi)){stop("{est.phi} must be a matrix or data.frame.")}
  if(!is.matrix(est.psi) & !is.data.frame(est.psi)){stop("{est.psi} must be a matrix or data.frame.")}
  if(!is.null(true.lambda)){if(!is.matrix(true.lambda) & !is.data.frame(true.lambda)){stop("{true.lambda} must be a matrix or data.frame.")}}
  if(!is.null(true.phi)){if(!is.matrix(true.phi) & !is.data.frame(true.phi)){stop("{true.phi} must be a matrix or data.frame.")}}
  if(align != "auto" & !is.numeric(align)){stop("{align} must be either 'auto' or a numeric vector.")}

  #------------------------
  # 2. Align factor scores
  #------------------------

  if(is.null(r)){r <- cor(data, use = "pair")}
  if(align == "auto"){
    if(!is.null(true.lambda)){
      new.order <- fungible::faAlign(true.lambda, est.lambda, est.phi)
      est.lambda <- new.order$F2
      est.phi <- new.order$Phi2
    }
  } else {
    est.lambda <- est.lambda[,order]
    est.phi <- est.phi[,order]
    est.phi <- est.phi[order,]
  }
  sigma <- est.lambda %*% est.phi %*% t(est.lambda) + est.psi

  #---------------------
  # 3. Warning messages
  #---------------------

  if(any(est.lambda > 1)){warning("Heywood case detected.")}
  if(min(eigen(r)$values) < .Machine$double.eps){
    warning(paste0("Observed correlation matrix is singular. Moore-Penrose generalized inverse of the observed correlation matrix has been done."))
    inv.r <- MASS::ginv(r)
  } else {
    inv.r <- solve(r)
  }
  if(min(eigen(sigma)$values) < .Machine$double.eps){
    warning(paste0("Reproduced correlation matrix is singular. Moore-Penrose generalized inverse of the reproduced correlation matrix has been done."))
    inv.sigma <- MASS::ginv(sigma)
  } else {
    inv.sigma <- solve(sigma)
  }

  #--------------------------
  # 4. Estimated reliability
  #--------------------------

  est.r.p <- sqrt(diag(est.phi %*% t(est.lambda) %*% inv.r %*% est.lambda %*% est.phi))
  est.r.p2 <- est.r.p^2
  est.r.minp <- 2 * est.r.p2 - 1
  est.r <- list(p = est.r.p, p2 = est.r.p2, minp = est.r.minp)

  est.sigma.p <- sqrt(diag(est.phi %*% t(est.lambda) %*% inv.sigma %*% est.lambda %*% est.phi))
  est.sigma.p2 <- est.sigma.p^2
  est.sigma.minp <- 2 * est.sigma.p2 - 1
  est.sigma <- list(p = est.sigma.p, p2 = est.sigma.p2, minp = est.sigma.minp)

  #---------------------
  # 5. True reliability
  #---------------------

  true.r <- true.sigma <- NULL
  if(!is.null(true.lambda) & !is.null(true.phi)){
    true.S <- true.lambda %*% true.phi
    est.S <- est.lambda %*% est.phi

    W.r <- inv.r %*% est.S
    C.r <- t(W.r) %*% r %*% W.r
    L.r <- diag(sqrt(diag(C.r)))
    ifelse(min(eigen(L.r)$values) < .Machine$double.eps, inv.L.r <- MASS::ginv(L.r), inv.L.r <- solve(L.r))
    true.r.p <- diag(t(true.S) %*% W.r %*% inv.L.r)
    true.r.p2 <- true.r.p^2
    true.r.minp <- 2 * true.r.p2 - 1
    true.r <- list(p = true.r.p, p2 = true.r.p2, minp = true.r.minp)

    W.s <- inv.sigma %*% est.S
    C.s <- t(W.s) %*% sigma %*% W.s
    L.s <- diag(sqrt(diag(C.s)))
    ifelse(min(eigen(L.s)$values) < .Machine$double.eps, inv.L.s <- MASS::ginv(L.s), inv.L.s <- solve(L.s))
    true.sigma.p <- diag(t(true.S) %*% W.s %*% inv.L.s)
    true.sigma.p2 <- true.sigma.p^2
    true.sigma.minp <- 2 * true.sigma.p2 - 1
    true.sigma <- list(p = true.sigma.p, p2 = true.sigma.p2, minp = true.sigma.minp)
  }

  #-------------------
  # 6. Export results
  #-------------------

  res <- list(est.r = est.r, est.sigma = est.sigma, true.r = true.r, true.sigma = true.sigma)
  class(res) <- "detFA"
  return(res)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detFA.R
\name{detFA}
\alias{detFA}
\title{Determinacy of factor score estimates}
\usage{
detFA(
  data = NULL,
  r = NULL,
  est.lambda,
  est.phi,
  est.psi,
  true.lambda = NULL,
  true.phi = NULL,
  align = "auto"
)
}
\arguments{
\item{data}{A \emph{N} individuals x \emph{J} observed variables \code{matrix} or \code{data.frame}. Missing values need to be coded as \code{NA}. It can be \code{NULL} if the correlation matrix (\code{r}) is provided. Default is \code{NULL}.}

\item{r}{A \emph{J} observed variables x \emph{J} observed variables correlation \code{matrix} or \code{data.frame}. It can be \code{NULL} if the \code{data} are provided. Default is \code{NULL}.}

\item{est.lambda}{A \emph{J} observed variables x \emph{K} latent variables estimated factor loading \code{matrix} or \code{data.frame}.}

\item{est.phi}{A \emph{K} latent variables x \emph{K} latent variables estimated factor correlation \code{matrix} or \code{data.frame}.}

\item{est.psi}{A \emph{J} observed variables x \emph{J} observed variables estimated residual correlation \code{matrix} or \code{data.frame}.}

\item{true.lambda}{A \emph{J} observed variables x \emph{K} latent variables generating factor loading \code{matrix} or \code{data.frame}. Required to compute the true determinacy of factor score estimates. Useful for simulation studies. Default is \code{NULL}.}

\item{true.phi}{A \emph{K} latent variables x \emph{K} latent variables generating factor correlation \code{matrix} or \code{data.frame}. Required to compute the true determinacy of factor score estimates. Useful for simulation studies. Default is \code{NULL}.}

\item{align}{Align factor scores? It can be either \code{'auto'} or a \code{numeric vector}. If \code{'auto'}, factor scores will be aligned based on \code{true.lambda} (not aligned if \code{true.lambda = NULL}). Otherwise, the \code{numeric vector} must identify the order for the factor scores. Default is \code{'auto'}.}
}
\value{
\code{detFA} returns an object of class \code{detFA}.
\describe{
\item{\code{est.r}}{The estimated determinacy based on the observed correlation matrix (\code{matrix}).}
\item{\code{est.sigma}}{The estimated determinacy based on the implied correlation matrix (\code{matrix}).}
\item{\code{true.r}}{The true determinacy based on the observed correlation matrix. Only if \code{true.lambda} and \code{true.phi}  have been provided (\code{matrix}).}
\item{\code{true.r}}{The true determinacy based on the implied correlation matrix. Only if \code{true.lambda} and \code{true.phi} have been provided (\code{matrix}).}
}
}
\description{
Computes the determinacy of factor score estimates. Apart from the estimated determinacy (using the implied matrices), the
true determinacy can be computed if the generating matrices are provided (e.g., in simulation studies). Both the observed correlation matrix
and the implied correlation matrix (Beauducel, 2011) are used to calculate the determinacy.
}
\examples{
\dontrun{
library(MBESS)
data(HS)
HS <- HS[, -c(1:8)]
colnames(HS) <- paste0("t", 1:26)
modHS <- data.frame(spatial = c(rep(1, 4), rep(0, 20), rep(1, 2)),
                    verbal = c(rep(0, 4), rep(1, 5), rep(0, 17)),
                    speed = c(rep(0, 9), rep(1, 4), rep(0, 13)),
                    memory = c(rep(0, 13), rep(1, 6), rep(0, 7)),
                    maths = c(rep(0, 19), rep(1, 5), rep(0, 2)), row.names = colnames(HS))
modHS.lav <- modFA(lambda = modHS, software = "lavaan", keep.names = TRUE)
efa_geomin <- EFA(model = modHS.lav$EFA, data = HS)
fd_efa <- detFA(HS, est.lambda = efa_geomin$lambda, est.phi = efa_geomin$phi, est.psi = efa_geomin$psi)
}
}
\references{
Beauducel, A. (2011). Indeterminacy of factor score estimates in slightly misspecified confirmatory factor models. \emph{Journal of Modern Applied Statistical Methods}, \emph{10}, 583-598. https://doi.org/10.22237/jmasm/1320120900
}

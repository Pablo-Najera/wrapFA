% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSEM.R
\name{BSEM}
\alias{BSEM}
\title{Bayesian structural equation modeling (BSEM)}
\usage{
BSEM(
  model,
  data,
  categorical = NULL,
  mplus.path = NULL,
  rm.files = FALSE,
  max.iter = 60000,
  n.cores = 1,
  n.chains = 2,
  suppressMessages = TRUE
)
}
\arguments{
\item{model}{A \code{list} containing two \code{character} elements: the user-specified BSEM model and the priors specification using the \emph{Mplus} syntax .}

\item{data}{A \emph{N} individuals x \emph{J} observed variables \code{matrix} or \code{data.frame}. Missing values need to be coded as \code{NA}.}

\item{categorical}{A \code{vector} containing the number or name of the observed variables that need to be treated as categorical. Default is \code{NULL}, which treates all observed variables as continuous.}

\item{mplus.path}{A \code{character} indicating the path where the the files generated by \emph{Mplus} will be stored, including the name of the files. Default is \code{NULL}, which stores the files in the current working directory with the name \emph{wrapBSEM}.}

\item{rm.files}{A \code{logical} indicating whether the stored files should be deleted. The default is \code{FALSE}.}

\item{max.iter}{A \code{double} indicating the maximum number of iterations. Default is 60000.}

\item{n.cores}{A \code{double} indicating the number of CPU processors to be used for the estimation of the model. Default is 1.}

\item{n.chains}{A \code{double} indicating the number of MCMC chains. Default is 2.}

\item{suppressMessages}{A \code{logical} indicating whether the messages from MplusAutomation should be silent. Default is \code{TRUE}.}
}
\value{
\code{BSEM} returns an object of class \code{wrapFA}.
\describe{
\item{\code{model.fit}}{Model fit indices (\code{vector}).}
\item{\code{lambda}}{The standardized estimated factor loading matrix (\code{matrix}).}
\item{\code{phi}}{The standardized estimated factor correlation matrix (\code{matrix}).}
\item{\code{psi}}{The standardized estimated error correlation matrix (\code{matrix}).}
\item{\code{lambda.p}}{The p-value associated to each free parameter in the estimated factor loading matrix (\code{matrix}).}
\item{\code{phi.p}}{The p-value associated to each free parameter in the estimated factor correlation matrix (\code{matrix}).}
\item{\code{psi.p}}{The p-value associated to each free parameter in the estimated error correlation matrix (\code{matrix}).}
\item{\code{out}}{The complete output provided by \emph{Mplus} (\code{list}).}
\item{\code{specifications}}{Specifications used to run the function (\code{list}).}
}
}
\description{
Fit a BSEM (BSEM; Muthén & Asparouhov, 2012) via \emph{MplusAutomation}.
}
\examples{
\dontrun{
BFI <- psych::bfi
BFI <- BFI[,1:25]
modBFI <- data.frame(agree = c(rep(1, 5), rep(0, 20)),
                       consc = c(rep(0, 5), rep(1, 5), rep(0, 15)),
                       extra = c(rep(0, 10), rep(1, 5), rep(0, 10)),
                       neuro = c(rep(0, 15), rep(1, 5), rep(0, 5)),
                       openn = c(rep(0, 20), rep(1, 5)), row.names = colnames(BFI))
modBFI.mpl <- modFA(lambda = modBFI, software = "mplus", keep.names = TRUE)
bsem <- BSEM(model = modBFI.mpl$BSEM, data = BFI, categorical = 1:25, n.cores = 2)
}
}
\references{
Muthén, B. & Asparouhov, T. (2012). Bayesian structural equation modeling: A more flexible representation of substantive theory. \emph{Psychological Methods}, \emph{17}, 313-335. https://doi.org/10.1037/a0026802
}

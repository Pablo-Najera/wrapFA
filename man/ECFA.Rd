% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECFA.R
\name{ECFA}
\alias{ECFA}
\title{EFA-based confirmatory factor analysis (ECFA)}
\usage{
ECFA(
  fit = NULL,
  model = NULL,
  data,
  method = "r2",
  r2.args = list(phi.grid = seq(0.7, 0.99, 0.01), selector = "BIC"),
  categorical = NULL,
  estimator = NULL,
  software = NULL,
  mimic.mplus = TRUE,
  mplus.path = NULL,
  rm.files = FALSE,
  max.iter = 10000,
  suppressMessages = TRUE
)
}
\arguments{
\item{fit}{A \code{wrapFA} object from \code{EFA()} or \code{BSEM()} functions. Default is \code{NULL}, which means that a user-specified model will be used.}

\item{model}{A \code{character} string containing the EFA model using either the \emph{Mplus} or the \emph{lavaan} syntax. Default is \code{NULL}, which means that the argument \code{fit} will be used.}

\item{data}{A \emph{N} individuals x \emph{J} observed variables \code{matrix} or \code{data.frame}. Missing values need to be coded as \code{NA}.}

\item{method}{A \code{character} indicating the method to detect the relevant loadings. It can be \code{"r2"} for the R-squared procedure (Citation, Year), \code{"nominal"} for statistically significant loadings with \emph{alpha} = 0.05, or \code{"bonferroni"} for statistically significant loadings using the Bonferroni correction. Default is \code{"r2"}.}

\item{r2.args}{A \code{list} of arguments about the R-squared procedure:
\describe{
\item{\code{phi.grid}}{A \code{numeric} vector indicating the phi parameter space to search. Default is from 0.70 to 0.99, in steps of 0.01.}
\item{\code{selector}}{A \code{character} indicating fit index to use to select the final model. It can be \code{"ChiSq/df"}, \code{"CFI"}, \code{"TLI"}, \code{"AIC"}, \code{"BIC"}, \code{"aBIC"}, \code{"RMSEA"}, \code{"SRMR"}, and \code{"AICC"}. Default is \code{"BIC"}.}
}}

\item{categorical}{A \code{vector} containing the number or name of the observed variables that need to be treated as categorical. Default is \code{NULL}, which treates all observed variables as continuous.}

\item{estimator}{A \code{character} indicating the estimator to be used to fit the model. It can be \code{"ML"}, \code{"MLM"}, \code{"MLMV"}, \code{"MLMVS"}, \code{"MLF"}, \code{"MLR"}, \code{"GLS"}, \code{"DWLS"}, \code{"WLS"}, \code{"WLSM"}, \code{"WLSMV"}, \code{"ULS"},\code{"ULSM"}, and \code{"ULSMV"}. Check \code{?lavaan::lavOptions} or the Mplus User's Guide for more information. Default is \code{NULL}, which uses \code{"ML"} if all observed variables are continuous and \code{"WLSMV"} if any observed variable is categorical.}

\item{software}{A \code{character} indicating the software to use for fitting the CFA model. It can be \code{"lavaan"} or \code{"mplus"}. Only required if design matrices are used in the \code{model} argument; the software is automatically detected if a syntax is used. Default is \code{NULL}.}

\item{mimic.mplus}{A \code{logical} indicating whether \emph{lavaan} should mimic \emph{Mplus} output. Only applicable to \emph{lavaan}. Default is \code{TRUE}.}

\item{mplus.path}{A \code{character} indicating the path where the the files generated by Mplus will be stored, including the name of the files. Default is \code{NULL}, which stores the files in the current working directory with the name \emph{wrapECFA}.}

\item{rm.files}{A \code{logical} indicating whether the stored files should be deleted. The default is \code{FALSE}.}

\item{max.iter}{A \code{double} indicating the maximum number of iterations. Default is 10000.}

\item{suppressMessages}{A \code{logical} indicating whether the messages from MplusAutomation should be silent. Default is \code{TRUE}.}
}
\value{
\code{ECFA} returns an object of class \code{wrapFA}.
\describe{
\item{\code{model.fit}}{Model fit indices (\code{vector}).}
\item{\code{lambda}}{The standardized estimated factor loading matrix (\code{matrix}).}
\item{\code{phi}}{The standardized estimated factor correlation matrix (\code{matrix}).}
\item{\code{psi}}{The standardized estimated error correlation matrix (\code{matrix}).}
\item{\code{lambda.p}}{The p-value associated to each free parameter in the estimated factor loading matrix (\code{matrix}).}
\item{\code{phi.p}}{The p-value associated to each free parameter in the estimated factor correlation matrix (\code{matrix}).}
\item{\code{psi.p}}{The p-value associated to each free parameter in the estimated error correlation matrix (\code{matrix}).}
\item{\code{MI}}{Information about the MI and SEPC of the fixed parameters (\code{matrix}).}
\item{\code{out}}{The complete output provided by either \emph{lavaan} or \emph{Mplus} (\code{list}).}
\item{\code{specifications}}{Specifications used to run the function (\code{list}).}
}
}
\description{
Fit an ECFA model (Nájera et al., 2022) via \emph{lavaan} or \emph{MplusAutomation}.
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
ecfa_r2 <- ECFA(fit = efa_geomin, data = HS, method = "r2", r2.args = list(phi.grid = seq(0.85, 0.95, 0.01), selector = "BIC"))
}
}
\references{
Nájera, P., Abad, F. J., & Sorrel, M. A. (2022). Is EFA always to be preferred? A systematic comparison of factor analytic techniques throughout the confirmatory-exploratory continuum. \emph{Manuscript submitted for publication}.
}

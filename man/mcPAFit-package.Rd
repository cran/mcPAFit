\name{mcPAFit-package}
\alias{mcPAFit-package}
\docType{package} 
\title{Estimating Preferential Attachment from a Single Network Snapshot by Markov Chain Monte Carlo}
\description{
A Markov chain Monte Carlo method is provided to estimate the preferential attachment function from a single network snapshot. Conventional methods require the complete information about the appearance order of all nodes and edges in the network. This package incorporates the appearance order into the state space and estimates it together with the preferential attachment function. Auxiliary variables are introduced to facilitate fast Gibbs sampling.

}
\details{
\tabular{ll}{
Package: \tab mcPAFit\cr
Type: \tab Package\cr
Version: \tab 0.1.3\cr
Date: \tab 2016-05-25\cr
License: \tab GPL-3 \cr
}
\itemize{
\item mcPAFit: estimates the Preferential Attachment function in a temporal complex network without temporal information.
\item create_sim_data: generate simulation data.
}
}
\author{
Thong Pham, Paul Sheridan, Hidetoshi Shimodaira.
Maintainer: Thong Pham \email{thongpham@thongpham.net}
}


\references{
1. Pham, T. and Sheridan, P. and Shimodaira, H. (2015). mcPAFit: Nonparametric Measurement of Preferential Attachment and Fitness from a Single Network Snapshot. Conference on Complex Systems 2015, September 2015, Arizona.
}

\keyword{package}
\keyword{Preferential Attachment function}
\keyword{Markov chain Monte Carlo}
\keyword{Gibbs sampling}
\keyword{missing timeline}
\examples{
\dontrun{
library("mcPAFit")
}
}

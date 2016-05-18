\name{mcPAFit-package}
\alias{mcPAFit-package}
\docType{package} 
\title{Estimating Node Fitness from a Single Network Snapshot by Markov Chain Monte Carlo}
\description{
A Markov chain Monte Carlo method is provided to estimate fitness from a single network snapshot. Conventional methods require the complete information about the appearance order of all nodes and edges in the network. This package incorporates the appearance order into the state space and estimate it together with node fitness. For estimation of fitness, auxiliary variables are introduced to facilitate fast Gibbs sampling for fitness.

If you use this package for your projects, please run the command: citation("mcPAFit") for a list of references you should cite.
}
\details{
\tabular{ll}{
Package: \tab mcPAFit\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2016-05-14\cr
License: \tab GPL-3 \cr
}
\itemize{
\item mcPAFit: estimates the Preferential Attachment function and fitness function in a temporal complex network.
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
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{MM algorithims}
\keyword{Bianconi-Barabasi model}
\keyword{fitness}
\keyword{maximum likelihood function}
\examples{
\dontrun{
library("mcPAFit")
}
}

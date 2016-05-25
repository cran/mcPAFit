\name{create_sim_data}
\alias{create_sim_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating simulated networks
}
\description{
A function generates a complex network from the Bianconi-Barabasi model and a random order to create a random network.
}
\usage{
create_sim_data(mode = 1, sat_at = 10, alpha = 0.5, beta = 1, N = 200)
}
\arguments{
\item{mode}{
Integer. Indicates the attachment function to be used in generating the network. If mode == 1, the attachment function is \eqn{A_k = k^\alpha}. If mode == 2, the attachment function is \eqn{A_k = min(k,sat_at)^\alpha}. If mode == 3, the attachment function is \eqn{A_k = \alpha log (k)^\beta}. Default value is \eqn{1}.
} 
\item{sat_at}{
Integer. This is the saturation position sat_at in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}.
}
\item{alpha}{
Numeric. If mode == 1, this is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. If mode == 2, this is the attachment exponenet in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}. If mode == 3, this is the alpha in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{beta}{
Numeric. This is the beta in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{N}{
Integer. Number of nodes in the final network.
}
}
\value{
An object of class mcPAFit.Sim containing the true network as well as the random-shuffled final snapshot together with the true timeline.
}
\author{
Thong Pham
\email{thongpham@thongpham.net}
}
\references{
1. Pham, T. and Sheridan, P. and Shimodaira, H. (2015). mcPAFit: Nonparametric Measurement of Preferential Attachment and Fitness from a Single Network Snapshot. Conference on Complex Systems 2015, September 2015, Arizona.
}
\examples{
library("mcPAFit")
data <- create_sim_data(N = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{Preferential Attachment function}
\keyword{Markov chain Monte Carlo}
\keyword{Gibbs sampling}
\keyword{missing timeline}


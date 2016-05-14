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
create_sim_data(shape = 5, rate = 5, N = 200)
}
\arguments{
  \item{shape}{
Numeric. Shape parameter for the true distribution of fitness.
}
\item{rate}{
Numeric. Rate parameter for the true distribution of fitness.
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
\keyword{Fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{Markov chain Monte Carlo}% __ONLY ONE__ keyword per line
\keyword{Gibbs sampling}
\keyword{Bianconi-Barabasi model}
\keyword{Fitness}


\name{simple_net}
\alias{simple_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating simulated networks
}
\description{
This function generates a complex network from the Barabasi-Albert model and a random order to create a random network.
}
\usage{
simple_net(time_step, 
           num_seed           = 2     , 
           p                  = 0.5   ,
           alpha              = 1     ,
           alpha_out          = 0   )
}
\arguments{
\item{time_step}{
Integer. Indicates the number of time-steps.
} 
\item{num_seed}{
Integer. Number of nodes in the initial network.
}
\item{p}{
Numeric. Proportion of number of nodes.
}
\item{alpha}{
Numeric. The attachment exponent for the in-degree.
}
\item{alpha_out}{
Numeric. The attachment exponent for the out-degree.
}
}
\value{
An three-column matrix containing the generated network.
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
data <- simple_net(100) # Time-steps = 100
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{Preferential Attachment function}
\keyword{Markov chain Monte Carlo}
\keyword{missing timeline}


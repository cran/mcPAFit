\name{mcPAFit}
\alias{mcPAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimates the Preferential Attachment function in a temporal complex network without temporal information.
}
\description{
This function takes an input file name and then output the MCMC results into files.  
}
\usage{
mcPAFit(file_name, burn_in = 1000, 
        needed_sample = 1000, 
        skip = 1, B = 1,
        step_s = 0.1, lambda = 1, G = 10000, s = 1,
        h_s_shape = 5, h_s_rate = 5,
        not_start_A = 1, file_A = "", not_start_f = 1,
        file_f = "", filename = "network", only_active_node = 1, random = 0)
}
\arguments{
  \item{file_name}{
String. Name of the file containing the network snapshot.
}
\item{burn_in}{
Integer. Number of burn-in steps. Default is 1000.
}
\item{needed_sample}{
Integer. Number of needed sample. Number of iterations = Skip * needed_sample + burn_in. Default is 1000.
}
\item{skip}{
Integer. Number of skips. Large skips reduce correleration between successive samples. Default is 1.
}
\item{B}{
Integer. Number of combinations of odd-sweep and even-sweep inside a complete sweep. Default is 1.
}
\item{step_s}{
Numeric. Step size for the Metropolis-Hastings random walk for s.
}
\item{lambda}{
Integer. Mean parameter of the Poisson prior of number of edges at a time.
}
  \item{G}{
Integer. Number of bins. Default value is 10000.
}
\item{s}{
Numeric. Parameter of the gamma distribution of PA.
}
\item{h_s_shape}{
Numeric. The shape of the prior for s.
}

\item{h_s_rate}{
Numeric. The rate of the prior for s.
}
\item{not_start_A}{
Integer. Indicates whether to use a custom PA function supplied by user. 1 is not to use. Default is 1. 
}
\item{file_A}{String. Name of the file containing the custom PA function in the case not_start_A == 0. Default is "".}
\item{not_start_f}{
Integer. Indicates whether to use a custom starting position for f. 1 is not to use. Default is 1. 
}
\item{file_f}{
String. Name of the file containing the custom fitness in the case not_start_f == 0. Default is "".
}
\item{filename}{String. Main part of names of the output files. Default is "network".}
\item{only_active_node}{Integer. If 1, then only nodes with degrees greater than 0 are considered. Default is 1.}
\item{random}{Integer. If 1, then the program randomizes the order of edges for starting value. Default is 0.}
}
\value{
The function writes directly to output files.
}
\author{
Thong Pham
\email{thongpham@thongpham.net}
}
\references{
1. Pham, T. and Sheridan, P. and Shimodaira, H. (2015). mcPAFit: Nonparametric Measurement of Preferential Attachment and Fitness from a Single Network Snapshot. Conference on Complex Systems 2015, September 2015, Arizona.
}
\examples{
library(mcPAFit)
N      <- 1000 #number of nodes in the network 
M      <- 500  #number of needed sample
lag   <-  1
start <-  100

data   <- create_sim_data(N = N, alpha = 1)
write.table(data$random[,1:2], row.names = FALSE,
            col.names = FALSE, file = "input.txt")
system.time(mcPAFit(file_name = "input.txt",burn_in = 0, s = 5, skip = 1, 
            needed_sample = M, G = 20,only_active_node = 1, random = 0))
# plot log posterior
to.read = file("log_pos_network.binary","rb")
log_pos = readBin(to.read, numeric(),n = M);
plot(log_pos,type = "l",xlab = "Iteration", ylab = "Log Posterior",cex.axis = 1.5, cex.lab = 1.5);
########

########
#plot b
to.read = file("b_out_network.binary","rb")
b = readBin(to.read, numeric(),n = M);
used  <- seq(start,M,1)
plot(b[used],type = "l", xlab = "Iteration", ylab = "b", cex.lab =2 , cex.axis =2);
acf(b[used])
mean(b[used])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{Preferential Attachment function}
\keyword{Markov chain Monte Carlo}
\keyword{Gibbs sampling}
\keyword{missing timeline}



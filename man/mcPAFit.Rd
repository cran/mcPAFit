\name{mcPAFit}
\alias{mcPAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate node fitness from a single network snapshot assuming the Bianconi-Barabasi model.
}
\description{
This function takes an input file name and then output the MCMC results into files.  
}
\usage{
mcPAFit(file_name, burn_in = 1000, 
        needed_sample = 1000, 
        skip = 1, B = 1,
        step_s = 0.1, lambda = 1, G = 10000, s = 5,
        h_s_shape = 5, h_s_rate = 5,
        not_start_A = 1, file_A = "", not_start_f = 1,
        file_f = "", filename = "network", only_active_node = 1)
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
Numeric. Parameter of the distribution of fitness.
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
N      <- 100 #number of nodes in the network 
M      <- 50  #number of needed sample
data   <- create_sim_data(N = N)
write.table(data$random[,1:2],row.names = FALSE,
            col.names = FALSE,file = "input.txt")
mcPAFit(file_name = "input.txt",burn_in = 0,needed_sample = M)
# plot log posterior
to.read = file("log_pos_network.binary","rb")
log_pos = readBin(to.read, numeric(),n = M);
plot(log_pos,type = "l",xlab = "Iteration", ylab = "Log Posterior",cex.axis = 1.5, cex.lab = 1.5);
########
start  <- 20 # discard all previous samples
########
#plot s
to.read = file("s_out_network.binary","rb")
s = readBin(to.read, numeric(),n = M);
used  <- seq(start,M,1)
plot(s[used],type = "l", xlab = "Iteration", ylab = "s", cex.lab =2 , cex.axis =2);
acf(s[used])
mean(s[used])
# get samples of fitness
to.read = file("f_out_network.binary","rb")
size_f <- N
f = matrix(readBin(to.read, numeric(),n = size_f*M,size = 8),nrow = M,ncol = size_f,byrow = TRUE);
f_mc   <- apply(f[start:M,],2,mean)
#### Node id ###
a <- read.table(file = "input.txt")
count <- 1
node_id  <- unique(as.vector(t(a)))
#plot true fitness versus estimated fitness
true        <- data$true$fitness
names(true) <- 1:max(node_id)
final_deg   <- table(data$random[,2])
ok          <- names(which(final_deg >= 2))
true        <- true[ok]
names(f_mc) <- node_id
f_mc        <- f_mc[as.character(ok)]
c_mc        <- lm(true ~ f_mc + 0)$coefficients[1]
f_temp      <- c_mc * f
colnames(f_temp) <- node_id
sd          <- apply(f_temp[start:M,ok],2,sd)
sd          <- sd/sqrt(dim(f_temp)[1])
upper       <- c_mc * f_mc + 2*sd
lower       <- c_mc * f_mc - 2*sd
plot(true,c_mc*f_mc,cex = 2,xlim = c(min(true,c_mc*f_mc), max(true,c_mc*f_mc)),
     ylim = c(min(true,c_mc*f_mc), max(true,c_mc*f_mc)),
     cex.lab = 1, cex.axis = 1,col = "blue", pch = 20, 
     xlab = " True fitness",
     ylab = "Estimated fitness")
abline(a = 0, b= 1,lwd = 2)
arrows(y0 = lower, x0 = true, y1 = upper, x1 = true, 
       code = 3,angle = 90, length = 0,
       col = "red",lwd = 1)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{MM algorithims}
\keyword{Bianconi-Barabasi model}
\keyword{fitness}
\keyword{maximum likelihood function}% __ONLY ONE__ keyword per line

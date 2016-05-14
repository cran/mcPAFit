#create simulated data
create_sim_data <- function(shape = 5, rate = 5, N = 200){
    data  <- GenerateNet(N = N, m = 1, mode = 1, alpha = 1,
                         shape = shape, 
                         rate = rate, num_seed = 2)
    order      <- sample(1:dim(data$graph)[1], dim(data$graph)[1],replace = FALSE)
    data_new   <- data$graph[order,1:2]
    data_new   <- cbind(data_new,0:(dim(data_new)[1] - 1))
    result     <- list(true = data, order = order, random = data_new)
    class(result) <- "mcPAFit.Sim"  
    return(result)
}


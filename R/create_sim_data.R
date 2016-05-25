#create simulated data
create_sim_data <- function(mode = 1, sat_at = 10, alpha = 0.5, beta = 1, N = 200){
    data  <- GenerateNet(N = N, m = 1, mode = mode, alpha = alpha, sat_at = sat_at, beta = beta,
                         num_seed = 2)
    order      <- sample(1:dim(data$graph)[1], dim(data$graph)[1],replace = FALSE)
    data_new   <- data$graph[order,1:2]
    data_new   <- cbind(data_new,0:(dim(data_new)[1] - 1))
    result     <- list(true = data, order = order, random = data_new)
    class(result) <- "mcPAFit.Sim"  
    return(result)
}


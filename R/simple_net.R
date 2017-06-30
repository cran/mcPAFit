# function to generate simulated network  
simple_net <-
function(time_step, 
         num_seed           = 2     , 
         p                  = 0.5   , # with probability p add a new node without any edge, with probability 1 - p add a new edge between two existing nodes
         alpha              = 1     ,
         alpha_out          = 0   
      ){
   # time_step: number of time-steps to grow the network

   if (num_seed < 2)
       stop("num_seed too small")
   if ((p < 0) || (time_step < 0))
       stop("The parameters must be non-negative")  

  
   num_of_row  <- num_seed + time_step 
                  
   edge_list   <- matrix(nrow = num_of_row, ncol = 3,0)
  
  
    
    edge_list_index      <- 1
    current_time_step    <- 0                 # current time_step
    for (n in 1:num_seed) {
         edge_list[edge_list_index,] <- c(n, -1, current_time_step);
         current_time_step           <- current_time_step + 1
         edge_list_index             <- edge_list_index + 1
    }
    #print(edge_list[1:(edge_list_index - 1),])
    
    degree                 <- rep(0,num_seed)
    names(degree)          <- as.integer(1:num_seed) 
    
    
    degree_out             <- rep(0,num_seed)
    names(degree_out)      <- as.integer(1:num_seed)

    P               <- degree
    P[P == 0]       <- 1
    seed.graph.size <- length(P)
    count           <- 0

    P_out             <- degree_out
    P_out[P_out == 0] <- 1
    sum_m             <- 0
  
    current_length_edge  <- dim(edge_list)[1] # current size of the edge_list matrix
    

    current_node <- 1:num_seed
    max_id       <- num_seed 
    current_time_step <- current_time_step - 1
    for (t in 1:time_step) {	
        P.sum             <- sum(P^alpha)
        node.weights      <- P^alpha/P.sum
        
        P_out.sum         <- sum(P_out^alpha_out)
        node_out.weights  <- P_out^alpha_out / P_out.sum
        
        current_time_step <- current_time_step + 1
        u                 <- runif(1)
        if (u < p)
            dice <- 1
        else dice <- 0 
        if (1 == dice) {
        ### adding a new node without any edge
           degree        <- c(degree,0)
           current_node  <- c(current_node,max_id + 1)
           max_id        <- max_id + 1
           names(degree) <- as.character(as.integer(current_node))
           
           degree_out    <- c(degree_out,0)
           names(degree_out) <- as.character(as.integer(current_node))
           
           from_node     <- max_id
           to_node       <- -1
           time_stamp    <- current_time_step
        } else {
          ### adding one new edge between existing nodes    
          to_node      <- sample(current_node, 1, prob = node.weights, replace = TRUE)     
          from_node    <- sample(current_node, 1, prob = node_out.weights, replace = TRUE)
          time_stamp   <- current_time_step
          # update the degree vector
          if (0 != length(to_node)) {
            temp    <- table(to_node)
            for(i in 1:length(temp)) { 
              num_edge            <- as.numeric(temp[i]) 
              node_name           <- as.character(as.integer(labels(temp[i])))
              degree[node_name]   <- degree[node_name] + num_edge # Update degrees.
            }
          }
          
          # update the degree_out vector
          if (0 != length(from_node)) {
            temp_out <- table(from_node)
            for(i in 1:length(temp_out)) { 
              num_edge_out              <- as.numeric(temp_out[i]) 
              node_name_out             <- as.character(as.integer(labels(temp_out[i])))
              degree_out[node_name_out] <- degree_out[node_name_out] + num_edge_out # Update degrees.
            }
          }
        }
        if (0 != length(to_node)) {
            edge_list[edge_list_index:(edge_list_index + length(to_node) - 1),] <- cbind(from_node,to_node,time_stamp)  
            edge_list_index <- edge_list_index + length(to_node)
        }
       P              <- degree
       P[degree == 0] <- 1  
       
       P_out          <- degree_out
       P_out[degree_out == 0] <- 1
   }
    if (edge_list_index < dim(edge_list)[1])
        edge_list <- edge_list[-(edge_list_index:dim(edge_list)[1]),]
  return(edge_list)
}

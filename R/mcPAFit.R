# function to estimate fitness from a single snapshot  2016-5-31 Thong Pham


mcPAFit <- function(file_name, burn_in = 1000, needed_sample = 1000, skip = 1, B = 1,
                    step_s = 0.1, lambda = 1, G = 10000, s = 5, h_s_shape = 5, h_s_rate = 5,
                    not_start_A = 1, file_A = "", not_start_f = 1,file_f = "", filename = "network") {
  m <- c(file_name, burn_in , needed_sample , skip , B ,
         step_s , lambda , G, s , h_s_shape , h_s_rate ,
         not_start_A, file_A, not_start_f, file_f, filename)
  .mcmc_running(m)      
}

# function to estimate fitness from a single snapshot  2016-5-31 Thong Pham


estimate_alpha <- function(net, net_type = "directed") {
  stats       <- get_my_statistics(net, net_type = net_type)
  center_k    <- stats$center_k
  center_k[1] <- 1
  result      <- .estimate_alpha_core(stats$sum_m_k,stats$n_tk,stats$m_t, center_k);
  return(result)
}


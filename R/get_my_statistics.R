# function to summarize statistics from a growing network 
get_my_statistics <-
function(net, net_type = "directed"){
  net_object <- as.PAFit_net(net, type = "directed")
  return(get_statistics(net_object))
}

.onUnload <- function (libpath) {
  library.dynam.unload("mcPAFit", libpath)
}

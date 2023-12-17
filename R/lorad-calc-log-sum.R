#' Calculate a sum on log scale
#' 
#' Calculates the (natural) log of a sum without leaving the log scale by factoring 
#' out the largest element.
#' 
#' @param logx Numeric vector in which elements are on log scale
#' @return The log of the sum of the (exponentiated) elements supplied in logx

lorad_calc_log_sum <- function(logx) {
  n <- length(logx)
  max_logx <- max(logx)
  sum_of_terms <- 0
  
  for (i in 1:n) {
    sum_of_terms <- sum_of_terms + exp(logx[[i]] - max_logx)
  }
  
  log_sum <- log(sum_of_terms) + max_logx
  log_sum
}





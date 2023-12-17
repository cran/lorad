#' Transforms training sample using training sample means and standard deviations
#' 
#' @param standardinfo List containing the log Jacobian of the standardization 
#'     transformation, the inverse square root matrix, the column means, and 
#'     rmax (the radial distance representing the edge of the working parameter space)
#' @param y Data frame containing a column for each transformed model parameter 
#'     in the estimation sample, with last column being the log kernel values
#' @returns A new data frame consisting of the standardized estimation sample with
#'     log kernel in last column
#'
lorad_standardize_estimation_sample <- function(standardinfo, y) {
    # Remove log kernel column because it's not a parameter
    last_col_num <- ncol(y)
    x <- as.matrix(y[,-last_col_num])

    # Getting dimensions of matrix
    n <- nrow(x)
    p <- ncol(x)

    logJ <- standardinfo[[1]]
    invsqrts <- standardinfo[[2]]
    col_means <- standardinfo[[3]]

    inversesqrts_matrix <- matrix(invsqrts, nrow=p, ncol=p)

    # Construct mean matrix m as a matrix with n rows, each a p-dimensional vector of parameter means
    mean_vec <- matrix(col_means, nrow=p, ncol=1)
    m <- matrix(rep(mean_vec, n), byrow = TRUE, ncol = p)

    # Perform the standardization
    z <- (x - m)%*%inversesqrts_matrix

    # Convert z to dataframe
    zdf <- as.data.frame(z)

    # Adding log Jacobian to log kernel to complete transformation
    logkern <- names(y)[last_col_num]
    zdf[logkern] <- y[,last_col_num] + logJ

    # Return the new data frame
    zdf 
}





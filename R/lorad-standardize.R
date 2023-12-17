#' Transforms unconstrained parameters to have the same location and scale
#' 
#' Standardizes parameters that have already been transformed (if necessary) to have
#' unconstrained support. Standardization involves subtracting the sample mean 
#' and dividing by the sample standard deviation. Assumes that the log posterior kernel
#' (i.e. the log of the unnormalized posterior) is the last column in the supplied
#' data frame.
#' 
#' @param df Data frame containing a column for each model parameter sampled 
#'     and a final column of log posterior kernel values
#' @param coverage Fraction of the training sample used to compute working parameter space
#' @returns List containing the log-Jacobian of the standardization transformation, 
#'     the inverse square root matrix, a vector of column means, and rmax (radial 
#'     distance to furthest point in working parameter space)
#'
lorad_standardize <- function(df, coverage) {  
    # Get dimensions 
    n <- nrow(df)
    p <- ncol(df) - 1  # last column is not a parameter

    # Create matrix x from df by dropping the log_kernel column
    last_col_num <- ncol(df)
    x <- as.matrix(df[,-last_col_num])

    # Create matrix m having n rows each equal to the p-dimensional mean vector
    meanvect <- colMeans(x)
    m <- matrix(rep(colMeans(x), n), byrow = TRUE, ncol = p)

    # Calculate the variance-covariance matrix
    s <- stats::cov(x)

    # Getting Eigen vectors and values
    lamb <- eigen(s)$value
    vec <- eigen(s)$vector

    # Calculate square root of s, inverse square root of s, and log of the Jacobian
    # for the standardization transformation
    if (p == 1) {
        sqrts <- matrix(sqrt(lamb), ncol=1)
        invsqrts <- matrix(1/sqrt(lamb), ncol=1)
        logJ <- 0.5*log(lamb)
    }
    else {
        sqrts <- vec%*%diag(sqrt(lamb))%*%t(vec)
        invsqrts <- vec%*%diag(1/sqrt(lamb))%*%t(vec)
        logJ <- as.numeric(determinant(sqrts,logarithm = TRUE)$modulus)
    }

    # Perform the standardization
    z <- (x - m)%*%invsqrts

    # Create a vector of norms
    norms <- vector(mode="numeric", length=n) # numeric()
    for (i in 1:n) {
        # Obtain ith vector of params
        v <- z[i,]

        # Compute the radius (norm)
        r <- sqrt(sum(v^2))

        # Set ith element of norms vector to r
        norms[i] <- r
    }

    # Convert z to a data frame
    zdf <- as.data.frame(z)

    # Append norms vector to z
    zdf["norm"] <- norms

    # Append log kernel vector to z
    zdf["logkernel"] <- df[,last_col_num] + logJ

    # Sort z by norm
    zdfsorted <- zdf[order(zdf$norm),]
    
    # Retain only the first n*coverage elements of sorted data frame
    cutoff <- ceiling(n*coverage) + 1
    zdf_cropped <- zdfsorted[1:cutoff,]

    # The value rmax is the largest radius (norm) found in that portion of the
    # training sample retained
    rmax <- zdf_cropped$norm[cutoff]
    
    # Return important info
    list(logJ, invsqrts, meanvect, rmax)
}





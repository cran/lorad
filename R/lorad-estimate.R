#' Calculates the LoRaD estimate of the marginal likelihood
#'
#' Provided with a data frame containing sampled paraneter vectors and a dictionary
#' relating column names to parameter types, returns a named character vector containing
#' the following quantities:
#' * logML (the estimated log marginal likelihood)
#' * nsamples (number of samples)
#' * nparams (length of each parameter vector)
#' * training_frac (fraction of samples used for training)
#' * tsamples (number of samples used for training)
#' * esamples (number of sampled used for etimation)
#' * coverage (nominal fraction of the estimation sampled used)
#' * esamplesused (number of estimation samples actually used for estimation)
#' * realized_coverage (actual fraction of estimation sample used)
#' * rmax (lowest radial distance: defines boundary of working parameter space)
#' * log_delta (volume under the unnormalized posterior inside working parameter space)
#'
#' @param params Data frame in which rows are sample points and columns are parameters, 
#'     except that last column holds the log posterior kernel
#' @param colspec Named character vector associating column names in params with 
#'     column specifications
#' @param training_frac Number between 0 and 1 specifying the training fraction
#' @param training_mode One of random, left, or right, specifying how training 
#'     fraction is chosen
#' @param coverage Number between 0 and 1 specifying fraction of training sample 
#'     used to compute working parameter space
#' @returns Named character vector of length 11.
#' @export 
#' @examples
#' normals <- rnorm(1000000,0,10)
#' prob_normals <- dnorm(normals,0,10,log=TRUE) 
#' proportions <- rbeta(1000000,1,2)
#' prob_proportions <- dbeta(proportions,1,2,log=TRUE)
#' lengths <- rgamma(1000000, 10, 1)
#' prob_lengths <- dgamma(lengths,10,1,log=TRUE)
#' paramsdf <- data.frame(
#'     normals,prob_normals,
#'     proportions, prob_proportions,
#'     lengths, prob_lengths)
#' columnkey <- c(
#'     "normals"="unconstrained", 
#'     "prob_normals"="posterior", 
#'     "proportions"="proportion", 
#'     "prob_proportions"="posterior", 
#'     "lengths"="positive", 
#'     "prob_lengths"="posterior")
#' results <- lorad_estimate(paramsdf, columnkey, 0.5, 'random', 0.1)
#' lorad_summary(results)
#' 
lorad_estimate <- function(params, colspec, training_frac, training_mode, coverage) {
  if (length(params) == 0) {
    stop("params are missing")
  }
  if (length(colspec) == 0) {
    stop("colspec has 0 length")
  }
    # Transform any parameters that are constrained and consolidate log kernel components
    # into a single column named log_kernel that includes Jacobian terms for transformations
    transform_df <- lorad_transform(params, colspec)
    
    nsamples <- nrow(transform_df)
    nparams <- ncol(transform_df) - 1
    tmode <- tolower(training_mode)

    # Specified training fraction must be numeric
    if (!is.numeric(training_frac)) {
      stop(sprintf("training fraction must be numeric"))
    }
    # Specified coverage fraction must be numeric
    if (!is.numeric(coverage)) {
      stop(sprintf("coverage fraction must be numeric"))
    }
    # Specified training fraction must lie between 0 and 1
    if (training_frac <= 0 || training_frac >= 1.0) {
        stop(sprintf("training fraction must be between 0 and 1 but is %g",training_frac))
    }
    # Specified coverage fraction must lie between 0 and 1 
    if (coverage <= 0 || coverage >= 1.0) {
      stop(sprintf("coverage fraction must be between 0 and 1 but is %g",coverage))
    }
	# Determine sites included in training sample and place remainder in estimation sample
    y <- floor(training_frac*nsamples)
    tinclsamples <- "randomly chosen"
    einclsamples <- "complement of training samples set"
    if (tmode == "random") {
        x <- 1:nsamples
        z <- sample(x, y)
    }
    else if (tmode == "left") {
        z <- 1:y
        tinclsamples <- sprintf("samples 1 to %d", y)
        einclsamples <- sprintf("samples %d to %d", y+1, nsamples)
    }
    else if (tmode == "right") {
        z <- y:nsamples
        tinclsamples <- sprintf("samples %d to %d", y, nsamples)
        einclsamples <- sprintf("samples 1 to %d", y)
    }
    else {
        stop(sprintf("Unknown training mode %s", training_mode))
    }
    
    # Partition transformed samples into training and estimation samples
    training_df <- transform_df[z,]
    estimation_df <- transform_df[-z,] 
    tsamples <- nrow(training_df)
    esamples <- nrow(estimation_df)
    
    # Store the log kernel values in a vector
    last_col_num <- ncol(estimation_df)
    log_post_kern <- as.array(estimation_df[,last_col_num])

    # Compute mean vector and inverse square root of covariance matrix
    # needed for transforming the estimation sample. standard_info is 
    # a list with elements logJ, invsqrts, colmeans(x), and rmax
    standard_info <- lorad_standardize(training_df, coverage)
    logj <- standard_info[[1]]
    
    # Now standardize the estimation sample using mean vector and 
    # standard deviation matrix estimated from the training sample
    df <- lorad_standardize_estimation_sample(standard_info, estimation_df)
 
    # Extract posterior kernel
    last_col_num <- ncol(estimation_df)
    x <- as.matrix(df[,-last_col_num])
    log_post_kern <- df[,last_col_num]
  
    # The number of parameters
    p <- ncol(x)

    # The maximum radius of any point in the training sample
    rmax <- standard_info[[4]]
    if (is.nan(rmax)) {
      stop("there is not enough parameter variation to estimate the marginal likelihood")
    }

    # The variance is 1 for standard normal
    sigma_sqr <- 1
  
    # Compute Delta, which is the integral from 0 to rmax of the marginal
    # distribution of radial vector lengths from a p-dimensional multivariate
    # standard normal distribution
    s <- p/2.0
    t <- rmax^2/(2.0*sigma_sqr)
    log_delta <- log(stats::pgamma(t, shape=s, scale=1))
    if (is.nan(log_delta)) {
      stop("there is not enough parameter variation to estimate the marginal likelihood")
    }

    # Calculate the normalizing constant for reference function (multivariate standard normal)
    sigma <- sqrt(sigma_sqr)
    log_mvnorm_constant <- .5*p*log(2.0*pi)+1.0*p*log(sigma)
    
    # Initialize variables
    log_ratios <- numeric()
    nestimation <- nrow(estimation_df)
  
    # Calculate sum of ratios, Using multivariate standard normal reference function 
    j <- 0.0
    for (i in 1:nestimation) {
        # Calculate norm for the ith sample
        v <- x[i,]
        r <- sqrt(sum(v^2))
        
        # If norm is in working parameter space, use it
        if (r <= rmax) {
            j <- j+1 
            log_kernel <- log_post_kern[i]
            log_reference <- -.5*sigma_sqr*r^2 - log_mvnorm_constant
            log_ratio <- log_reference - log_kernel
            log_ratios[j] <- log_ratio
        }
    }
 
    esamplesused <- j
    realized_coverage <- j/nestimation
  
    if (length(log_ratios)==0){
        warning(sprintf("No estimation samples were within the working parameter space (rmax=%g)",rmax))
        stop()
    }
  
    log_sum_ratios <- lorad_calc_log_sum(log_ratios)

    #Calculate LoRaD estimate of maximum likelihood
    log_marginal_likelihood <- log_delta - (log_sum_ratios - log(nestimation))
    
    if (is.nan(log_marginal_likelihood)) {
      stop("there is not enough parameter variation to estimate the marginal likelihood")
    }
        
    results <- c(  
        "logML"=log_marginal_likelihood,
        "nsamples"=nsamples,
        "nparams"=nparams,
        "training_frac"=training_frac,
        "coverage"=coverage,
        "tsamples"=tsamples,
        "esamples"=esamples,
        "esamplesused"=esamplesused,
        "realized_coverage"=realized_coverage,
        "rmax"=rmax,
        "log_delta"=log_delta)
    results
}

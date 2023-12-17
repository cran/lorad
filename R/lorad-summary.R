#' Summarize output from [lorad_estimate()]
#' 
#' @param results Named character vector returned from [lorad_estimate()]
#' @returns String containing a summary of the supplied `results` object
#' @export 
#' @examples
#' normals <- rnorm(1000000,0,10)
#' prob_normals <- dnorm(normals,0,10,log=TRUE) 
#' paramsdf <- data.frame(normals,prob_normals)
#' columnkey <- c("normals"="unconstrained", "prob_normals"="posterior")
#' results <- lorad_estimate(paramsdf, columnkey, 0.5, 'left', 0.1)
#' lorad_summary(results)
lorad_summary <- function(results) {
    summary <- paste("This is lorad",utils::packageVersion("lorad"),"\n") 
    summary <- paste(summary, sprintf("   Parameter sample comprises %d sampled points\n", results[["nsamples"]]))
    if (results[["nparams"]] == 1) {
        summary <- paste(summary, "   Each sample point is a vector of 1 parameter value\n")
    } else {
        summary <- paste(summary, sprintf("   Each sample point is a vector of %d parameter values\n", results[["nparams"]]))
    }
    summary <- paste(summary, sprintf("   Training fraction is %g\n", results[["training_frac"]]), sep='\n')
    summary <- paste(summary, sprintf("   Coverage specified is %g\n", results[["coverage"]]))
    summary <- paste(summary, "\nPartitioning samples into training and estimation:\n")
    summary <- paste(summary, sprintf("   Sample size is %d\n",results[["nsamples"]]))
    summary <- paste(summary, sprintf("   Training sample size is %d\n", results[["tsamples"]]))
    summary <- paste(summary, sprintf("   Estimation sample size %d\n", results[["esamples"]]))
    summary <- paste(summary, "\nProcessing training sample...\n")
    summary <- paste(summary, sprintf("   Lowest radial distance is %.9f\n",results[["rmax"]]))
    summary <- paste(summary, sprintf("   Log Delta %.9f\n",results[["log_delta"]]))
    summary <- paste(summary, "\nProcessing estimation sample...\n")
    summary <- paste(summary, sprintf("   Number of samples used is %d\n",results[["esamplesused"]]))
    summary <- paste(summary, sprintf("   Nominal coverage specified is %f\n",results[["coverage"]]))
    summary <- paste(summary, sprintf("   Realized coverage is %f\n",results[["realized_coverage"]]))
    summary <- paste(summary, sprintf("   Log marginal likelihood is %f\n",results[["logML"]]))
    message(cat(summary))
}

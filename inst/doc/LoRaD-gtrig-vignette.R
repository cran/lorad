## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(lorad)

## -----------------------------------------------------------------------------
# Dimensions of gtrigsamples
dim(gtrigsamples)

# Column names of gtrigsamples
colnames(gtrigsamples)

## -----------------------------------------------------------------------------
# Create a named vector holding the column specifications
colspec <- c("Iteration"                  = "iteration", 
             "Posterior"                  = "posterior", 
             "Likelihood"                 = "ignore", 
             "Prior"                      = "ignore",
             "alpha"		                  = "positive",          
             "edge_length_proportions.1." = "simplex",	
             "edge_length_proportions.2." = "simplex",	
             "edge_length_proportions.3." = "simplex",	
             "edge_length_proportions.4." = "simplex",	
             "edge_length_proportions.5." = "simplex",	
             "edge_length_proportions.6." = "simplex",	
             "edge_length_proportions.7." = "simplexfinal",	
             "edgelens.1."                = "ignore",               
             "edgelens.2."                = "ignore",               
             "edgelens.3."                = "ignore",	              
             "edgelens.4."                = "ignore",	              
             "edgelens.5."                = "ignore",	              
             "edgelens.6."                = "ignore",             
             "edgelens.7."                = "ignore",	              
             "er.1."                      = "simplex",                    
             "er.2."                      = "simplex",	                    
             "er.3."                      = "simplex",                    
             "er.4."                      = "simplex",	                    
             "er.5."                      = "simplex",                    
             "er.6."                      = "simplexfinal",	                    
             "pi.1."                      = "simplex",                    
             "pi.2."                      = "simplex",	                    
             "pi.3."                      = "simplex",                    
             "pi.4."                      = "simplexfinal",	                    
             "pinvar"                     = "proportion",
             "site_rates.1."              = "ignore",
             "site_rates.2."              = "ignore",
             "site_rates.3."              = "ignore",
             "site_rates.4."              = "ignore",
             "tree_length"                = "positive")

## -----------------------------------------------------------------------------
results <- lorad_estimate(gtrigsamples, colspec, 0.5, "left", 0.1)
lorad_summary(results)


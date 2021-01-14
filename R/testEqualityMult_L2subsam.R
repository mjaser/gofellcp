

### Function to perform the test for equality for a given sample ###

testEqualityMult_L2subsam <- function(w, b, nb) {

  # Computes the L2 test statistic and
  # uses a subsampling bootstrap to get the empirical p-value.
  #
  # Input:
  #   w: numeric matrix or data frame with d columns (dimensions) and n rows (sample size).
  #   b: number of observations for the subsample.
  #   nb: number of bootstrap replications.
  #
  # Output:
  #   List containing the value of the test statistic and the empirical p-value of the test.
  #   $statistic: numeric vector of length one giving the value of the test statistic.
  #   $pvalue: numeric value giving the empirical p-value of the test.


  # Load functions
  source("compD.R")


  # Create list for output
  out <- list(statistic=NA, pvalue=NA)

  # Save sample size n, dimension d, number of pairs dd
  n <- dim(w)[1]
  d <- dim(w)[2]
  dd <- d * (d - 1) / 2

  # Compute pseudo-observations
  w_pseudo <- copula::pobs(w)

  # Compute D (vector of differences) for pseudo-observations of the sample
  D <- compD(w_pseudo)

  # Compute value of the test statistic (without variance)
  stat <- n * t(D) %*% (D)

  # Subsampling to compute empirical p-value
  # Vector to store bootstrap statistics
  stat_subsam <- rep(NA, nb)

  # Bootstrap: nb replications
  for (k in 1:nb) {

    # Generate bootstrap subsample from the data (without replacement)
    w_subsam <- w[sample(1:n, b, replace = FALSE), ]

    # Compute pseudo-observations
    w_subsam_pseudo <- copula::pobs(w_subsam)

    # Compute D for pseudo-observations of bootstrap subsample
    D_subsam <- compD(w_subsam_pseudo)

    # Compute value of the test statistic (without variance)
    stat_subsam[k] <- b * t(D_subsam - D) %*% (D_subsam - D)

  }

  # Compute the empirical p-value
  pval <- mean((stat_subsam > stat[1, 1]))

  # Store the values in the list for the output
  out$statistic <- stat
  out$pvalue <- pval

  # Return output list
  return(out)

}


# ### Examples
# library(copula)
#
# ### Gauss copula
#
# # Simulate copula data
# cop <- ellipCopula("normal", param=c(0.5, 0.6, 0.7) , dim=3, dispstr="un")
# w <- rCopula(500, cop)
#
# # Test function
# testEqualityMult_L2subsam(w, 50, 300)
#
#
# ### Frank copula
#
# # Simulate copula data
# cop <- frankCopula(3, dim=3)
# w <- rCopula(500, cop)
#
# # Test function
# testEqualityMult_L2subsam(w, 50, 300)

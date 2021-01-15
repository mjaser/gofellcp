

### Function to compute the test statistic for a given sample ###

testEqualityMult <- function(z) {

  # Computes the test statistic for a multivariate random sample
  #
  # Input:
  #   w: numeric matrix or data frame with d columns (dimensions) and n rows (sample size).
  #
  # Output:
  #   List containing the value of the test statistic and the p-value of the test.
  #   Under the null hypothesis, this test statistic has a standard normal distribution.
  #   $statistic: numeric vector of length one giving the value of the test statistic.
  #   $pvalue: numeric value giving the p-value of the test.


  # Create list for output
  out <- list(statistic=NA)

  # Save sample size n, and dimension d
  n <- dim(z)[1]
  d <- dim(z)[2]

  # Create two auxiliary matrices
  # to estimate Kendall's tau and Blomqvist's beta for each pair of coordinates
  FIRST  <- c()
  SECOND <- c()

  for(i in 1:(d-1)) {
    for(j in (i+1):d) {

      FIRST  <- cbind(FIRST, z[, i])
      SECOND <- cbind(SECOND, z[, j])

    }
  }

  # Number of pairs of coordinates
  dd <- dim(FIRST)[2]   # d(d-1)/2

  # Vector to store the differences
  # between the estimates ofKendall's tau and Blomqvist's beta
  D <- numeric(dd)

  # Auxiliary Matrix to store values needed
  # for the estimation of the covariance matrix Sigma
  A <- matrix(nrow=n, ncol=2*dd)

  # For each pair of coordinates
  for(i in 1:dd) {

    U <- FIRST[, i]
    V <- SECOND[, i]

    # Estimate Blomqvist's beta
    # Note that for copula data the median of the margins is 0.5
    hatBeta <- mean(sign(U - 0.5) * sign(V - 0.5))

    # Estimate Kendall's tau
    hatTau <- cor(x=U, y=V, method="kendall")

    # Compute the difference for this pair of coordinates
    D[i] <- hatBeta - hatTau


    # Vector of empirical copula
    Cn <- numeric(n)
    for (j in 1:n) {
      Cn[j] <- sum((U <= U[j]) & (V <= V[j])) / n
    }

    # Vector of h1
    h1 <-	1 - 2 * U - 2 * V + 4 * Cn


    # Store values for the estimation of Sigma
    A[, i] <- sign(U - 0.5) * sign(V - 0.5) - hatBeta
    A[, dd+i]   <- 2 * h1 - mean(2 * h1)

  }

  # Compute the estimator of the covariance matrix Sigma
  hatS  <- t(A) %*% A / n

  # Compute the Jacobian matrix of Phi for the delta method
  J <- cbind(diag(dd), -diag(dd))

  # Compute the estimated covariance matrix for the Wald-type statistic
  hatV <-  J %*% hatS %*% t(J)

  # Compute value of the Wald-type statistic
  stat <- n * t(D) %*% solve(hatV)  %*% (D)

  # Compute the p-value
  pval <- 1 - pchisq(stat, df=dd)

  # Store the values in the list for the output
  out$statistic <- stat
  out$pvalue <- pval

  return(out)

}


# ### Examples
# library(copula)
#
# ### Gauss copula
#
# # Simulate copula data
# cop <- ellipCopula("normal", param=c(0.5, 0.6, 0.7) , dim=3, dispstr="un")
# z <- rCopula(1000, cop)
#
# # Test function
# testEqualityMult(z)
#
#
# ### Frank copula
#
# # Simulate copula data
# cop <- frankCopula(3, dim=3)
# z <- rCopula(1000, cop)
#
# # Test function
# testEqualityMult(z)

# Functions needed for the testing procedure for ellipticity:
# - Test for symmetry
# - Test for radial symmetry
# - Test for equality of Kendall's tau and Blomqvist's beta
# plus Examples

# Load necessary packages
#library(copula)



##########################################################################################################################



### Test for symmetry
testSymmetry2D <- function(w) {

  # Tests a two dimensional random sample for symmetry.
  #
  # Input:
  #   w: numeric matrix or data frame with two columns.
  #
  # Output:
  #   List containing the value of the test statistic.
  #   Under the null hypothesis, this test statistic has a standard normal distribution.
  #   $statistic: numeric vector of length one giving the value of the test statistic.


  # Create list for output
  out <- list(statistic = NA)

  # Manipulate data in order to get two random samples w1 and w2
  n <- nrow(w)

  ind1 <- (w[, 1] - w[, 2]) <= rep(0, n)
  ind2 <- (w[, 1] - w[, 2]) >  rep(0, n)

  w1 <- w[ind1, ]
  w2 <- w[ind2, ]

  nn <- min(nrow(w1), nrow(w2))

  w1 <- w1[1:nn, ]
  w2 <- w2[1:nn, ]

  ind1 <- rbinom(nn, 1, 0.5)
  ind2 <- rbinom(nn, 1, 0.5)

  w1[ind1 == 0, ] <- w1[ind1 == 0, 2:1]
  w2[ind2 == 0, ] <- w2[ind2 == 0, 2:1]


  # Estimate Kendall's tau for both samples
  hatTau1 <- corKendall(w1)[1, 2]
  hatTau2 <- corKendall(w2)[1, 2]


  # Compute the test statistic
  T <- hatTau2 - hatTau1


  # Compute the estimated variance of T using
  # Var(hatTau)=Var(2h1(u,v)) and h1(u,v)=1-2u-2v+4Cn(u,v)

  # Original sample
  u <- w[, 1]
  v <- w[, 2]

  # Vector of empirical copula for original sample
  Cn <- numeric(n)
  for (i in 1:n) {
    Cn[i] <- sum((u <= u[i]) & (v <= v[i])) / n
  }

  # Vector of h1 for original sample
  h1 <-	1 - 2 * u - 2 * v + 4 * Cn

  # Compute the standard deviation of the test statistic
  sigma <- sqrt(2 * var(2 * h1))


  # Compute the transformed test statistic
  stat <-  sqrt(n/2) * T / sigma


  # Store the value of the transformed test statistic in the list for the output
  out$statistic <- stat

  return(out)

}



##########################################################################################################################



### Test for radial symmetry
testRadialSymmetry2D <- function(w) {

  # Tests a two dimensional random sample for radial symmetry.
  #
  # Input:
  #   w: numeric matrix or data frame with two columns.
  #
  # Output:
  #   List containing the value of the test statistic.
  #   Under the null hypothesis, this test statistic has a standard normal distribution.
  #   $statistic: numeric vector of length one giving the value of the test statistic.


  # Create list for output
  out <- list(statistic = NA)

  # Manipulate data in order to get two random samples
  n <- nrow(w)

  ind1 <- (w[, 1] + w[, 2]) <= rep(1, n)
  ind2 <- (w[, 1] + w[, 2]) >  rep(1, n)

  w1 <- w[ind1, ]
  w2 <- w[ind2, ]

  nn <- min(nrow(w1), nrow(w2))

  w1 <- w1[1:nn, ]
  w2 <- w2[1:nn, ]

  ind1 <- rbinom(nn, 1, 0.5)
  ind2 <- rbinom(nn, 1, 0.5)

  w1[ind1 == 0, ] <- 1 - w1[ind1 == 0, ]
  w2[ind2 == 0, ] <- 1 - w2[ind2 == 0, ]


  # Estimate Kendall's tau for both samples
  hatTau1 <- corKendall(w1)[1, 2]
  hatTau2 <- corKendall(w2)[1, 2]


  # Compute the test statistic
  T <- hatTau2 - hatTau1


  # Compute the estimated variance of T using
  # Var(hatTau)=Var(2h1(u,v)) and h1(u,v)=1-2u-2v+4Cn(u,v)

  # Original and reflected sample
  u <- w[, 1]
  v <- w[, 2]
  ur <- 1 - u
  vr <- 1 - v

  # Vector of empirical copula for original sample
  Cn <- numeric(n)
  for (i in 1:n) {
    Cn[i] <- sum((u <= u[i]) & (v <= v[i])) / n
  }

  # Vector of h1 for original sample
  h1 <-	1 - 2 * u - 2 * v + 4 * Cn


  # Vector of empirical copula for reflected sample
  Cnr <- numeric(n)
  for (i in 1:n) {
    Cnr[i] <- sum((u <= ur[i]) & (v <= vr[i])) / n
  }

  # Vector of h1 for reflected sample
  h1r <-	1 - 2 * ur - 2 * vr + 4 * Cnr


  # Compute the standard deviation of the test statistic
  sigma <- sqrt(2 * var(2 * h1))


  # Compute the transformed test statistic
  stat <-  sqrt(n/2) * T / sigma


  # Store the value of the transformed test statistic in the list for the output
  out$statistic <- stat

  return(out)

}



##########################################################################################################################



### Test for equality of Kendall's tau and Blomqvist's beta
testEquality2D <- function(w) {

  # Tests a two dimensional random sample for equality of Kendall's tau and Blomqvist's beta.
  # (intrinsic property of elliptical distributions)
  #
  # Input:
  #   w: numeric matrix or data frame with two columns.
  #
  # Output:
  #   List containing the value of the test statistic.
  #   Under the null hypothesis, this test statistic has a standard normal distribution.
  #   $statistic: numeric vector of length one giving the value of the test statistic.


  # Create list for output
  out <- list(statistic=NA)


  # Store given random sample in two vectors
  u <- w[, 1]
  v <- w[, 2]
  n <- length(u)


  # Estime Blomquvist's beta
  hatBeta <- mean(sign(u - 0.5) * sign(v - 0.5)) # for copula data the median of the margins is 0.5

  # Estimate Kendall's tau
  hatTau <- corKendall(w)[1, 2]

  # Compute the test statistic
  T <- hatBeta - hatTau


  # Compute the estimated variance of T using
  # Var(hatTau)=Var(2h1(u,v)) and h1(u,v)=1-2u-2v+4Cn(u,v)

  # Vector of empirical copula
  Cn <- numeric(n)
  for (i in 1:n) {
    Cn[i] <- sum((u <= u[i]) & (v <= v[i])) / n
  }

  # Vector of h1
  h1 <-	1 - 2 * u - 2 * v + 4 * Cn

  # Compute the estimated covariance matrix C
  A <- cbind(sign(u - 0.5) * sign(v - 0.5) - hatBeta, (2 * h1) - mean(2 * h1))
  C <- t(A) %*% A / n

  # Compute Df(t_0) for the delta method
  Df <- array(0, c(1, 2))

  Df[1, 1] <-  1
  Df[1, 2] <- -1

  # Compute the standard deviation of the test statistic
  std <- sqrt(Df %*% C %*% t(Df))


  # Compute the transformed test statistic
  stat <- sqrt(n) * T / std

  # Store the value of the transformed test statistic in the list for the output
  out$statistic <- stat

  return(out)

}




##########################################################################################################################



### Testing procedure for ellipticity
testProEllip <- function(u, alpha) {

  # Testing procedure for ellipticity:
  # 1) Performs tests for symmetry, radial symmetry and equality of Blomqvist's beta and Kendall's tau.
  # 2) Checks for Ellipticity with Bonferroni method.
  #
  # Input:
  #   data: dataframe with two columns containing the bivariate copula data set for the tests.
  #   alpha: significance level for the testing procedure.
  #
  # Output:
  #   List containing:
  #     reject: the rejection decisions for the testing procedure.
  #     pvalues: the p-values of the individual tests.
  #     copula_data: the copula data resulting from the transformation with the empirical cdf.


  # Create list for output
  out <- list(stat=NA, pvalues=NA, reject=NA)


  # Test for symmetry
  test_stat_sym <-  testSymmetry2D(u)$statistic
  p_sym  <- stats::pnorm(- abs(test_stat_sym)) + (1 - stats::pnorm(abs(test_stat_sym))) #1 - pchisq(test_stat^2, df=1)
  reject_sym <- as.numeric(p_sym < alpha/3)

  # Test for radial symmetry
  test_stat_radsym <-  testRadialSymmetry2D(u)$statistic
  p_radsym  <- stats::pnorm(- abs(test_stat_radsym)) + (1 - stats::pnorm(abs(test_stat_radsym)))
  reject_radsym <- as.numeric(p_radsym < alpha/3)

  # Test for equality
  test_stat_equal <-  testEquality2D(u)$statistic
  p_equal  <- stats::pnorm(- abs(test_stat_equal)) + (1 - stats::pnorm(abs(test_stat_equal)))
  reject_equal <- as.numeric(p_equal < alpha/3)

  # Testing procedure for ellipticity using Bonferroni
  reject_ellip <- NA
  if(reject_sym | reject_radsym | reject_equal) {
    # reject H_0
    reject_ellip <- 1
  } else {
    # accept H_0
    reject_ellip <- 0
  }

  out$stat <- data.frame(test_stat_sym, test_stat_radsym, test_stat_equal)
  out$pvalues <- data.frame(p_sym, p_radsym, p_equal)
  out$reject <- data.frame(reject_sym, reject_radsym, reject_equal, reject_ellip)

  return(out)

}



##########################################################################################################################



# ### Examples
#
# library(copula)
#
# # set significance level
# alpha <- 0.05
#
# # create copula object
# cop <- archmCopula("clayton", param=iTau(archmCopula("clayton", ), tau=0.5), dim=2)
# #cop <- ellipCopula("normal", param=iTau(ellipCopula(), tau=0.5), dim=2)
#
#
# # simulate copula data of given sample size
# w <- rCopula(1000, cop)
#
# # get pseudo observations
# #w <- pobs(w)
#
#
#
# ### perform the test for symmetry
# testSymmetry2D(w)
# test_stat_sym <-  testSymmetry2D(w)$statistic
#
# # pvalue
# (p_sym <- 1 - pchisq(test_stat_sym^2, df=1))
#
# # reject H_0^s?
# as.numeric(p_sym < alpha/3)
#
#
#
# ### perform the test for radial symmetry
# testRadialSymmetry2D(w)
# test_stat_radsym <-  testRadialSymmetry2D(w)$statistic
#
# # pvalue
# (p_radsym <- 1 - pchisq(test_stat_radsym^2, df=1))
#
# # reject H_0^r?
# as.numeric(p_radsym < alpha/3)
#
#
# ### perform the testing procedure for ellipticity
# testProEllip(w, alpha)


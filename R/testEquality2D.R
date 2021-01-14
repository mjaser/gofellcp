

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
  #hatTau <- cor(x=u, y=v, method="kendall")
  hatTau <- copula::corKendall(w)[1, 2]

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



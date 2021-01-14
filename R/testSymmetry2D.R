

testSymmetry2D <- function(w, method) {

  # Tests a two dimensional random sample for symmetry.
  #
  # Input:
  #   w: numeric matrix or data frame with two columns.
  #   method: character string specifying the method used for the data manipulation.
  #           One of "mixt" (mixture) or "refl" (reflection).
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

  if (method == "mixt") {

    ind1 <- stats::rbinom(nn, 1, 0.5)
    ind2 <- stats::rbinom(nn, 1, 0.5)

    w1[ind1 == 0, ] <- w1[ind1 == 0, 2:1]
    w2[ind2 == 0, ] <- w2[ind2 == 0, 2:1]

  } else if (method == "refl") {

    w1 <- rbind(w1, w1[, 2:1])
    w2 <- rbind(w2, w2[, 2:1])

  }


  # Estimate Kendall's tau for both samples
  #hatTau1 <- cor(x = w1[, 1], y = w1[, 2], method="kendall")
  #hatTau2 <- cor(x = w2[, 1], y = w2[, 2], method="kendall")
  hatTau1 <- copula::corKendall(w1)[1, 2]
  hatTau2 <- copula::corKendall(w2)[1, 2]


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
  # sigma does not depend on the method
  sigma <- sqrt(2 * stats::var(2 * h1))


  # Compute the transformed test statistic
  stat <-  sqrt(n/2) * T / sigma


  # Store the value of the transformed test statistic in the list for the output
  out$statistic <- stat

  return(out)

}



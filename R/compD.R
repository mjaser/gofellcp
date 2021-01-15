
compD <- function(w) {

  # Computes the statistic D_n (vector of beta-tau for all pairs of coordinates)
  # for the pseudo-observations of a given sample
  #
  # Input:
  #   w: numeric matrix or data frame with d columns (dimensions) and n rows (sample size).
  #
  # Output:
  #   numeric vector of length d*(d-1)/2
  #   containing the differences tau-beta for all pairs of coordinates.


  # Dimension
  d <- ncol(w)

  # Number of pairs of coordinates
  dd <- d * (d - 1) / 2

  # Create two auxiliary matrices
  # to estimate Kendall's tau and Blomqvist's beta for each pair of coordinates
  FIRST  <- c()
  SECOND <- c()

  for (i in 1:(d-1)) {
    for (j in (i+1):d) {

      FIRST  <- cbind(FIRST, w[, i])
      SECOND <- cbind(SECOND, w[, j])

    }
  }

  # Vector to store estimates of pairwise Blomqvist's beta
  hatBeta <- numeric(dd)

  # For each pair of coordinates
  for(i in 1:dd) {

    # Select data for i-th pair of coordinates
    U <- FIRST[, i]
    V <- SECOND[, i]

    # Estimate Blomqvist's beta with the empirical copula Cn
    hatBeta[i] <- 4 * copula::C.n(c(0.5, 0.5), cbind(U, V), smoothing = "none") - 1

  }

  # Estimate (fast) pairwise Kendall's tau
  hatTau <- copula::corKendall(w)

  # Compute the difference for this pair of coordinates
  D <- hatBeta - t(hatTau)[lower.tri(t(hatTau), diag = FALSE)]

  # Return vector of differences D
  return(D)

}




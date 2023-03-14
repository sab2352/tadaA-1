#' Bayesian FDR
#' @description Bayesian FDR control.
#' @param BF A sorted vector of BFs (in decreasing order).
#' @param pi0 The prior probability that the null model is true.
#' @param alpha The FDR target. Default is 0.05.
#' @return the q-value of each BF, and the number of findings with q below alpha.
#' @export
#' @examples NULL
Bayesian_FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  # [BF]: a sorted vector of BFs (in decreasing order)
  # [pi0]: the prior probability that the null model is true
  # [alpha]: the FDR target
  # [Return]: the q-value of each BF, and the number of findings with q below alpha.
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model

  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i

  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}


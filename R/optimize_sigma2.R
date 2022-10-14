
#' Optimize for sigma2
optimize_sigma2 <- function(R2, S2) {
  interval.max <- max((R2 - S2))
  if (interval.max <= 0)
    return(0)
  opt.res <- optimize(function(v) {
    sum(log(v + S2)) + sum(R2 / (v + S2))
  }, interval = c(0, interval.max), tol = sqrt(.Machine$double.eps))
  return(opt.res$minimum)
}

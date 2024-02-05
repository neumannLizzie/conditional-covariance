#' @export
covest <- function(x, z, znew, h, cest_old, sumK_old, mx) {
  .cov_est(x, z, znew, h, cest_old, sumK_old, mx)
}

#' @export
meanest <- function(x, z, znew, h, mean_old, sumK_old) {
  .mean_est(x, z, znew, h, mean_old, sumK_old)
}

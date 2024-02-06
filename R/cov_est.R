#' @name covest
#' @title Nadaraya-Watson type kernel estimation of the conditional covariance
#' @description Estimation of the conditional covariance.
#' @param x Numeric Matrix. An \R function of a discrete univariate Cumulative
#'  Distribution Function taking a numeric first argument and returning a
#'  numeric vector of the same length. For example, \code{\link{pbinom}} for
#'  the Binomial or \code{\link{ppois}} for the Poisson distribution. See
#'  \link{Distributions} for other standard distributions.
#' @param z Numeric. Enumerator of rational approximation of reference value
#' \code{k}.
#' @param znew Integer. Enumerator of rational approximation of control limit
#' value \code{h}.
#' @param h Numeric. Denominator of rational approximation of reference value.
#' @param cest_old Numeric Matrix. Randomization probability. If the CUSUM statistic is
#' equal to the threshold \code{h}, an control chart alarm is triggered with
#' probability \code{gamma}.
#' @param sumK_old Numeric. Switch for activating randomization in order to allow
#' continuous ARL control.
#' @param mx Numeric Matrix. Head start value as integer multiple of \code{1/m}.
#' Should be an element of \code{0:hm}.
#' @return Returns a single value which is the Average Run Length.
#' @details TODO: Equations how the CUSUMs are implemented
#' @author Lizzie Neumann
#' @references Neumann et al. (2024).
#'  Confounder-adjusted Covariances of System Outputs and Applications to Structural Health Monitoring.
#'  \emph{PREPRINT}, \strong{XX}, pp. 1--26.
#'
#'  Yin et al. (2010).
#'  Nonparametric Covariance Model. \emph{Statistica Sinica} \strong{20}, 469-–79.
#' @examples
#' \dontrun{
#' A   <- 500
#' mu0 <- 1
#' m   <- 20
#' km1 <- 1.3*m
#' km2 <- 2.2*m
#' CDF <- function(q) ppois(q, lambda=mu0)
#' ## calibrate Bi-CUSUM
#' res1 <- bicusum_discrete_crit(CDF, mu0, A, km1, km2, m)
#' ## ARLs of single CUSUM
#' C1 <- cusum_discrete_arl(CDF, km1, res1$hm1, m, gamma=res1$g1, rando=TRUE)
#' C2 <- cusum_discrete_arl(CDF, km2, res1$hm2, m, gamma=res1$g2, rando=TRUE)
#' ## ARL of Bi-CUSUM with randomization
#' CC <- bicusum_discrete_arl(CDF, km1, km2, res1$hm1, res1$hm2, m, rando=TRUE, g1=res1$g1, g2=res1$g2)
#' cbind(round(cbind(C1, C2, CC), 3))
#' ## ARL of Bi-CUSUM without randomization
#' bicusum_discrete_arl(CDF, km1, km2, res1$hm1, res1$hm2, m)
#' }
#' @export
covest <- function(x, z, znew, h, cest_old, sumK_old, mx) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_matrix(x)
  checkmate::assert_vector(z)
  checkmate::assert_numeric(znew, len = 1, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, add = arg_checks)
  checkmate::assert_matrix(cest_old)
  checkmate::assert_numeric(sumK_old, len = 1, add = arg_checks)
  checkmate::assert_matrix(mx)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .cov_est(x, z, znew, h, cest_old, sumK_old, mx)
}

#' @name meanest
#' @title Zero-state ARL of discrete CUSUM charts
#' @description Compute the ARL of CUSUM.
#' @param x Numreic Matrix An \R function of a discrete univariate Cumulative
#'  Distribution Function taking a numeric first argument and returning a
#'  numeric vector of the same length. For example, \code{\link{pbinom}} for
#'  the Binomial or \code{\link{ppois}} for the Poisson distribution. See
#'  \link{Distributions} for other standard distributions.
#' @param z Numeric Vector. Enumerator of rational approximation of reference value
#' \code{k}.
#' @param znew Numeric. Enumerator of rational approximation of control limit
#' value \code{h}.
#' @param h Numeric. Denominator of rational approximation of reference value.
#' @param mean_old Numeric. Randomization probability. If the CUSUM statistic is
#' equal to the threshold \code{h}, an control chart alarm is triggered with
#' probability \code{gamma}.
#' @param sumK_old Boolean. Switch for activating randomization in order to allow
#' continuous ARL control.
#' @param mx Integer. Head start value as integer multiple of \code{1/m}.
#' Should be an element of \code{0:hm}.
#' @return Returns a single value which is the Average Run Length.
#' @details TODO: Equations how the CUSUMs are implemented
#' @author Lizzie Neumann
#' @references Brook and Evans (1972).
#'  An approach to the probability distribution of CUSUM run length.
#'  \emph{Biometrika}, \strong{59}, pp. 539--549
#' @examples
#' \dontrun{
#' A   <- 500
#' mu0 <- 1
#' m   <- 20
#' km1 <- 1.3*m
#' km2 <- 2.2*m
#' CDF <- function(q) ppois(q, lambda=mu0)
#' ## calibrate Bi-CUSUM
#' res1 <- bicusum_discrete_crit(CDF, mu0, A, km1, km2, m)
#' ## ARLs of single CUSUM
#' C1 <- cusum_discrete_arl(CDF, km1, res1$hm1, m, gamma=res1$g1, rando=TRUE)
#' C2 <- cusum_discrete_arl(CDF, km2, res1$hm2, m, gamma=res1$g2, rando=TRUE)
#' ## ARL of Bi-CUSUM with randomization
#' CC <- bicusum_discrete_arl(CDF, km1, km2, res1$hm1, res1$hm2, m, rando=TRUE, g1=res1$g1, g2=res1$g2)
#' cbind(round(cbind(C1, C2, CC), 3))
#' ## ARL of Bi-CUSUM without randomization
#' bicusum_discrete_arl(CDF, km1, km2, res1$hm1, res1$hm2, m)
#' }
#' @export
meanest <- function(x, z, znew, h, mean_old, sumK_old) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_matrix(x)
  checkmate::assert_vector(z)
  checkmate::assert_numeric(znew, len = 1, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, add = arg_checks)
  checkmate::assert_vector(mean_old)
  checkmate::assert_numeric(sumK_old, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .mean_est(x, z, znew, h, mean_old, sumK_old)
}

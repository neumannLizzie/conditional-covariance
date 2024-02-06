#' @name covest
#' @title Nadaraya-Watson type kernel estimation of the conditional covariance
#' @description Estimation of the conditional covariance using a Nadaraya-Watson 
#' type kernel estimatior.
#' @param x Numeric Matrix. An nxp matrix containing the data (of, e.g., sensor 
#' measurements or features) where the columns refer to p different sensors and 
#' n is the number of observations. 
#' @param z Numeric Vector. Containing the measured confounder (e.g., temperature).
#' @param znew Integer. The value for which the conditional covariance is to be estimated.
#' @param h Numeric. The bandwidths, i.e. the smoothing parameter.
#' @param cest_old Numeric Matrix. Randomization probability. If the CUSUM statistic is
#' equal to the threshold \code{h}, an control chart alarm is triggered with
#' probability \code{gamma}.
#' @param sumK_old Numeric. Switch for activating randomization in order to allow
#' continuous ARL control.
#' @param mx Numeric Matrix. Head start value as integer multiple of \code{1/m}.
#' Should be an element of \code{0:hm}.
#' @return Returns a matrix with p rows and p+1 columns. 
#' The first row of the p+1 column contains the sum of the kernel (sumK).
#' The first pxp matrix is the covariance matrix at value znew.
#' @details TODO: Equations how the cond cov is implemented
#' @author Lizzie Neumann
#' @references Neumann et al. (2024).
#'  Confounder-adjusted Covariances of System Outputs and Applications to Structural Health Monitoring.
#'  \emph{PREPRINT}, \strong{XX}, pp. 1--26.
#'
#'  Yin et al. (2010).
#'  Nonparametric Covariance Model. \emph{Statistica Sinica} \strong{20}, 469-–79.
#' @examples
#' \dontrun{
#' x <- cbind(c(1,3,5),c(3,2,6),c(4,1,5))
#' z <- c(5.1,3.1,4.3)
#' znew   <- 3.2
#' h <- 1.3
#' cest_old <- matrix(0,ncol=3,nrow=3)
#' sumK_old <- 0
#' mx <- cbind(c(1.2,2.9,3.8),c(3.2,2.2,1.1))
#' 
#' covest::covest(x=x,z=z,znew=znew,h=h,cest_old=cest_old, sumK_old=sumK_old, mx=mx)
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
#' @title Nadaraya-Watson type kernel estimation of the conditional mean
#' @description Estimation of the conditional mean using a Nadaraya-Watson type kernel estimatior.
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
#' @references Neumann et al. (2024).
#'  Confounder-adjusted Covariances of System Outputs and Applications to Structural Health Monitoring.
#'  \emph{PREPRINT}, \strong{XX}, pp. 1--26.
#'
#'  Yin et al. (2010).
#'  Nonparametric Covariance Model. \emph{Statistica Sinica} \strong{20}, 469-–79.
#' @examples
#' \dontrun{
#' x <- cbind(c(1,3,5),c(3,2,6),c(4,1,5))
#' z <- c(5.1,3.1,4.3)
#' znew   <- 3.2
#' h <- 1.3
#' mean_old <- rep(0,3)
#' sumK_old <- 0
#' 
#' covest::meanest(x=x,z=z,znew=znew,h=h,mean_old=cest_mean_oldold, sumK_old=sumK_old)
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

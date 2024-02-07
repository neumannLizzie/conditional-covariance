#' @name covest
#' @title Nadaraya-Watson type kernel estimation of the conditional covariance
#' @description Estimation of the conditional covariance using a Nadaraya-Watson
#' type kernel estimatior.
#' @param x Numeric Matrix. An nxp matrix containing the data (of, e.g., sensor
#' measurements or features) where the columns refer to p different sensors and
#' n is the number of observations.
#' @param z Numeric Vector. Containing the measured confounder (e.g., temperature).
#' @param znew Numeric. The value for which the conditional covariance is to be estimated.
#' @param h Numeric. The bandwidths, i.e. the smoothing parameter.
#' @param cest_old Numeric Matrix. The conditional covariance at value znew.
#' Only not zero if the data set is split and the conditional covariance is estimated
#' and updated piece by piece or new data is added.
#' @param sumK_old Numeric. Sum of kernels. Only not zero if the data set is split
#' and the conditional covariance is estimated and updated piece by piece or new data is added.
#' @param mx Numeric Matrix. The estimated conditional mean, a nxp matrix where
#' each row contains the conditional mean for each value of z.
#' @return Returns a matrix with p rows and p+1 columns.
#' The pxp matrix is the covariance matrix at value znew.
#' The first element of the p+1 column contains the sum of the kernel (sumK).
#' @details The used kernel is the Gaussian kernel which equals a normal density with
#' mean 0 and chosen bandwidth h.
#' K = dnorm(z - znew, 0, h, false)
#'
#' The conditional covariance is estimated as follows
#' cest = (K*t(x - mx)*(x - mx)  + sumK_old*cest_old) / (sum(K) + sumK_old) where
#' mx is the conditional mean estimated using, e.g., meanest.
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
#' @param x Numeric Matrix. An nxp matrix containing the data (of, e.g., sensor
#' measurements or features) where the columns refer to p different sensors and
#' n is the number of observations.
#' @param z Numeric Vector. Containing the measured confounder (e.g., temperature).
#' @param znew Integer. The value for which the conditional covariance is to be estimated.
#' @param h Numeric. The bandwidths, i.e. the smoothing parameter.
#' @param mean_old Numeric Vector. A 1xp vector, containing the conditional mean at value znew.
#' Only not zero if the data set is split and the conditional mean is estimated
#' and updated piece by piece or new data is added.
#' @param sumK_old Numeric. The sum of kernels. Only not zero if the data set is
#' split and the conditional mean is estimated and updated piece by piece or new data is added.
#' @return Returns a matrix where the first column contains the conditional mean
#' at value znew and the first element in the second column is the sum of kernels (sumK).
#' @details The used kernel is the Gaussian kernel which equals a normal density with
#' mean 0 and chosen bandwidth h.
#' K = dnorm(z - znew, 0, h, false)
#'
#' The conditional mean is estimated as follows
#' mest = (K*t(x) + sumK_old*mean_old) / (sum(K) + sumK_old).
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

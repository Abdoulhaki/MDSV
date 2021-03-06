#' CAC 40
#'
#' Data from Oxford-Man institute of the CAC 40 from date "2000-01-03" to date "2021-02-24". 
#' Using close to close prices, log-retuns and realized variances are computed using formulas 
#' \eqn{r_t = log(P_{t})-log(P_{t-1})} and \eqn{rv_t = \sum_i r_{t,i}^2} where \eqn{r_{t,i}} stands 
#' for the intra-day log-returns of the date t and \eqn{P_t} stands for the price at day t.
#'
#' @docType data
#'
#' @usage data(cac40)
#'
#' @format Matrix of two columns where the first column is the log-returns and the second the realized variances.
#'
#' @keywords datasets
#'
#' @source \href{https://realized.oxford-man.ox.ac.uk/}{Oxford-Man Realized Library}
#'
#' @examples
#' \dontrun{
#' data(cac40)
#' }

"cac40"
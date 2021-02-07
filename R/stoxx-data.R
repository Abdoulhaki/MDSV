#' Arabidopsis QTL data on gravitropism
#'
#' Data from a QTL experiment on gravitropism in
#' Arabidopsis, with data on 162 recombinant inbred lines (Ler x
#' Cvi). The outcome is the root tip angle (in degrees) at two-minute
#' increments over eight hours.
#'
#' @docType data
#'
#' @usage data(sp500)
#'
#' @format Matrix of two columns where the first column is the log-returns and the second the realized variances.
#'
#' @keywords datasets
#'
#' @source \href{https://realized.oxford-man.ox.ac.uk/}{Oxford-Man Realized Library}
#'
#' @examples
#' @dontrun{
#' data(sp500)
#' dates <- names(sp500)[1:nrow(sp500)]
#' }

"sp500"
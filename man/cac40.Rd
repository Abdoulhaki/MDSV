% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cac40-data.R
\docType{data}
\name{cac40}
\alias{cac40}
\title{CAC 40}
\format{
Matrix of two columns where the first column is the log-returns and the second the realized variances.
}
\source{
\href{https://realized.oxford-man.ox.ac.uk/}{Oxford-Man Realized Library}
}
\usage{
data(cac40)
}
\description{
Data from Oxford-Man institute of the CAC 40 from date "2000-01-03" to date "2021-02-24". 
Using close to close prices, log-retuns and realized variances are computed using formulas 
\eqn{r_t = log(P_{t})-log(P_{t-1})} and \eqn{rv_t = \sum_i r_{t,i}^2} where \eqn{r_{t,i}} stands 
for the intra-day log-returns of the date t and \eqn{P_t} stands for the price at day t.
}
\examples{
\dontrun{
data(cac40)
}
}
\keyword{datasets}

#' @title MDSV Filtering
#' @description Method for filtering the MDSV model on log-retruns and realized variances (uniquely or jointly).
#' @param N An integer designing the number of components for the MDSV process
#' @param K An integer designing the number of states of each MDSV process component
#' @param data A univariate or bivariate data matrix. Can only be a matrix of 1 or 2 columns. If data has 2 columns, the first one has to be the log-returns and the second the realized variances.
#' @param para A vector of parameters use for the MDSV filtering on \code{data}. For more informations see Details. 
#' @param ModelType An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.
#' @param LEVIER if \code{TRUE}, estime the MDSV model with leverage.
#' 
#' @return A list consisting of: 
#' \itemize{
#'     \item ModelType : type of model to be filtered.
#'     \item LEVIER : wheter the filter take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item data : data use for the filtering.
#'     \item dates : vector or names of data designing the dates.
#'     \item estimates : input parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item Levier : numeric vector representing the leverage effect at each date. \code{Levier} is 1 when no leverage is detected
#'     \item filtred_proba : matrix containing the filtred probabilities \eqn{P(C_t=c_i | x_1,\dots, x_t)} of the Markov Chain.
#'     \item smoothed_proba : matrix containing the smoothed probabilities \eqn{P(C_t=c_i | x_1,\dots, x_T)} of the Markov Chain.
#'     \item Marg_loglik : marginal log-likelihood corresponding to the log-likelihood of log-returns. This is only return when \eqn{ModelType = 2}.
#'     \item VaR95 : 5\% Value-at-Risk compute empirically.
#'     \item VaR99 : 1\% Value-at-Risk compute empirically.
#' }
#' 
#' @details 
#' The MDSVfilter help to simply filter a set of data with a predefined set of parameters. Like in \code{\link{MDSVfit}}, the 
#' likelihood calculation is performed in \code{C++} through the \pkg{Rcpp} package. 
#' According to Augustyniak et al., (2021), the para consist on a vector of :
#' \itemize{
#'     \item{\eqn{\omega}}{ probability of success of the stationnary distribution of the Markov chain. (\eqn{0 < \omega < 1}).}
#'     \item{\eqn{a}}{ highest persistance of the component of the MDSV process. (\eqn{0 < a < 1}).}
#'     \item{\eqn{b}}{ decroissance rate of the persistances. (\eqn{1 < b}).}
#'     \item{\eqn{\sigma^2}}{ unconditionnal variance of the MDSV process.}
#'     \item{\eqn{\nu_0}}{ states defined parameter. (\eqn{0 < \nu_0 < 1}).}
#'     \item{\eqn{\gamma}}{ parameter of the realized variances innovation. This parameter is required only if \eqn{ModelType = 1} or \eqn{ModelType = 2}. (\eqn{0 < \gamma}).}
#'     \item{\eqn{\xi}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
#'     \item{\eqn{\phi}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
#'     \item{\eqn{\delta1}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
#'     \item{\eqn{\delta2}}{ parameter of realized variances equation in joint estimation. This parameter is required only if \eqn{ModelType = 2}.}
#'     \item{\eqn{l}}{ leverage effect parameter. This parameter is required only if \eqn{LEVIER = TRUE}. (\eqn{0 < l}).}
#'     \item{\eqn{\theta}}{ leverage effect parameter. This parameter is required only if \eqn{LEVIER = TRUE}. (\eqn{0 < \theta < 1}).}
#' }
#' The leverage effect is taken into account according to the FHMV model (see Augustyniak et al., 2019). While filtering an
#' univariate realized variances data, log-returns are required to add leverage effect. 
#' AIC and BIC are computed using the formulas : 
#' \itemize{
#' \item{AIC : }{\eqn{L - k}}
#' \item{BIC : }{\eqn{L - (k/2)*log(n)}}
#' }
#' where \eqn{L} is the log-likelihood, \eqn{k} is the number of parameters and \eqn{n} the number of observations in the dataset.
#' The \link[base]{class} of the output of this function is \code{\link{MDSVfilter}}. This class has a \link[base]{summary}, 
#' \link[base]{print} and \link[base]{plot} \link[utils]{methods} to summarize, print and plot the results. See 
#' \code{\link{summary.MDSVfilter}}, \code{\link{print.MDSVfilter}} and \code{\link{plot.MDSVfilter}} for more details.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' @references 
#' Augustyniak, M., Dufays, A., & Maoude, K.H.A. (2021). Multifractal Discrete Stochastic Volatility.
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, bootstrap forecasting \code{\link{MDSVboot}}, simulation \code{\link{MDSVsim}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @examples 
#' \dontrun{
#' # MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
#' data(sp500)        # Data loading
#' N         <- 2     # Number of components
#' K         <- 3     # Number of states
#' ModelType <- 0     # Univariate log-returns
#' LEVIER    <- FALSE # No leverage effect
#' 
#' # Model estimation
#' out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = ModelType, LEVIER = LEVIER)
#' # Model filtering
#' para      <-out_fit$estimates # parameter
#' out_filter<- MDSVfilter(K = K, N = N, data = sp500, para = para, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out_filter)
#' # Plot
#' plot(out_filter)
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' data(nasdaq)      # Data loading
#' N         <- 3    # Number of components
#' K         <- 3    # Number of states
#' ModelType <- 2    # Joint log-returns and realized variances
#' LEVIER    <- TRUE # No leverage effect
#' 
#' para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72, 
#'               xi = -0.5, varphi = 0.93, delta1 = 0.93, delta2 = 0.04, shape = 2.10,
#'               l = 0.78, theta = 0.876)
#' # Model filtering
#' out       <- MDSVfilter(K = K, N = N,data = nasdaq, para = para, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out)
#' # Plot
#' plot(out)
#' }
#' 
#' @import Rcpp
#' @export
MDSVfilter<-function(N,K,data,para,ModelType=0,LEVIER=FALSE){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVfilter(): input N and K must be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVfilter(): input N and K must be integer!")
  }else if(K<2){
    stop("MDSVfit(): input K must be greater than one!")
  }else if(N<1){
    stop("MDSVfit(): input N must be positive!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVfilter(): input ModelType must be numeric!")
  }else if(!(ModelType %in% c(0,1,2))){
    stop("MDSVfilter(): input ModelType must be 0, 1 or 2!")
  }
  
  if(!is.logical(LEVIER)) {
    stop("MDSVfilter(): input LEVIER must be logical!")
  }
  
  if ( (!is.numeric(data)) || (!is.matrix(data))  ) {
    stop("MDSVfilter(): input data must be numeric matrix!")
  }
  
  if(!is.numeric(para)) {
    stop("MDSVfilter(): input para must be numeric!")
  }else if(!is.vector(para)) {
    stop("MDSVfilter(): input para must be vector!")
  }else if(((!LEVIER) & (ModelType==0) & !(length(para)==5)) ||
       ((!LEVIER) & (ModelType==1) & !(length(para)==6)) ||
       ((!LEVIER) & (ModelType==2) & !(length(para)==10)) ||
       ((LEVIER) & (ModelType==0) & !(length(para)==7)) ||
       ((LEVIER) & (ModelType==1) & !(length(para)==8)) ||
       ((LEVIER) & (ModelType==2) & !(length(para)==12))){
    stop("MDSVfilter(): incorrect input para!")
  }
  
  if((para[1]>1) || (para[1]<0)) {
    stop("MDSVfilter(): input para[omega] must be between 0 and 1!")
  }else if((para[2]>1) || (para[2]<0)) {
    stop("MDSVfilter(): input para[a] must be between 0 and 1!")
  }else if((para[3]<=1)) {
    stop("MDSVfilter(): input para[b] must be greater than 1!")
  }else if((para[4]<=0)) {
    stop("MDSVfilter(): input para[sigma] must be greater than 0!")
  }else if((para[5]>1) || (para[5]<0)) {
    stop("MDSVfilter(): input para[v0] must be between 0 and 1!")
  }else if((ModelType==1) & (para[6]<=0)) {
    stop("MDSVfilter(): input para[shape] must be greater than 0!")
  }else if((ModelType==2) & (para[10]<=0)){
    stop("MDSVfilter(): input para[shape] must be greater than 0!")
  }else if(LEVIER){
    if(ModelType==0){
      if(para[6]<0){
        stop("MDSVfilter(): input para[l] must be greater than 0!")
      }else if((para[7]>1) || (para[7]<0)){
        stop("MDSVfilter(): input para[theta_l] must be between 0 and 1!")
      }
    }else if(ModelType==1){
      if(para[7]<=0){
        stop("MDSVfilter(): input para[l] must be greater than 0!")
      }else if((para[8]>1) || (para[8]<0)){
        stop("MDSVfilter(): input para[theta_l] must be between 0 and 1!")
      }
    }else if(ModelType==2){
      if(para[11]<=0){
        stop("MDSVfilter(): input para[l] must be greater than 0!")
      }else if((para[12]>1) || (para[12]<0)){
        stop("MDSVfilter(): input para[theta_l] must be between 0 and 1!")
      }
    }
  }
  
  k <- ncol(data)
  T <- nrow(data)
  
  if(!is.null(names(data))) {
    dates <- as.Date(names(data)[1:T])
  }else {
    dates<- 1:T
  }
  
  if((k==1) & (!(ModelType == 0) & (!((ModelType == 1) & (LEVIER == FALSE))))){
    stop("MDSVfilter(): improperly data dimensioned matrices!")
  }
  
  if((k==1) & (ModelType == 1) & sum(data<0)>0 ){
    stop("MDSVfilter(): data must be positive!")
  } 
  if((k==1) & (ModelType == 1)) data <- matrix(c(rep(1,nrow(data)),data),nrow(data),2)
  
  if(k==2) if(sum(data[,2]<0)>0 ){
    stop("MDSVfilter(): data second colomn must be positive!")
  }
  
  l<-logLik2(ech=data, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N)
  if(!(ModelType==1)){
    pi_0 <- l$w_hat
    sig<- volatilityVector(para=para,N=N,K=K)
    if(LEVIER){
      Levier<-levierVolatility(ech=data[((T-200):T),1],para=para,Model_type=ModelType)$`levier`
      sig<-sig*Levier
    }
  }
  
  ### Results
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)       vars <- c(vars,"l","theta")
  
  names(para)<-vars
  if(N==1) para<-para[(vars[!(vars=='b')])]
  
  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  
  out<-list(ModelType      = Model_type, 
            LEVIER         = LEVIER, 
            N              = N, 
            K              = K,
            data           = data,
            dates          = dates,
            estimates      = para, 
            LogLikelihood  = -l$loglik,
            AIC            = -l$loglik-length(para), 
            BIC            = -l$loglik-0.5*length(para)*log(T),
            Levier         = l$Levier,
            filtred_proba  = l$filtred_proba,
            smoothed_proba = l$smoothed_proba)
  
  if(ModelType==2) out<-c(out,list(Marg_loglik = l$Marg_loglik))
  
  if(!(ModelType==1)){
    out<-c(out,list(VaR95 = qmist2n(0.05,sigma=sig,p=pi_0),
                    VaR99 = qmist2n(0.01,sigma=sig,p=pi_0)))
  }
  
  class(out) <- "MDSVfilter"
  
  return(out)
}

pmist2n <- function(y,sigma,p){
  return(sum(p*pnorm(y, 0, sqrt(sigma))))
}#pmist2n

qmist2n <- function(q,sigma,p){
  # the min minmax below is computed to supply a range to the solver
  # the solution must be between the min and max
  # quantile of the mixed distributions
  minmax <- range(qnorm(q,0,sqrt(sigma)))
  uniroot(function(x) pmist2n(x,sigma=sigma,p=p)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}

#' @title Summarize and print MDSV Filtering
#' @description Summary and print methods for the class \link{MDSVfilter} as returned by the function \code{\link{MDSVfilter}}.
#' @param object An object of class \link{MDSVfilter}, output of the function \code{\link{MDSVfilter}}.
#' @param x An object of class \link{summary.MDSVfilter}, output of the function \code{\link{summary.MDSVfilter}}
#' or class \link{MDSVfilter} of the function \code{\link{MDSVfilter}}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be filtered.
#'     \item LEVIER : wheter the filter take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item data : data use for the filtering.
#'     \item dates : vector or names of data designing the dates.
#'     \item estimates : input parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item Levier : numeric vector representing the leverage effect at each date. \code{Levier} is 1 when no leverage is detected
#'     \item filtred_proba : matrix containing the filtred probabilities \eqn{\mathbb{P}(C_t=c_i\mid x_1,\dots,\x_t)} of the Markov Chain.
#'     \item smoothed_proba : matrix containing the smoothed probabilities \eqn{\mathbb{P}(C_t=c_i\mid x_1,\dots,\x_T)} of the Markov Chain.
#'     \item Marg_loglik : marginal log-likelihood corresponding to the log-likelihood of log-returns. This is only return when \eqn{ModelType = 2}.
#'     \item VaR95 : 5% Value-at-Risk compute empirically.
#'     \item VaR99 : 1% Value-at-Risk compute empirically.
#' }
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @export

"summary.MDSVfilter" <- function(object, ...){
  stopifnot(class(object) == "MDSVfilter")
  
  out<-c(object,...)
  
  class(out) <- "summary.MDSVfilter"
  return(out)
}


#' @rdname summary.MDSVfilter
#' @export

"print.summary.MDSVfilter" <- function(x, ...){
  stopifnot(class(x) == "summary.MDSVfilter")
  
  cat("=================================================\n")
  cat(paste0("================  MDSV Filtering ================\n"))
  cat("=================================================\n\n")
  # cat("Conditional Variance Dynamique \n")
  # cat("------------------------------------------------- \n")
  cat(paste0("Model   : MDSV(",x$N,",",x$K,")\n"))
  cat(paste0("Data    : ",x$ModelType,"\n"))
  cat(paste0("Leverage: ",x$LEVIER,"\n\n"))
  
  cat("Optimal Parameters \n")
  cat("------------------------------------------------- \n")
  cat(paste(paste(names(x$estimates),round(x$estimates,6),sep=" \t: "),collapse = "\n"))
  cat("\n\n")
  cat(paste0("LogLikelihood \t: ", round(x$LogLikelihood,2),"\n"))
  if(x$ModelType == "Joint log-return and realized variances"){
    cat(paste0("Marginal LogLikelihood : ", round(x$Marg_loglik,2),"\n\n"))
  } else{
    cat("\n")
  }
    
  cat("Information Criteria \n")
  cat("------------------------------------------------- \n")
  cat(paste0("AIC \t: ", round(x$AIC,2),"\n"))
  cat(paste0("BIC \t: ", round(x$BIC,2),"\n\n"))
  
  if(!(x$ModelType == "Univariate realized variances")){
    cat("Value at Risk \n")
    cat("------------------------------------------------- \n")
    cat(paste0("95%  \t: ", round(x$VaR95,6),"\n"))
    cat(paste0("99%  \t: ", round(x$VaR99,6),"\n"))
  }
  
  invisible(x)
}

#' @rdname summary.MDSVfilter
#' @export
"print.MDSVfilter" <- function(x, ...) {
  stopifnot(class(x) == "MDSVfilter")
  print(summary(x, ...))
}

#' @title Plot MDSV Filtering
#' @description Plot methods for the class \link{MDSVfilter} as returned by the function \code{\link{MDSVfilter}}.
#' @param x An object of class \link{MDSVfilter}, output of the function \code{\link{MDSVfilter}}
#' or class \link{plot.MDSVfilter} of the function \code{\link{plot.MDSVfilter}}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item V_t : smoothed volatilities taking leverage effect into accound when existing.
#'     \item data : data use for the filtering.
#'     \item dates : vector or names of data designing the dates.
#'     \item ModelType : type of model to be filtered.
#'     \item ... : further arguments passed to the function.
#' }
#' 
#' @importFrom graphics par
#' @export
"plot.MDSVfilter" <- function(x, ...) {
  stopifnot(class(x) == "MDSVfilter")
  
  if(x$ModelType == "Univariate log-return")                   ModelType <- 0
  if(x$ModelType == "Univariate realized variances")           ModelType <- 1
  if(x$ModelType == "Joint log-return and realized variances") ModelType <- 2
  
  para      <- x$estimates
  sig       <- sqrt(volatilityVector(para,K=x$K,N=x$N))
  LEVIER    <- x$LEVIER
  
  proba_lis <- x$smoothed_proba
  data      <- as.matrix(x$data)
  n<-nrow(data)
  s<-numeric(n)
  for(i in 1:n) s[i]<-which.max(proba_lis[,i])
  
  if(LEVIER){
    Levier<-levierVolatility(ech=data[,1],para=para,Model_type=ModelType)$`Levier`
    V_t<-sig[s]*Levier
  }else{
    V_t<-sig[s]
  }
  
  x              <- list(V_t       = V_t,
                         data      = x$data,
                         dates     = x$dates,
                         ModelType = ModelType)
  
  x              <- c(x, list(... = ...))
  
  class(x)       <- "plot.MDSVfilter"
  x
}

#' @rdname plot.MDSVfilter
#' @export
"print.plot.MDSVfilter" <- function(x, ...) {
  
  V_t         <- x$V_t
  ModelType   <- x$ModelType
  data        <- as.matrix(x$data)
  dates       <- x$dates
  
  .pardefault <- do.call("par", list(mar=c(6,6,4,4)))
  
  if(ModelType==0){
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2.3,2.3), respect = FALSE)
    tmp           <- c(list(x = dates,
                            y = V_t,
                            type = "l",
                            main = "Filtred Volatilities",
                            ylab = "Filtred Volatilities",
                            xaxt='n',
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(0, 4.1, 4.1, 2.1)))
    do.call("plot", tmp)
    
    tmp           <- c(list(x = dates,
                            y = data[,1],
                            type = "l",
                            ylab = "Log-returns",
                            xlab="Date",
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(4.1, 4.1, 0, 2.1)))
    do.call("plot", tmp)
    
  }else if(ModelType==1){
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2.3,2.3), respect = FALSE)
    tmp           <- c(list(x = dates,
                            y = V_t,
                            type = "l",
                            main = "Filtred Volatilities",
                            ylab = "Filtred Volatilities",
                            xaxt='n',
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(0, 4.1, 4.1, 2.1)))
    do.call("plot", tmp)
    
    tmp           <- c(list(x = dates,
                            y = sqrt(data[,2]),
                            type = "l",
                            ylab = "Realized Volatilities",
                            xlab="Date",
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(4.1, 4.1, 0, 2.1)))
    do.call("plot", tmp)
  }else{
    layout(matrix(1:3, ncol = 1), widths = 1, heights = c(2.3,2,2.3), respect = FALSE)
    tmp           <- c(list(x = dates,
                            y = V_t,
                            type = "l",
                            main = "Filtred Volatilities",
                            ylab = "Filtred Volatilities",
                            xaxt='n',
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(0, 4.1, 4.1, 2.1)))
    do.call("plot", tmp)
    
    tmp           <- c(list(x = dates,
                            y = sqrt(data[,2]),
                            type = "l",
                            ylab = "Realized Volatilities",
                            xaxt='n',
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(0, 4.1, 0, 2.1)))
    do.call("plot", tmp)
    
    tmp           <- c(list(x = dates,
                            y = data[,1],
                            type = "l",
                            ylab = "Log-returns",
                            xlab="Date",
                            ...  = ...), x[-(1:4)])
    do.call("par", list(mar=c(4.1, 4.1, 0, 2.1)))
    do.call("plot", tmp)
    
  }
  
  oldw <- getOption("warn")
  par(.pardefault)
  options(warn = oldw)
  
  invisible(x)
}

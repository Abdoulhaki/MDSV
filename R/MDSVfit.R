#' @title MDSV Fitting
#' @description Method for fitting the MDSV model on log-retruns and realized variances (uniquely or jointly).
#' @param N An integer designing the number of components for the MDSV process
#' @param K An integer designing the number of states of each MDSV process component
#' @param data A univariate or bivariate data matrix. Can only be a matrix of 1 or 2 columns. If data has 2 columns, the first one has to be the log-returns and the second the realized variances.
#' @param ModelType An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.
#' @param LEVIER if `TRUE`, estime the MDSV model with leverage.
#' @param ... further options for the \code{\link{solnp}} solver of the \pkg{Rsolnp} package.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item{ModelType}{type of model to be fitted.}
#'     \item{LEVIER}{wheter the fit take the leverage effect into account or not.}
#'     \item{N}{number of components for the MDSV process.}
#'     \item{K}{number of states of each MDSV process component.}
#'     \item{estimates}{estimated parameters.}
#'     \item{LogLikelihood}{log-likelihood of the model on the data.}
#'     \item{AIC}{Akaike Information Criteria of the model on the data.}
#'     \item{BIC}{Bayesian Information Criteria of the model on the data.}
#'     \item{data}{data use for the fitting.}
#' }
#' 
#' @details 
#' The MDSV optimization routine set of feasible starting points which are used to initiate the MDSV recursion. The 
#' likelihood calculation is performed in \code{C++} through the \pkg{Rcpp} package. The optimization is perform using
#' the solnp solver of the \pkg{Rsolnp} package and additional options can be supply to the fonction. 
#' The leverage effect is taken into account according to the FHMV model (see Augustyniak et al., 2019). While fitting an
#' univariate realized variances data, log-returns are required to add leverage effect.
#' AIC and BIC are computed using the formulas : 
#' \itemize{
#' \item{AIC : }{\eqn{L - k}}
#' \item{BIC : }{\eqn{L - (k/2)*log(n)}}
#' }
#' where \eqn{L} is the log-likelihood, \eqn{k} is the number of parameters and \eqn{n} the number of observations in the dataset.
#' The \link[base]{class} of the output of this function is \code{MDSVfit}. This class has a \link[base]{summary}, 
#' \link[base]{print} and \link[base]{plot} \link[utils]{methods} to summarize, print and plot the results. See 
#' \code{\link{summary.MDSVfit}}, \code{\link{print.MDSVfit}} and \code{\link{plot.MDSVfit}} for more details.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' 
#' @seealso For filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @examples 
#' \dontrun{
#' # MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
#' data      <- data(sp500)  # Data loading
#' N         <- 2           # Number of components
#' K         <- 3           # Number of states
#' ModelType <- 0           # Univariate log-returns
#' LEVIER    <- FALSE       # No leverage effect
#' 
#' # Model estimation
#' out       <- MDSVfit(K = K, N = N, data = donne, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out)
#' # Plot
#' plot(out,"both")
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' data      <- data(nasdaq)  # Data loading
#' N         <- 3             # Number of components
#' K         <- 3             # Number of states
#' ModelType <- 2             # Joint log-returns and realized variances
#' LEVIER    <- TRUE          # No leverage effect
#' 
#' # Model estimation
#' out       <- MDSVfit(K = K, N = N, data = donne, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out)
#' # Plot
#' plot(out,"nic")
#' }

#' @export
#' @importFrom Rsolnp solnp 
MDSVfit<-function(N,K,data,ModelType=0,LEVIER=FALSE,...){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVfit(): input N and K must be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVfit(): input N and K must be integer!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVfit(): input ModelType must be numeric!")
  }else if(!(ModelType %in% c(0,1,2))){
    stop("MDSVfit(): input ModelType must be 0, 1 or 2!")
  }
  
  if(!is.logical(LEVIER)) {
    stop("MDSVfit(): input LEVIER must be logical!")
  }
  
  if ( (!is.numeric(data)) || (!is.matrix(data))  ) {
    stop("MDSVfit(): input data must be numeric matrix!")
  }
  
  k <- ncol(data)
  T <- nrow(data)
  
  if((k==1) & (!(ModelType == 0) & (!((ModelType == 1) & (LEVIER == FALSE))))){
    stop("MDSVfit(): improperly data dimensioned matrices!")
  }
  
  if((k==1) & (ModelType == 1) & sum(data<0)>0 ){
    stop("MDSVfit(): data must be positive!")
  } 
  if((k==1) & (ModelType == 1)) data <- matrix(c(rep(1,nrow(data)),data),nrow(data),2)
  
  if(k==2) if(sum(data[,2]<0)>0 ){
    stop("MDSVfit(): data second colomn must be positive!")
  }
  
  ### Some constants
  ctrl <- list(... = ...)
  if(!("TOL" %in% names(ctrl))) ctrl<-c(ctrl, list(TOL=1e-15))
  if(!("trace" %in% names(ctrl))) ctrl<-c(ctrl, list(trace=0))
  
  para<-c(0.52,0.99, 2.77,sqrt(var(data[,1])),0.72)
  if(ModelType==1) para <- c(para,2.10)
  if(ModelType==2) para <- c(para,-0.5,	0.93,	-0.09,	0.04,	2.10)
  if(LEVIER)       para <- c(para,0.78,0.87568)
  para_tilde <- natWork(para=para,LEVIER=LEVIER,Model_type=ModelType)
  t<-0
  
  while(!is.numeric(logLik(ech=donne,para_tilde=para_tilde,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N))){
    if(t>1000) stop("MDSVfit(): Fail to find an optimization starting point")
    
    para_tilde <- para_tilde+mean(para_tilde)*sample(c(-1,1),length(para_tilde),T)
    
    t<-t+1
  }
  
  
  oldw <- getOption("warn")
  options(warn = -1)
  opt<-try(solnp(pars=para_tilde,fun=logLik,ech=data,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
  
  if (class(opt) =='try-error'){
    stop("MDSVfit(): Optimization ERROR")
  }
  options(warn = oldw)
  
  ### Results
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)        vars <- c(vars,"l","theta")
  
  params<-workNat(para=opt$pars,LEVIER=LEVIER,Model_type=ModelType)
  names(params)<-vars
  
  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  if(N==1) params<-params[(vars[!(vars=='b')])]
  
  out<-list(ModelType     = Model_type, 
            LEVIER        = LEVIER, 
            N             = N, 
            K             = K, 
            estimates     = params, 
            LogLikelihood = -as.numeric(opt$values[length(opt$values)]),
            AIC           = -as.numeric(opt$values[length(opt$values)])-length(params), 
            BIC           = -as.numeric(opt$values[length(opt$values)])-0.5*length(params)*log(T),
            data          = data)
  
  class(out) <- "MDSVfit"
  
  return(out)
}

#' @title Summarize, print and plot MDSV Fitting
#' @description Summary, print and plot methods for the class `MDSVfit` as returned by the function \link{MDSVfit}.
#' @param object An object of class `MDSVfit`, output of the function \code{\link{MDSVfit}}.
#' @param x An object of class `summary.MDSVfit`, output of the function \code{\link{summary.MDSVfit}},
#' class `MDSVfit` of the function \code{\link{MDSVfit}} or `plot.MDSVfit` of the function \code{\link{plot.MDSVfit}}.
#' @param plot.type A character designing the type of plot. `dis` for the stationnary distribution of the volatilities,
#'  `nic` for the New Impact Curve (see. Engle and Ng, 1993).
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item{ModelType}{type of model to be fitted.}
#'     \item{LEVIER}{wheter the fit take the leverage effect into account or not.}
#'     \item{N}{number of components for the MDSV process.}
#'     \item{K}{number of states of each MDSV process component.}
#'     \item{estimates}{estimated parameters.}
#'     \item{LogLikelihood}{log-likelihood of the model on the data.}
#'     \item{AIC}{Akaike Information Criteria of the model on the data.}
#'     \item{BIC}{Bayesian Information Criteria of the model on the data.}
#'     \item{data}{data use for the fitting.}
#'     \item{...}{further arguments passed to the function.}
#' }
#'
#' @details 
#' `dis` as argument `plot.type` lead to plot the stationnary distribution of the Markov chain process MDSV. The leverage
#' effect is not took into account for that plot.
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @references 
#' Engle, R. F., & Ng, V. K. (1993). Measuring and testing the impact of news on volatility. 
#' \emph{The journal of finance}, 48(5), 1749-1778. \url{ https://doi.org/10.1111/j.1540-6261.1993.tb05127.x}
#' 
#' @export
"summary.MDSVfit" <- function(object, ...){
  stopifnot(class(object) == "MDSVfit")
  
  out<-c(object,...)
  
  class(out) <- "summary.MDSVfit"
  return(out)
}


#' @rdname summary.MDSVfit
#' @export
"print.summary.MDSVfit" <- function(x, ...){
  stopifnot(class(x) == "summary.MDSVfit")
  
  cat("=================================================\n")
  cat(paste0("=================  MDSV fitting =================\n"))
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
  cat(paste0("LogLikelihood : ", round(x$LogLikelihood,2),"\n\n"))
  
  cat("Information Criteria \n")
  cat("------------------------------------------------- \n")
  cat(paste0("AIC \t: ", round(x$AIC,2),"\n"))
  cat(paste0("BIC \t: ", round(x$BIC,2),"\n"))
  
  invisible(x)
}


#' @rdname summary.MDSVfit
#' @export
"print.MDSVfit" <- function(x, ...) {
  stopifnot(class(x) == "MDSVfit")
  print(summary(x, ...))
}


#' @rdname summary.MDSVfit
#' @importFrom graphics par
#' @export
"plot.MDSVfit" <- function(x, plot.type = c("dis", "nic"), ...) {
  stopifnot(class(x) == "MDSVfit")
  stopifnot(prod(plot.type %in% c("dis", "nic"))==1)
  
  para      <- x$estimates
  if(is.null(para["b"])) {
    para        <- c(para,1)
    names(para) <- c(names(para),"b")
    para        <- para[c("omega","a","b",names(para)[3:length(para)])]
  }
  
  sig       <- sqrt(volatilityVector(para,K=x$K,N=x$K))
  prob      <- probapi(para["omega"],K=x$K,N=x$K)
  
  x              <- list(ModelType   = x$ModelType,
                         N           = x$N,
                         K           = x$K,
                         LEVIER      = x$LEVIER,
                         estimates   = para,
                         sig         = sig,
                         prob        = prob,
                         data        = x$data,
                         plot.type   = plot.type)
  
  if (is.null(match.call()$mfrow)) {
    nrow         <- 1
    ncol         <- 1
    if(length(plot.type)==2) ncol <- 2
    
    x            <- c(x, list(mfrow = c(nrow, ncol)))
  } 
  
  x              <- c(x, list(... = ...))
  
  class(x)       <- "plot.MDSVfit"
  x
}


#' @rdname summary.MDSVfit
#' @export
"print.plot.MDSVfit" <- function(x, ...) {
  
  if(x$ModelType=="Univariate log-return")                   ModelType <- 0
  if(x$ModelType=="Univariate realized variances")           ModelType <- 1
  if(x$ModelType=="Joint log-return and realized variances") ModelType <- 2
  
  N           <- x$N
  K           <- x$K
  LEVIER      <- x$LEVIER
  para        <- x$estimates
  sig         <- x$sig
  prob        <- x$prob
  data        <- as.matrix(x$data)
  plot.type   <- x$plot.type
  
  do.call("par", c(x[-(1:9)], list(...)))
  
  if("dis" %in% plot.type){
    temp<-aggregate(x=prob,by=list(round(sig,4)),FUN="sum")
    
    tmp           <- c(list(x = temp[,1],
                            y = temp[,2],
                            type = "l",
                            main = "Density plot : Stationnary distribution of the volatilities",
                            xlab = "volatilities",
                            ylab = "probabilities",
                            ...  = ...), x[-(1:10)])
    do.call("plot", tmp)
  }
  if("nic" %in% plot.type){
    temp      <- logLik2(ech=data, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N)
    proba_lis <- temp$smoothed_proba
    
    n<-nrow(data)
    s<-numeric(n)
    for(i in 1:n) s[i]<-which.max(proba_lis[,i])
    
    if(LEVIER){
      Levier<-levierVolatility(ech=data[,1],para=para,Model_type=ModelType)$`Levier`
      V_t<-sig[s]*Levier
    }else{
      V_t<-sig[s]
    }
    
    ind<-order(data[1:(n-1),1])
    lo <- loess(V_t[ind]~data[ind,1])
    
    tmp           <- c(list(x=data[ind,1],
                            y=predict(lo),
                            type = "l",
                            main = "New Impact Curve",
                            xlab = "log-returns (rtm1)",
                            ylab = "volatilities (Vt)",
                            ...  = ...), x[-(1:10)])
    do.call("plot", tmp)
    
  }
  
  par(mfrow = c(1, 1))
  
  invisible(x)
}

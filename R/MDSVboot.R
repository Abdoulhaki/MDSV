#' @title MDSV forecasting via Bootstrap
#' @description Method for forecasting volatility using the MDSV model on log-retruns and realized variances (uniquely or jointly).
#' @param fit An object of \link[base]{class} \code{MDSVfilter} obtain after fitting the model using \code{\link{MDSVfilter}}.
#' @param n.ahead An integer designing the forecast horizon.
#' @param n.bootpred An integer designing the number of simulation based re-fits the model. Not relevant for one horizon forecast or for non-leverage type model.
#' @param rseed An integer use to initialize the random number generator for the resampling with replacement method (if not supplied take randomly).
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item estimates : estimated parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item data : data use for the fitting.
#'     \item dates : vector or names of data designing the dates.
#'     \item n.ahead : integer designing the forecast horizon.
#'     \item n.bootpred : integer designing the number of simulation based re-fits used to generate the parameter distribution.
#'     \item rt_sim : matrix of log-returns forecast simulated where the row stand for the simulations and the columns for the horizon.
#'     \item rt2 : vector of mean by column of the square of rt_sim.
#'     \item rvt_sim : matrix of realized variances forecast simulated where the row stand for the simulations and the columns for the horizon.
#'     \item rvt : vector of mean by column of rvt_sim.
#' }
#' 
#' @details 
#' The MDSVboot perform the forecasting of the model a different horizon. The forecasting is based on a close form where the estimation does not
#' involve leverage effect (see Hamilton, 1989. chapter 22 for hidden markov model forecasting). But to take into account the leverage effect,
#' the forecasting is perform by a bootstrap analysis. The innovations are bootstrapped using a standard normal distribution. This process 
#' is performed in \code{C++} through the \pkg{Rcpp} package. The leverage effect is taken into account according to the FHMV model 
#' (see Augustyniak et al., 2019). For the univariate realized variances forecasting, log-returns are required to add leverage effect. 
#'
#' The \link[base]{class} of the output of this function is \code{\link{MDSVboot}}. This class has a \link[base]{summary} and 
#' \link[base]{print} \link[utils]{methods} to summarize and print the results. See 
#' \code{\link{summary.MDSVboot}} and \code{\link{print.MDSVboot}} for more details.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, simulation \code{\link{MDSVsim}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @examples 
#' \dontrun{
#' # MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
#' data(sp500)         # Data loading
#' N         <- 2      # Number of components
#' K         <- 3      # Number of states
#' ModelType <- 0      # Univariate log-returns
#' LEVIER    <- FALSE  # No leverage effect
#' 
#' # Model estimation
#' out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = ModelType, LEVIER = LEVIER)
#' # Model forecasting (no bootstrapp is need as no leverage)
#' para      <-out_fit$estimates # parameter
#' out       <- MDSVboot(fit = out_fit, n.ahead = 100, rseed = 125)
#' # Summary
#' summary(out)
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' data(nasdaq)       # Data loading
#' N         <- 3     # Number of components
#' K         <- 3     # Number of states
#' ModelType <- 2     # Joint log-returns and realized variances
#' LEVIER    <- TRUE  # No leverage effect
#' 
# Model estimation
#' out_fit   <- MDSVfit(K = K, N = N, data = nasdaq, ModelType = ModelType, LEVIER = LEVIER)
#' # Model bootstrap forecasting
#' out       <- MDSVboot(fit = out_fit, n.ahead = 100, n.bootpred = 10000, rseed = 349)
#' # Summary
#' summary(out)
#' 
#' }
#' 
#' @import Rcpp
#' @export
#' @importFrom mhsmm sim.mc
MDSVboot<-function(fit,n.ahead=100,n.bootpred=500,rseed=NA){
  stopifnot(class(fit) %in% c("MDSVfit","MDSVfilter"))
  
  if ( (!is.numeric(n.ahead)) || (!is.numeric(n.bootpred)) ) {
    stop("MDSVboot(): inputs must be numeric!")
  }else if(!(n.ahead%%1==0) || !(n.bootpred%%1==0)){
    stop("MDSVboot(): input n.ahead and n.bootpred must be integer!")
  }else if((n.ahead<1) || (n.bootpred<1)){
    stop("MDSVfit(): input n.ahead and n.bootpred must be positive!")
  }
  
  if ((!is.numeric(rseed)) & !is.na(rseed)) {
    print("MDSVboot() WARNING: input rseed must be numeric! rseed set to random")
    rseed <- sample.int(.Machine$integer.max,1)
  }else if(is.numeric(rseed)){ 
    if(!(rseed%%1==0)){
      rseed <- floor(rseed)
      print(paste0("MDSVboot() WARNING: input rseed must be an integer! rseed set to ",rseed))
    }
    set.seed(rseed)
  }
  
  if(fit$ModelType == "Univariate log-return")                   ModelType <- 0
  if(fit$ModelType == "Univariate realized variances")           ModelType <- 1
  if(fit$ModelType == "Joint log-return and realized variances") ModelType <- 2
  para      <- fit$estimates
  N         <- fit$N
  K         <- fit$K
  LEVIER    <- fit$LEVIER
  data      <- as.matrix(fit$data)
  
  k <- ncol(data)
  T <- nrow(data)
  
  if(!is.null(names(data))) {
    dates <- as.Date(names(data)[1:T])
  }else {
    dates<- 1:T
  }
  
  out<-c(fit, list(dates      = dates,
                   n.ahead    = n.ahead,
                   n.bootpred = n.bootpred))
  
  l<-logLik2(ech=data, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N, t=T)
  
  pi_0 <- l$w_hat
  sig  <- volatilityVector(para=para,K=K,N=N)
  matP <- P(para=para,K=K,N=N)
  
  MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(n.ahead,n.bootpred)),n.ahead,n.bootpred,byrow=FALSE)) #simulation of Markov chain
  
  z_t<-matrix(rnorm(n.bootpred*n.ahead),nrow=n.bootpred,ncol=n.ahead)
  
  if(LEVIER){
    Levier     <- rep(1,n.bootpred)%*%t(levierVolatility(ech=data[((T-200):T),1],para=para,Model_type=ModelType)$`Levier`)
    sim        <- R_hat( H          = n.ahead,
                         ech        = data[((T-200):T),1],
                         MC_sim     = MC_sim,
                         z_t        = z_t,
                         Levier     = Levier,
                         sig        = sig,
                         para       = para,
                         Model_type = ModelType,
                         N          = N)
    if(!(ModelType==1)){
      rt2        <- sim$`rt2`
      rt_sim     <- sim$`rt_sim`
      rt2        <- rt2[(ncol(rt2)-n.ahead+1):ncol(rt2)]
      rt_sim     <- rt_sim[,(ncol(rt_sim)-n.ahead+1):ncol(rt_sim)]
      out        <- c(out,list(rt2 = rt2, rt_sim = rt_sim))
      if(ModelType==2){
        rvt      <- sim$`rvt`
        rvt_sim  <- sim$`rvt_sim`
        rvt      <- rvt[(ncol(rvt)-n.ahead+1):ncol(rvt)]
        rvt_sim  <- rvt_sim[,(ncol(rvt_sim)-n.ahead+1):ncol(rvt_sim)]
        out      <- c(out,list(rvt = rvt, rvt_sim = rvt_sim))
      }
    }else{
      sim1       <- sim$LevierMat
      sim        <- sim$Levier
      rvt        <- (f_sim(n.ahead,sig,pi_0,matP)$`rt2`)*(sim[(length(sim)-n.ahead+1):length(sim)])
      rvt_sim    <- (f_sim(n.ahead,sig,pi_0,matP)$`rt2`)*(sim1[,(ncol(sim1)-n.ahead+1):ncol(sim1)])
      out        <- c(out,list(rvt = rvt, rvt_sim = rvt_sim))
    }
  }else{
    if(!(ModelType==2)){
      sim     <- f_sim(n.ahead,sig,pi_0,matP)
      rvt     <- rt2 <- sim$`rt2`
      if(ModelType==0) {
        out   <- c(out,list(rt2 = rt2))
      }else{
        out   <- c(out,list(rvt = rvt))
      }
    }else{
      xi      <- para[6];
      varphi  <- para[7];
      delta1  <- para[8];
      delta2  <- para[9];
      shape   <- para[10];
      sim     <- f_sim(n.ahead,sig,pi_0,matP,varphi,xi,shape,delta1,delta2)
      rt2     <- sim$`rt2`
      rvt     <- sim$`rvt`
      out     <- c(out,list(rt2 = rt2, rvt = rvt))
    }
  }
  
  class(out) <- "MDSVboot"
  
  return(out)
}

#' @title Summarize and print MDSV bootstrap forecasting
#' @description Summary and print methods for the class \link{MDSVboot} as returned by the function \code{\link{MDSVboot}}.
#' @param object An object of class \link{MDSVboot}, output of the function \code{\link{MDSVboot}}.
#' @param x An object of class \link{summary.MDSVboot}, output of the function \code{\link{summary.MDSVboot}}
#' or class \link{MDSVboot} of the function \code{\link{MDSVboot}}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item estimates : estimated parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item data : data use for the fitting.
#'     \item dates : vector or names of data designing the dates.
#'     \item n.ahead : integer designing the forecast horizon.
#'     \item n.bootpred : integer designing the number of simulation based re-fits used to generate the parameter distribution.
#'     \item rt_sim : matrix of log-returns forecast simulated where the row stand for the simulations and the columns for the horizon.
#'     \item rt2 : vector of mean by column of the square of rt_sim.
#'     \item rvt_sim : matrix of realized variances forecast simulated where the row stand for the simulations and the columns for the horizon.
#'     \item rvt : vector of mean by column of rvt_sim.
#' }
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @export

"summary.MDSVboot" <- function(object, ...){
  stopifnot(class(object) == "MDSVboot")
  
  out<-c(object,...)
  
  
  class(out) <- "summary.MDSVboot"
  return(out)
}

f.present <- function(X,thr=100){
  tmp<-apply(X[,(ncol(X)-thr+1):ncol(X)],2,quantile)
  Y <- data.frame(min    = tmp[1,],
                  q.25   = tmp[2,],
                  mean   = apply(X[,(ncol(X)-thr+1):ncol(X)],2,mean),
                  median = tmp[3,],
                  q.75   = tmp[4,],
                  max    = tmp[5,])
  return(as.matrix(Y))
}

#' @rdname summary.MDSVboot
#' @export
"print.summary.MDSVboot" <- function(x, ...){
  stopifnot(class(x) == "summary.MDSVboot")
  
  cat("=================================================\n")
  cat(paste0("==========  MDSV Bootstrap Forecasting ==========\n"))
  cat("=================================================\n\n")
  # cat("Conditional Variance Dynamique \n")
  # cat("------------------------------------------------- \n")
  cat(paste0("Model       : MDSV(",x$N,",",x$K,")\n"))
  cat(paste0("Data        : ", x$ModelType,"\n"))
  cat(paste0("Leverage    : ", x$LEVIER,"\n"))
  cat(paste0("n.ahead     : ", x$n.ahead,"\n"))
  cat(paste0("Date (T[0]) : ", x$dates[length(x$dates)],"\n\n"))
  
  n.ahead <- x$n.ahead
  
  if(x$ModelType == "Univariate log-return")                   ModelType <- 0
  if(x$ModelType == "Univariate realized variances")           ModelType <- 1
  if(x$ModelType == "Joint log-return and realized variances") ModelType <- 2
  
  if(x$LEVIER){
    if(!(ModelType == 1)){
      rt_sim <- f.present(X = x$rt_sim, thr = n.ahead)
      rownames(rt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      cat("Log-returns (summary) : \n")
      print(head(round(rt_sim,6), min(n.ahead, 10)))
      if(n.ahead>10)  cat("......................... \n")
      
      if(ModelType == 0){
        rt2_sim <- f.present(X = (x$rt_sim)^2, thr = n.ahead)
        rownames(rt2_sim) <- paste("t", 1:n.ahead, sep="+")
        
        cat(paste0("\n","Sigma (summary) : \n"))
        print(head(round(sqrt(rt2_sim),6), min(n.ahead, 10)))
        if(n.ahead>10)  cat("......................... \n")
      }
    }
    
    if(!(ModelType == 0)){
      rvt_sim <- f.present(X = x$rvt_sim, thr = n.ahead)
      rownames(rvt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      if(ModelType == 2) cat("\n")
      cat("Realized Variances (summary) : \n")
      print(head(round(rvt_sim,6), min(n.ahead, 10)))
      if(n.ahead>0)  cat("......................... \n")
    }
  }else{
    if(!(ModelType == 1)){
      rt_sim <- x$rt2
      names(rt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      cat("Sigma (summary) : \n")
      print(head(round(sqrt(rt_sim),6), min(n.ahead, 10)))
      if(n.ahead>10)  cat("......................... \n")
    }
    
    if(!(ModelType == 0)){
      rvt_sim <- x$rvt
      names(rvt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      if(ModelType == 2) cat("\n")
      cat("Realized Variances (summary) : \n")
      print(head(round(rvt_sim,6), min(n.ahead, 10)))
      if(n.ahead>0)  cat("......................... \n")
    }
  }
  
  invisible(x)
}

#' @rdname summary.MDSVboot
#' @export
"print.MDSVboot" <- function(x, ...) {
  stopifnot(class(x) == "MDSVboot")
  print(summary(x, ...))
}

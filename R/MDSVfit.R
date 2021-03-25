#' @title MDSV Fitting
#' @description Method for fitting the MDSV model on log-retruns and realized variances (uniquely or jointly).
#' @param N An integer designing the number of components for the MDSV process
#' @param K An integer designing the number of states of each MDSV process component
#' @param data A univariate or bivariate data matrix. Can only be a matrix of 1 or 2 columns. If data has 2 columns, the first one has to be the log-returns and the second the realized variances.
#' @param ModelType An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.
#' @param LEVIER if \code{TRUE}, estime the MDSV model with leverage.
#' @param start.pars List of staring parameters for the optimization routine. These are not usually required unless the optimization has problems converging.
#' @param ... Further options for the \code{\link{solnp}} solver of the \pkg{Rsolnp} package.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item convergence : 0 if convergence, 1 if not.
#'     \item estimates : estimated parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item data : data use for the fitting.
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
#' The \link[base]{class} of the output of this function is \code{\link{MDSVfit}}. This class has a \link[base]{summary}, 
#' \link[base]{print} and \link[base]{plot} \link[utils]{methods} to summarize, print and plot the results. See 
#' \code{\link{summary.MDSVfit}}, \code{\link{print.MDSVfit}} and \code{\link{plot.MDSVfit}} for more details.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' 
#' @seealso For filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}}, simulation \code{\link{MDSVsim}} and rolling estimation and forecast \code{\link{MDSVroll}}.
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
#' out       <- MDSVfit(K = K, N = N, data = sp500, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out)
#' # Plot
#' plot(out,c("dis","nic"))
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' data(nasdaq)      # Data loading
#' N         <- 3    # Number of components
#' K         <- 3    # Number of states
#' ModelType <- 2    # Joint log-returns and realized variances
#' LEVIER    <- TRUE # No leverage effect
#' 
#' # Model estimation
#' out       <- MDSVfit(K = K, N = N, data = nasdaq, ModelType = ModelType, LEVIER = LEVIER)
#' # Summary
#' summary(out)
#' # Plot
#' plot(out,"nic")
#' }

#' @export
#' @import Rcpp
#' @importFrom Rsolnp solnp gosolnp
MDSVfit<-function(N,K,data,ModelType=0,LEVIER=FALSE,start.pars=list(),...){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVfit() ERROR: input N and K must be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVfit() ERROR: input N and K must be integer!")
  }else if(K<2){
    stop("MDSVfit() ERROR: input K must be greater than one!")
  }else if(N<1){
    stop("MDSVfit() ERROR: input N must be positive!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVfit() ERROR: input ModelType must be numeric!")
  }else if(!(ModelType %in% c(0,1,2))){
    stop("MDSVfit() ERROR: input ModelType must be 0, 1 or 2!")
  }
  
  if(!is.logical(LEVIER)) {
    stop("MDSVfit(): input LEVIER must be logical!")
  }
  
  if ( (!is.numeric(data)) || (!is.matrix(data))  ) {
    stop("MDSVfit() ERROR: input data must be numeric matrix!")
  }
  
  k <- ncol(data)
  T <- nrow(data)
  
  if((k==1) & (!(ModelType == 0) & (!((ModelType == 1) & (LEVIER == FALSE))))){
    stop("MDSVfit() ERROR: improperly data dimensioned matrices!")
  }
  
  if((k==1) & (ModelType == 1) & sum(data<0)>0 ){
    stop("MDSVfit() ERROR: data must be positive!")
  } 
  if((k==1) & (ModelType == 1)) data <- matrix(c(rep(1,nrow(data)),data),nrow(data),2)
  
  if(k==2) if(sum(data[,2]<0)>0 ){
    stop("MDSVfit() ERROR: data second colomn must be positive!")
  }
  
  ### Some constants
  ctrls <- list(... = ...)
  ctrl  <- NULL
  if(!("control" %in% names(ctrls))) ctrl<-ctrls$control
  if(!("TOL" %in% names(ctrls))){
    ctrl<-c(ctrl, list(TOL=1e-15))
  }else{
    if(!is.numeric(ctrls$TOL)){
      ctrl<-c(ctrl, list(TOL=1e-15))
    }else{
      ctrl<-c(ctrl, list(TOL=ctrls$TOL))
    }
  }
  if(!("trace" %in% names(ctrls))){
    ctrl<-c(ctrl, list(trace=0))
  }else{
    if(!is.numeric(ctrls$trace)){
      ctrl<-c(ctrl, list(trace=0))
    }else{
      ctrl<-c(ctrl, list(trace=ctrls$trace))
    }
  }
  
  vars<-c("omega","a","b","sigma","v0")
  LB<-rep(-10,5)
  UB<-rep(10,5)
  if(ModelType==1) {
    vars <- c(vars,"shape")
    LB<-c(LB,-10)
    UB<-c(UB,10)
  }else if(ModelType==2) {
    vars <- c(vars,"xi","varphi","delta1","delta2","shape")
    LB<-c(LB,rep(-2,4),-10)
    UB<-c(UB,rep(2,4),10)
  }
  if(LEVIER) {
    vars <- c(vars,"l","theta")
    LB<-c(LB,-10,-10)
    UB<-c(UB,10,10)
  }
  
  if("LB" %in% names(ctrls)){
    if((length(ctrls$LB)==length(vars)) & is.numeric(ctrls$LB)){
      LB<-ctrls$LB
    }else{
      print("MDSVfit() WARNING: Incorrect Lower Bound! set to default.")
    }
  }
  
  if("UB" %in% names(ctrls)){
    if((length(ctrls$UB)==length(vars)) & is.numeric(ctrls$UB)){
      UB<-ctrls$UB
    }else{
      print("MDSVfit() WARNING: Incorrect Upper Bound! set to default.")
    }
  }
  
  if("n.restarts" %in% names(ctrls)){
    if(is.numeric(ctrls$n.restarts)){
      n.restarts<-ctrls$n.restarts
      if(!(n.restarts%%1==0)){
        print("MDSVfit() WARNING: Incorrect n.restarts! set to 1.")
        n.restarts<-1
      }
    }else{
      print("MDSVfit() WARNING: Incorrect n.restarts! set to 1.")
      n.restarts<-1
    }
  }else{
    n.restarts<-1
  }
  
  if("n.sim" %in% names(ctrls)){
    if(is.numeric(ctrls$n.sim)){
      n.sim<-ctrls$n.sim
      if(!(n.sim%%1==0)){
        print("MDSVfit() WARNING: Incorrect n.sim! set to 200.")
        n.sim<-200
      }
    }else{
      print("MDSVfit() WARNING: Incorrect n.sim! set to 200.")
      n.sim<-200
    }
  }else{
    n.sim<-200
  }
  
  if("cluster" %in% names(ctrls)){
    cluster <- ctrls$cluster
    if ( (!is.null(cluster)) & (!("cluster" %in% class(cluster))) ) {
      print("MDSVroll() WARNING: input cluster must be a cluster object the package parallel or set to NULL! cluster set to NULL")
      cluster<-NULL
    }
  }else{
    cluster<-NULL
  }
  
  tmp <- c(0.52,0.85, 2.77,sqrt(var(data[,1])),0.72)
  if(ModelType==1) tmp <- c(tmp,2.10)
  if(ModelType==2) tmp <- c(tmp,-1.5,	0.72,	-0.09,	0.04,	2.10)
  if(LEVIER)       tmp <- c(tmp,1.5,0.87568)
  names(tmp)           <- vars
  
  if(length(start.pars)){
    if(!is.null(names(start.pars))){
      para <- start.pars[vars]
      if(!(sum(is.na(para))==0)){
        print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
        nam_       <- names(para)[is.na(para)]
        para[nam_] <- tmp[nam_]
      }
    }else{
      if(length(start.pars)==length(vars)){
        para        <- start.pars
        names(para) <- vars
        
        if((para["omega"]>1) || (para["omega"]<0) || (para["a"]>1) || (para["a"]<0) || (para["b"]<=1) ||
           (para["sigma"]<=0) || (para["v0"]>1) || (para[v0]<0)) {
          print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
          para <- NULL
        }else if(((ModelType==1) & (para["shape"]<=0)) || ((ModelType==2) & (para[10]<=0))){
          print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
          para <- NULL
        }else if(LEVIER){
          if(ModelType==0){
            if( (para[6]<0) || (para[7]>1) || (para[7]<0) ){
              print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
              para <- NULL
            }
          }else if(ModelType==1){
            if((para[7]<=0) || (para[8]>1) || (para[8]<0)){
              print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
              para <- NULL
            }
          }else if(ModelType==2){
            if((para[11]<=0) || (para[12]>1) || (para[12]<0)){
              print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
              para <- NULL
            }
          }
        }
        
      }else{
        print("MDSVfit() WARNING: Incorrect start.pars! set to default.")
        para <- NULL
      }
    }
  }else{
    para <- NULL
  }
  
  oldw        <- getOption("warn")
  options(warn = -1)
  if(!is.null(para)){
    para_tilde <- natWork(para=para,LEVIER=LEVIER,Model_type=ModelType)
    opt<-try(solnp(pars=para_tilde,fun=logLik,ech=data,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
  }else{
    opt<-try(gosolnp(pars=NULL,fun=logLik,ech=data,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl,
                     LB=LB,UB=UB,n.restarts=n.restarts,n.sim=n.sim,cluster=cluster),silent=T)
  }
  options(warn = oldw)
  
  if(class(opt) =='try-error'){
    stop("MDSVfit() ERROR: Fail to converge!")
  }else{
    params<-workNat(opt$pars, LEVIER=LEVIER, Model_type=ModelType)
    names(params)<-vars
  } 
  
  if((round(params["omega"],5)==0) || (round(params["omega"],5)==1) || (round(params["a"],5)==0) ||
           (round(params["a"],5)==1) ||(round(params["v0"],5)==0) ||(round(params["v0"],5)==1) ||
           (round(params["b"],5)==1) || ( 0 %in% round(params["a"]^(params["b"]^c(0:N)),5) ) ||
           ( LEVIER & ((round(params["theta"],5)==0) || (round(params["theta"],5)==1))) ||
           ( LEVIER & (round(params["l"],5)==0)) ){
    print("MDSVfit() WARNING: Fail to converge! Return the best results found.")
  }
  
  ### Results
  
  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  if(N==1) params<-params[(vars[!(vars=='b')])]
  
  out<-list(ModelType     = Model_type, 
            LEVIER        = LEVIER, 
            N             = N, 
            K             = K, 
            convergence   = opt$convergence,
            estimates     = params, 
            LogLikelihood = -as.numeric(opt$values[length(opt$values)]),
            AIC           = -as.numeric(opt$values[length(opt$values)])-length(params), 
            BIC           = -as.numeric(opt$values[length(opt$values)])-0.5*length(params)*log(T),
            data          = data)
  
  class(out) <- "MDSVfit"
  
  return(out)
}


#' @title Summarize, print and plot MDSV Fitting
#' @description Summary, print and plot methods for the class \link{MDSVfit} as returned by the function \code{\link{MDSVfit}}.
#' @param object An object of class \link{MDSVfit}, output of the function \code{\link{MDSVfit}}.
#' @param x An object of class \link{summary.MDSVfit}, output of the function \code{\link{summary.MDSVfit}},
#' class \link{MDSVfit} of the function \code{\link{MDSVfit}} or class \link{plot.MDSVfit} of the function \code{\link{plot.MDSVfit}}.
#' @param plot.type A character designing the type of plot. \code{dis} for the stationnary distribution of the volatilities,
#'  \code{nic} for the New Impact Curve (see. Engle and Ng, 1993).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item convergence : 0 if convergence, 1 if not.
#'     \item estimates : estimated parameters.
#'     \item LogLikelihood : log-likelihood of the model on the data.
#'     \item AIC : Akaike Information Criteria of the model on the data.
#'     \item BIC : Bayesian Information Criteria of the model on the data.
#'     \item data : data use for the fitting.
#'     \item ... : further arguments passed to the function.
#' }
#'
#' @details 
#' \code{dis} as argument \code{plot.type} lead to plot the stationnary distribution of the Markov chain process MDSV. The leverage
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
  if(1-x$convergence){
    convergence <- "Convergence."
  }else{
    convergence <- "No Convergence. Retrun the best result."
  }
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
  cat(paste0("Convergence : ",convergence,"\n"))
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
  if(x$ModelType == "Univariate realized variances") if("nic" %in% plot.type){
    print("print.MDSVfit() WARNING: New Impact Curve is not plot for realized variances")
    plot.type <- plot.type[plot.type != "nic"]
  }
  
  para      <- x$estimates
  if(is.null(para["b"])) {
    para        <- c(para,1)
    names(para) <- c(names(para),"b")
    para        <- para[c("omega","a","b",names(para)[3:length(para)])]
  }
  
  sig       <- sqrt(volatilityVector(para,K=x$K,N=x$N))
  prob      <- probapi(para["omega"],K=x$K,N=x$N)
  
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
                            main = "Density plot : Stationnary \n distribution of the volatilities",
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

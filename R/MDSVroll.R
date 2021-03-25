#' @title MDSV Rolling estimates, volatility forecast and backtesting
#' @description Method for creating rolling estimates and volatility forecast from MDSV models with option for refitting every n periods 
#' with parallel functionality. The rolling estimate can be done using univariate log-returns or realized variances or using joint log-returns and realized variances.
#' @param N An integer designing the number of components for the MDSV process
#' @param K An integer designing the number of states of each MDSV process component
#' @param data A univariate or bivariate data matrix. Can only be a matrix of 1 or 2 columns. If data has 2 columns, the first one has to be the log-returns and the second the realized variances.
#' @param ModelType An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.
#' @param LEVIER if \code{TRUE}, estime the MDSV model with leverage.
#' @param n.ahead An integer designing the forecast horizon.
#' @param n.bootpred An integer designing the number of simulation based re-fits the model. Not relevant for one horizon forecast or for non-leverage type model.
#' @param forecast.length An integer designing the length of the total forecast for which out of sample data from the dataset will be used for testing.
#' @param refit.every Determines every how many periods the model is re-estimated.
#' @param refit.window Whether the refit is done on an expanding window including all the previous data or a moving window where all previous 
#' data is used for the first estimation and then moved by a length equal to refit.every (unless the window.size option is used instead).
#' @param window.size If not NULL, determines the size of the moving window in the rolling estimation, which also determines the first point used.
#' @param calculate.VaR Whether to calculate forecast Value at Risk during the estimation.
#' @param VaR.alpha The Value at Risk tail level to calculate.
#' @param cluster A cluster object created by calling makeCluster from the parallel package. If it is not NULL, then this will be used for parallel estimation of the refits (remember to stop the cluster on completion).
#' @param rseed An integer use to initialize the random number generator for the resampling with replacement method (if not supplied take randomly).
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item ModelType : type of models fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item n.ahead : integer designing the forecast horizon.
#'     \item forecast.length : length of the total forecast for which out of sample data from the dataset will be used for testing.
#'     \item refit.every : Determines every how many periods the model is re-estimated.
#'     \item refit.window : Whether the refit is done on an expanding window including all the previous data or a moving window where all previous 
#' data is used for the first estimation and then moved by a length equal to refit.every (unless the window.size option is used instead).
#'     \item window.size : If not NULL, determines the size of the moving window in the rolling estimation, which also determines the first point used.
#'     \item calculate.VaR : Whether to calculate forecast Value at Risk during the estimation.
#'     \item VaR.alpha : The Value at Risk tail level to calculate.
#'     \item cluster : A cluster object created by calling makeCluster from the parallel package.
#'     \item data : data use for the fitting.
#'     \item dates : vector or names of data designing the dates.
#'     \item estimates : matrix of all the parameters estimates at each date.
#'     \item prevision : matrix of all prevision made a each date.
#' }
#' 
#' @details 
#' This is a wrapper function for creating rolling estimates and volatility forecasts using MDSV models, and optionally calculating the Value at Risk 
#' at specified levels. The argument refit.every determines every how many periods the model is re-estimated. Given a dataset of length n, 
#' it is possible to set how many periods from the end to use for out of sample forecasting (using the forecast.length option). 
#' For rolling 1-ahead forecasts and forecasts without leverage effect, no bootstrap in done and then n.bootpred is not required.
#' However, the Value-at-Risk baskesting is done for 1-ahead forecasts. A very important part of the function is performed 
#' in \code{C++} through the \pkg{Rcpp} package. The leverage effect is taken into account according to the FHMV model 
#' (see Augustyniak et al., 2019). For the univariate realized variances forecasting, log-returns are required to add leverage effect.
#' When cluster is not NULL, the estimations and forecasting are perform with parallel functionalilty which  is entirely based 
#' on the \pkg{parallel} package, and it is up to the user to pass a cluster object, and then stop it once the routine is completed.
#' The \link[base]{class} of the output of this function is \code{MDSVroll}. This class has a \link[base]{summary}, \link[base]{print} and 
#' \link[base]{plot} \link[utils]{methods} to summarize, print and plot the results. See 
#' \code{\link{summary.MDSVroll}}, \code{\link{print.MDSVroll}} and \code{\link{plot.MDSVroll}} for more details.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' @references 
#' JamesHamilton: D.(1994), time series analysis, 1994.
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, simulation \code{\link{MDSVsim}} and bootstrap forecasting \code{\link{MDSVboot}}.
#' 
#' @examples 
#' \dontrun{
#' # MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
#' data(sp500)         # Data loading
#' 
#' N                <- 2            # Number of components
#' K                <- 3            # Number of states
#' ModelType        <- 0            # Univariate log-returns
#' LEVIER           <- FALSE        # No leverage effect
#' n.ahead          <- 100          # Forecast horizon
#' forecast.length  <- 756          # rolling forecast length
#' refit.every      <- 63           # Period to re-estimate the model
#' refit.window     <- "recursive"  # No leverage effect
#' calculate.VaR    <- TRUE
#' VaR.alpha        <- c(0.01, 0.05, 0.1)
#' cluster          <- parallel::makeCluster(parallel::detectCores()[1]-1)
#' rseed            <- 125
#' 
#' # rolling forecasts
#' out<-MDSVroll(N=N, K=K, data=sp500, ModelType=ModelType, LEVIER=LEVIER, n.ahead = n.ahead, 
#'             forecast.length = forecast.length, refit.every = refit.every, refit.window = refit.window, 
#'             window.size=NULL,calculate.VaR = calculate.VaR, VaR.alpha = VaR.alpha, cluster = cluster, rseed = rseed)
#' parallel::stopCluster(cluster)
#' # Summary
#' summary(out, VaR.test=TRUE, Loss.horizon = c(1,5,10,25,50,75,100), Loss.window = 756)
#' # plot
#' plot(out, plot.type=c("VaR","sigma","dens"))
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' data(nasdaq)       # Data loading
#' 
#' # rolling forecasts
#' out<-MDSVroll(N=3, K=3, data=nasdaq, ModelType=2, LEVIER=TRUE, n.ahead = 10, forecast.length = 100, 
#'             refit.every = 25, refit.window = "recursive", window.size = 1000, 
#'             calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05), cluster = NULL, rseed = NA)
#' # Summary
#' summary(out, VaR.test=TRUE, Loss.horizon = c(1,2,4,9), Loss.window = 50)
#' 
#' }
#' 
#' @export
#' @import Rcpp
#' @import KScorrect
#' @importFrom mhsmm sim.mc
#' @importFrom Rsolnp solnp gosolnp
#' @importFrom foreach foreach %dopar%
#' @import doSNOW
MDSVroll<-function(N, K, data, ModelType=0, LEVIER=FALSE, n.ahead = 1, n.bootpred = 10000, forecast.length = 500, 
                   refit.every = 25, refit.window = "recursive", window.size = NULL, 
                   calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, rseed = NA, ...){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVroll(): input N and K must all be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVfit(): input N and K must all be integer!")
  }else if(K<2){
    stop("MDSVfit(): input K must be greater than one!")
  }else if(N<1){
    stop("MDSVfit(): input N must be positive!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVroll(): input ModelType must be numeric!")
  }else if(!(ModelType %in% c(0,1,2))){
    stop("MDSVfit(): input ModelType must be 0, 1 or 2!")
  }
  
  if ( (!is.numeric(data)) || (!is.matrix(data))  ) {
    stop("MDSVroll(): input data must be numeric matrix!")
  }
  
  k <- ncol(data)
  T <- nrow(data)
  
  if(!is.null(names(data))) {
    dates <- as.Date(names(data)[1:T])
  }else if(!is.null(rownames(data))){
    dates <- as.Date(rownames(data)[1:T])
  }else{
    dates <- 1:T
  }
  
  if((k==1) & (!(ModelType == 0) & (!((ModelType == 1) & (LEVIER == FALSE))))){
    stop("MDSVroll(): improperly data dimensioned matrices!")
  }
  
  if((k==1) & (ModelType == 1) & sum(data<0)>0 ){
    stop("MDSVroll(): data must be positive!")
  } 
  if((k==1) & (ModelType == 1)) data <- matrix(c(rep(1,T),data),T,2)
  
  if(k==2) if(sum(data[,2]<0)>0 ){
    stop("MDSVroll(): data second colomn must be positive!")
  }

  if ( (!is.numeric(n.ahead)) || (!is.numeric(forecast.length)) || 
       (!is.numeric(refit.every)) || (!is.numeric(n.bootpred)) ) {
    stop("MDSVroll(): inputs n.ahead, forecast.length, refit.every and n.bootpred must all be numeric!")
  }else if(!(n.ahead%%1==0) || !(forecast.length%%1==0) || !(refit.every%%1==0) || !(n.bootpred%%1==0)){
    stop("MDSVfit(): input n.ahead, forecast.length, refit.every and n.bootpred must all be integer!")
  }else if((n.ahead<1) || (forecast.length<1) || (refit.every<1) || (n.bootpred<1)){
    stop("MDSVfit(): input n.ahead, forecast.length, refit.every and n.bootpred must all be positive!")
  }
  
  if(forecast.length < refit.every){
    print("MDSVroll() WARNING: input forecast.length must be greater or equal to refit.every!")
    refit.every <- forecast.length
  }
  
  if((!is.logical(LEVIER)) || (!is.logical(calculate.VaR))){
    stop("MDSVroll(): input LEVIER and calculate.VaR must all be logical!")
  }
  
  if((ModelType == 1) & (calculate.VaR)){
    print("MDSVroll() WARNING: VaR are compute only for log-returns! calculate.VaR set to FALSE!")
    calculate.VaR <- FALSE
  }
  
  if ( !(refit.window %in% c("recursive", "moving")) ) {
    stop("MDSVroll(): input refit.window must be recursive or moving!")
  }
  
  if ( (!is.null(window.size)) & (!is.numeric(window.size)) ) {
    stop("MDSVroll(): input window.size must be numeric or set to NULL!")
  }
  if(is.numeric(window.size)) if(!(window.size%%1==0)){
    stop("MDSVfit(): input window.size must be integer!")
  }else if(window.size<1){
    stop("MDSVfit(): input window.size must be positive!")
  }
  
  if ( (!is.null(cluster)) & (!("cluster" %in% class(cluster))) ) {
    print("MDSVroll() WARNING: input cluster must be a cluster object the package parallel or set to NULL! cluster set to NULL")
    cluster<-NULL
  }
  
  if(((is.null(window.size)) & (refit.window == "moving"))){
    print("MDSVroll() WARNING: input window.size must be numeric for moving refit! compute by the algorithm")
    window.size <- (T-forecast.length) - 1
  }
  if(!(refit.window == "moving")){
    window.size <- (T-forecast.length) - 1
  }
  
  if ( (!is.numeric(VaR.alpha)) ) {
    stop("MDSVroll(): input VaR.alpha must be numeric!")
  }else if(!prod(VaR.alpha > 0) & !prod(VaR.alpha < 1)){
    stop("MDSVroll(): input VaR.alpha must be between 0 and 1!")
  }
  
  if ((!is.numeric(rseed)) & !is.na(rseed)) {
    print("MDSVroll() WARNING: input rseed must be numeric! rseed set to random")
    rseed <- sample.int(.Machine$integer.max,1)
  }else if(is.numeric(rseed)){ 
    if(!(rseed%%1==0)){
      rseed <- floor(rseed)
      print(paste0("MDSVroll() WARNING: input rseed must be an integer! rseed set to ",rseed))
    }
    set.seed(rseed)
  }
  
  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  
  model<-expand.grid(date = dates[((T-forecast.length+1):T)],rt=1,rvt=1,model="MDSV",N=N,K=K,Levier=LEVIER,ModelType=Model_type)
  vars<-c('predict_loglik', 'loglik', 'AIC', 'BIC',"omega","a","b","sigma","v0")
  
  if(ModelType == 0) {
    model$rvt <- NULL
    model$rt  <- data[((T-forecast.length+1):T),1]
  }else if(ModelType == 1) {
    vars      <- c(vars,"shape")
    model$rt  <- NULL
    model$rvt <- data[((T-forecast.length+1):T),2]
  }else if(ModelType == 2){
    vars      <- c(vars,'Marg_loglik', 'AICm', 'BICm',"varphi","xi","delta1","delta2","shape")
    model$rt  <- data[((T-forecast.length+1):T),1]
    model$rvt <- data[((T-forecast.length+1):T),2]
  }
  if(LEVIER) vars <- c(vars,"l","theta")
  if((calculate.VaR) & (!(ModelType == 1)))
    vars <- c(vars, paste0("VaR",100*(1-VaR.alpha)), paste0("I",100*(1-VaR.alpha)))
  
  model_prev<-model
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  model_rt2 <-expand.grid(date = dates[((T-forecast.length+1):T)])
  vars      <- paste0("r2_for",1:n.ahead)
  model_add <- matrix(0, nrow=nrow(model_rt2), ncol=length(vars))
  colnames(model_add) <- vars
  model_rt2 <- cbind(model_rt2, model_add)
  
  model_rvt <-expand.grid(date = dates[((T-forecast.length+1):T)])
  vars      <- paste0("RV_for",1:n.ahead)
  model_add <- matrix(0, nrow=nrow(model_rvt), ncol=length(vars))
  colnames(model_add) <- vars
  model_rvt <- cbind(model_rvt, model_add)
  
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
  
  tmp <- c(0.52,0.99, 2.77,sqrt(var(data[,1])),0.72)
  if(ModelType==1) tmp <- c(tmp,2.10)
  if(ModelType==2) tmp <- c(tmp,-1.5,	0.72,	-0.09,	0.04,	2.10)
  if(LEVIER)       tmp <- c(tmp,1.5,0.87568)
  names(tmp)           <- vars
  
  if("start.pars" %in% names(ctrls)){
    if((length(ctrls$start.pars)==length(vars)) & is.numeric(ctrls$start.pars)){
      start.pars <- ctrls$start.pars
    }else{
      start.pars <- NULL
    }
  }else{
    start.pars <- NULL
  }
  
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
  
  if(!is.null(para)){
    para_tilde <- natWork(para=para,LEVIER=LEVIER,Model_type=ModelType)
  }
  
  ### Estimation
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- (T-forecast.length) - window.size
  opt<-NULL
  if(is.null(cluster)){
    cat("Estimation step : \n")
    
    pb <- txtProgressBar(min=0, max = forecast.length-1, style = 3)
    
    for(t in 0:(forecast.length-1)){
      ech    <- data[strt:(T-forecast.length+t),]
      
      oldw <- getOption("warn")
      options(warn = -1)
      if(t %in% update_date){
        
        if(!is.null(para)){
          if(!is.null(opt)) para_tilde<-opt$pars
          opt<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
        }else{
          opt<-try(gosolnp(pars=NULL,fun=logLik,ech=ech,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl,
                           LB=LB,UB=UB,n.restarts=n.restarts,n.sim=n.sim,cluster=cluster),silent=T)
        }
        if (class(opt) =='try-error'){
          stop("MDSVroll(): Optimization ERROR")
        }
        para <- workNat(opt$pars,LEVIER=LEVIER,Model_type=ModelType)
        
        if(refit.window == "moving") strt <- strt + refit.every
      }
      options(warn = oldw)
      
      model[t+1,vars] <- round(para,5)
      if(!(ModelType == 1)){
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
        
        pi_0 <- l$w_hat
        sig  <- volatilityVector(para=para,N=N,K=K)
        if(LEVIER){
          Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
          sig    <- sig*Levier
        }
        
        if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
          va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
          if(class(va) =='try-error') {
            va <- try(qmixnorm(VaR.alpha[iter], rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
            if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
          }else{
            model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
          }
        }
        
      }else{
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
      }
      
      model[t+1,"loglik"]         <- -as.numeric(opt$values[length(opt$values)])
      model[t+1,"predict_loglik"] <- l$Pred_loglik
      model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
      model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
      
      if(ModelType == 2){
        model[t+1,"Marg_loglik"]  <- l$Marg_loglik
        model[t+1,'AICm']         <- model[(t+1),"Marg_loglik"]-length(para)
        model[t+1,'BICm']         <- model[(t+1),"Marg_loglik"]-length(para)*log(nrow(ech))/2
      }
      
      setTxtProgressBar(pb, t)
    }
    close(pb)
    
  }else{
    registerDoSNOW(cluster)
    cat("Estimation step 1: \n")
    
    pb <- txtProgressBar(min= 0, max = length(update_date[!(update_date==forecast.length)]), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
      Y<- foreach(t=update_date[!(update_date==forecast.length)], .export=c("solnp", "gosolnp"), 
                .packages =c("Rcpp","RcppArmadillo","RcppEigen","Rsolnp"), .combine = rbind, .options.snow = opts) %dopar% { 

      if(refit.window == "moving") strt <- (T-forecast.length) - window.size + t*refit.every
      ech    <- data[strt:(T-forecast.length+t),]
      if(!is.null(opt)) para_tilde<-opt$pars
      oldw <- getOption("warn")
      options(warn = -1)

      if(!is.null(para)){
        if(!is.null(opt)) para_tilde<-opt$pars
        opt<-solnp(pars=para_tilde,fun=logLik,ech=ech,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl)
      }else{
        opt<-gosolnp(pars=NULL,fun=logLik,ech=ech,Model_type=ModelType,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl,
                         LB=LB,UB=UB,n.restarts=n.restarts,n.sim=n.sim)
      }
      options(warn = oldw)
      
      para <- workNat(opt$pars,LEVIER=LEVIER,Model_type=ModelType)
       
     model[t+1,vars] <- round(para,5)
     if(!(ModelType == 1)){
       l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
       
       pi_0 <- l$w_hat
       sig  <- volatilityVector(para=para,N=N,K=K)
       if(LEVIER){
         Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
         sig    <- sig*Levier
       }
       
       if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
         va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
         if(class(va) =='try-error') {
           va <- try(qmixnorm(VaR.alpha[iter], rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
           if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
         }else{
           model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
         }
       }
     }else{
       l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
     }
     
     model[t+1,"loglik"]         <- -as.numeric(opt$values[length(opt$values)])
     model[t+1,"predict_loglik"] <- l$Pred_loglik
     model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
     model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
     
     if(ModelType == 2){
       model[t+1,"Marg_loglik"]  <- l$Marg_loglik
       model[t+1,'AICm']         <- model[(t+1),"Marg_loglik"]-length(para)
       model[t+1,'BICm']         <- model[(t+1),"Marg_loglik"]-length(para)*log(nrow(ech))/2
     }
     model[t+1,]
    }
      close(pb)
      
    model[update_date[!(update_date==forecast.length)]+1,]<-Y
    
    strt <- (T-forecast.length) - window.size
    cat("Estimation step 2: \n")
    pb <- txtProgressBar(min=0, max = (forecast.length-1), style = 3)
    for(t in 0:(forecast.length-1)){
      ech    <- data[strt:(T-forecast.length+t),]
      
      if(t %in% update_date){
        para <- unlist(model[t+1,vars])
        if(refit.window == "moving") strt <- strt + refit.every
        next
      }
      
      model[t+1,vars] <- round(para,5)
      if(!(ModelType == 1)){
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
        
        pi_0 <- l$w_hat
        sig  <- volatilityVector(para=para,N=N,K=K)
        if(LEVIER){
          Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
          sig    <- sig*Levier
        }
        
        if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
          va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
          if(class(va) =='try-error') {
            va <- try(qmixnorm(VaR.alpha[iter], rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
            if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
          }else{
            model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
          }
        }
      }else{
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
      }
      
      model[t+1,"loglik"]         <- -l$loglik
      model[t+1,"predict_loglik"] <- l$Pred_loglik
      model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
      model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
      
      if(ModelType == 2){
        model[t+1,"Marg_loglik"]  <- l$Marg_loglik
        model[t+1,'AICm']         <- model[(t+1),"Marg_loglik"]-length(para)
        model[t+1,'BICm']         <- model[(t+1),"Marg_loglik"]-length(para)*log(nrow(ech))/2
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  }
  
  if((calculate.VaR) & (!(ModelType == 1))){
    for(iter in 1:length(VaR.alpha)){
      model[,paste0('I',100*(1-VaR.alpha[iter]))] <- (model[,'rt'] < model[,paste0('VaR',100*(1-VaR.alpha[iter]))])
    }
  }
  
  out<-list(N               = N,
            K               = K,
            ModelType       = ModelType,
            LEVIER          = LEVIER,
            n.ahead         = n.ahead,
            forecast.length = forecast.length, 
            refit.every     = refit.every, 
            refit.window    = refit.window,
            window.size     = window.size, 
            calculate.VaR   = calculate.VaR,
            VaR.alpha       = VaR.alpha,
            cluster         = cluster,
            data            = data,
            dates           = dates,
            estimates       = model)
  
  ### Prevision
  vars<-NULL
  if(!(ModelType==1)) vars <- c(vars,paste0("rt2p",1:n.ahead))
  if(!(ModelType==0)) vars <- c(vars,paste0("rvtp",1:n.ahead))
  model_add <- matrix(0, nrow=nrow(model_prev), ncol=length(vars))
  colnames(model_add) <- vars
  model_prev <- cbind(model_prev, model_add)
  
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)       vars <- c(vars,"l","theta")
  strt                  <- (T-forecast.length) - window.size
  
  if(is.null(cluster)){
    
    cat("Prevision step : \n")
    pb <- txtProgressBar(min=1, max = forecast.length, style = 3)
    
    for(t in 1:nrow(model)){
      ech    <- data[strt:(T-forecast.length+t-1),]
      
      para <- unlist(model[t, vars])
      l    <- logLik2(ech=ech, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N, t=nrow(ech))
      
      pi_0 <- l$w_hat
      if(t %in% update_date+1){
        if(refit.window == "moving") strt <- strt + refit.every
        sig  <- volatilityVector(para=para,K=K,N=N)
        matP <- P(para=para,K=K,N=N)
      }
    
      MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(n.ahead,n.bootpred)),n.ahead,n.bootpred,byrow=FALSE)) #simulation of Markov chain
      z_t    <- matrix(rnorm(n.bootpred*n.ahead),nrow=n.bootpred,ncol=n.ahead)
      
      if(LEVIER){
        Levier     <- rep(1,n.bootpred)%*%t(levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`Levier`)
        sim        <- R_hat( H          = n.ahead,
                             ech        = ech[((nrow(ech)-200):nrow(ech)),1],
                             MC_sim     = MC_sim,
                             z_t        = z_t,
                             Levier     = Levier,
                             sig        = sig,
                             para       = para,
                             Model_type = ModelType,
                             N          = N)
        if(!(ModelType==1)){
          rt2                                    <- sim$`rt2`
          model_prev[t,paste0("rt2p",1:n.ahead)] <- rt2[(ncol(rt2)-n.ahead+1):ncol(rt2)]
          
          if(ModelType==2){
            rvt                                    <- sim$`rvt`
            model_prev[t,paste0("rvtp",1:n.ahead)] <- rvt[(ncol(rvt)-n.ahead+1):ncol(rvt)]
          }
        }else{
          sim                                    <- sim$Levier
          model_prev[t,paste0("rvtp",1:n.ahead)] <- (f_sim(n.ahead,sig,pi_0,matP)$`rt2`)*(sim[(length(sim)-n.ahead+1):length(sim)])
        }
      }else{
        if(!(ModelType==2)){
          sim     <- f_sim(n.ahead,sig,pi_0,matP)
          if(ModelType==0) {
            model_prev[t,paste0("rt2p",1:n.ahead)]  <- sim$`rt2`
          }else{
            model_prev[t,paste0("rvtp",1:n.ahead)] <- sim$`rt2`
          }
        }else{
          xi      <- para[6];
          varphi  <- para[7];
          delta1  <- para[8];
          delta2  <- para[9];
          shape   <- para[10];
          sim     <- f_sim(n.ahead,sig,pi_0,matP,varphi,xi,shape,delta1,delta2)
          model_prev[t,paste0("rt2p",1:n.ahead)]  <- sim$`rt2`
          model_prev[t,paste0("rvtp",1:n.ahead)] <- sim$`rvt`
        }
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  }else{
    registerDoSNOW(cluster)
    cat("Prevision step : \n")
    pb <- txtProgressBar(min=1, max = forecast.length, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    model_prev <- foreach(t=1:nrow(model), .export=c("sim.mc"), .combine = "rbind",.options.snow=opts) %dopar% {
      if(refit.window == "moving") strt <- (T-forecast.length) - window.size + (t-1)*refit.every
      ech    <- data[strt:(T-forecast.length+t-1),]
              
      para <- unlist(model[t, vars])
      l<-logLik2(ech=ech, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N,t=nrow(ech))
      
      pi_0 <- l$w_hat
      if(t %in% update_date+1){
        if(refit.window == "moving") strt <- strt + refit.every
        sig  <- volatilityVector(para=para,K=K,N=N)
        matP <- P(para=para,K=K,N=N)
      }
      
      MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(n.ahead,n.bootpred)),n.ahead,n.bootpred,byrow=FALSE)) #simulation of Markov chain
      z_t<-matrix(rnorm(n.bootpred*n.ahead),nrow=n.bootpred,ncol=n.ahead)
      
      if(LEVIER){
        Levier     <- rep(1,n.bootpred)%*%t(levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`Levier`)
        sim        <- R_hat( H          = n.ahead,
                             ech        = ech[((nrow(ech)-200):nrow(ech)),1],
                             MC_sim     = MC_sim,
                             z_t        = z_t,
                             Levier     = Levier,
                             sig        = sig,
                             para       = para,
                             Model_type = ModelType,
                             N          = N)
        if(!(ModelType==1)){
          rt2                                   <- sim$`rt2`
          model_prev[t,paste0("rt2p",1:n.ahead)] <- rt2[(ncol(rt2)-n.ahead+1):ncol(rt2)]
          
          if(ModelType==2){
            rvt                                    <- sim$`rvt`
            model_prev[t,paste0("rvtp",1:n.ahead)] <- rvt[(ncol(rvt)-n.ahead+1):ncol(rvt)]
          }
        }else{
          sim                                    <- sim$Levier
          model_prev[t,paste0("rvtp",1:n.ahead)] <- (f_sim(n.ahead,sig,pi_0,matP)$`rt2`)*(sim[(length(sim)-n.ahead+1):length(sim)])
        }
      }else{
        if(!(ModelType==2)){
          sim     <- f_sim(n.ahead,sig,pi_0,matP)
          if(ModelType==0) {
            model_prev[t,paste0("rt2p",1:n.ahead)]  <- sim$`rt2`
          }else{
            model_prev[t,paste0("rvtp",1:n.ahead)] <- sim$`rt2`
          }
        }else{
          xi      <- para[6];
          varphi  <- para[7];
          delta1  <- para[8];
          delta2  <- para[9];
          shape   <- para[10];
          sim     <- f_sim(n.ahead,sig,pi_0,matP,varphi,xi,shape,delta1,delta2)
          model_prev[t,paste0("rt2p",1:n.ahead)]  <- sim$`rt2`
          model_prev[t,paste0("rvtp",1:n.ahead)] <- sim$`rvt`
        }
      }
      model_prev[t,]
    }
    close(pb)
  }
  
  out <- c(out,list(prevision = model_prev))
  
  class(out) <- "MDSVroll"
  
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

g<-function(vector){
  Sortie<-matrix(NA,2,2)
  for(i in c(0,1)) for(j in c(0,1)) Sortie[i+1,j+1]<-sum((vector[-1]==j)*(vector[-length(vector)]==i))
  if(vector[1]==0) Sortie[1,1]<-Sortie[1,1]+1
  if(vector[1]==1) Sortie[1,2]<-Sortie[1,2]+1
  colnames(Sortie)<-c(0,1)
  rownames(Sortie)<-c(0,1)
  return(Sortie)
}


#' @title Summarize and print MDSV Rolling estimates, volatility forecast and backtesting
#' @description Summary and print methods for the class \link{MDSVroll} as returned by the function \code{\link{MDSVroll}}.
#' @param object An object of class \link{MDSVroll}, output of the function \code{\link{MDSVroll}}.
#' @param x An object of class \link{summary.MDSVroll}, output of the function \code{\link{summary.MDSVroll}}
#' or class \link{MDSVroll} of the function \code{\link{MDSVroll}}.
#' @param VaR.test Whether to perform Value at Risk forecast backtesting.
#' @param Loss.horizon Horizon to summary the forecasts (cummulative and marginal).
#' @param Loss.window Window on which the forecasts are summarized.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item ModelType : type of models fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item n.ahead : integer designing the forecast horizon.
#'     \item forecast.length : length of the total forecast for which out of sample data from the dataset will be used for testing.
#'     \item refit.every : Determines every how many periods the model is re-estimated.
#'     \item refit.window : Whether the refit is done on an expanding window including all the previous data or a moving window where all previous 
#' data is used for the first estimation and then moved by a length equal to refit.every (unless the window.size option is used instead).
#'     \item window.size : If not NULL, determines the size of the moving window in the rolling estimation, which also determines the first point used.
#'     \item calculate.VaR : Whether to calculate forecast Value at Risk during the estimation.
#'     \item VaR.alpha : The Value at Risk tail level to calculate.
#'     \item cluster : A cluster object created by calling makeCluster from the parallel package.
#'     \item data : data use for the fitting.
#'     \item dates : vector or names of data designing the dates.
#'     \item estimates : matrix of all the parameters estimates at each date.
#'     \item prevision : matrix of all prevision made a each date.
#'     \item VaR.test : Whether to perform Value at Risk forecast backtesting.
#'     \item Loss.horizon : Horizon to summary the forecasts.
#'     \item Loss.window : Window on which the forecasts are summarized.
#'     \item Loss : Matrice containing the forecasts summary.
#' }
#' 
#' @details 
#' The \code{\link{summary.MDSVroll}} function compute the Root Mean Square Error, the Mean Average Error and the Quasi-Likehood 
#' error to summarize the forecasts. Those loss functions are compute for cummulative (by horizon) and marginal forecasts. 
#' For univariate realized variances model and joint log-returns and realized variances model, the loss functions are computed for
#' the realized variances and for the univariate log-returns model, the loss functions are computed for the log-returns.
#' For the Value-at-Risk basktest, the unconditionnal coverage test (see. Kupiec), the independance test (see Christoffersen) and the 
#' conditional coverage test (see Christoffersen and ) are performed.
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @export

"summary.MDSVroll" <- function(object, VaR.test=TRUE, Loss.horizon = c(1,5,10,25,50,75,100), Loss.window = 756, ...){
  stopifnot(class(object) == "MDSVroll")
  
  if(!is.logical(VaR.test)){
    stop("summary.MDSVroll(): input VaR.test must be logical!")
  }else if(VaR.test){
    if(!(object$calculate.VaR)){
      print("summary.MDSVroll() WARNING: Unable to perform VaR.test because object$calculate.VaR = FALSE! set VaR.test to FALSE")
      VaR.test<-FALSE
    }else if(object$ModelType == 1){
      print("summary.MDSVroll() WARNING: VaR is only for log-returns! set VaR.test to FALSE")
      VaR.test<-FALSE
    }
  } 
  
  if((!is.numeric(Loss.horizon)) || (!is.numeric(Loss.window))){
    stop("summary.MDSVroll(): input Loss.horizon and Loss.window must all be numeric!")
  }
  
  if(Loss.window > object$forecast.length){
    print("summary.MDSVroll() WARNING: Loss.window must be less than the forecast.length! set Loss.window to forecast.length")
    Loss.window <- object$forecast.length
  }
  
  ModelType    <- object$ModelType
  n.ahead      <- object$n.ahead
  if(sum(Loss.horizon>n.ahead)>0){
    Loss.horizon <- Loss.horizon[Loss.horizon<=n.ahead]
    print("summary.MDSVroll() WARNING: input Loss.horizon must be a less than n.ahead of the MDSVroll() object!")
  }
  
  Loss <- object$prevision[,1:7]
  if(ModelType == 2) Loss<-cbind(Loss,ModelType=object$prevision[,8])
  vars<-NULL
  if(!(ModelType==1)) vars<-c(vars,paste0("R_for_",Loss.horizon),paste0("R_tru_",Loss.horizon),
            paste0("R_for_m_",Loss.horizon),paste0("R_tru_m_",Loss.horizon))
  if(!(ModelType==0)) vars<-c(vars,paste0("RV_for_",Loss.horizon),paste0("RV_tru_",Loss.horizon),
                              paste0("RV_for_m_",Loss.horizon),paste0("RV_tru_m_",Loss.horizon))
  Loss_add <- matrix(0, nrow=nrow(Loss), ncol=length(vars))
  colnames(Loss_add) <- vars
  Loss <- cbind(Loss, Loss_add)
  
  if(!(ModelType==0)){
    for_RV_var <- rep(0,length(Loss.horizon))
    RV_var <- rep(0,length(Loss.horizon))
    
    for(t in 1:nrow(Loss)){
      rvt_sim <- object$prevision[t,paste0("rvtp",1:n.ahead)]
      for(k in 1:length(Loss.horizon)) for_RV_var[k] <- sum(rvt_sim[1:Loss.horizon[k]])
      names(for_RV_var) <- paste0("RV_for_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(for_RV_var)]<-for_RV_var
      for(k in 1:length(Loss.horizon)) for_RV_var[k] <- rvt_sim[Loss.horizon[k]]
      names(for_RV_var) <- paste0("RV_for_m_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(for_RV_var)]<-for_RV_var
      
      ech_rv <- object$prevision[t:nrow(Loss),"rvt"]
      for(k in 1:length(Loss.horizon)) RV_var[k] <- sum(ech_rv[1:Loss.horizon[k]])
      names(RV_var) <- paste0("RV_tru_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(RV_var)]<-RV_var
      for(k in 1:length(Loss.horizon)) RV_var[k] <- ech_rv[Loss.horizon[k]]
      names(RV_var) <- paste0("RV_tru_m_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(RV_var)]<-RV_var
    }
  }
  if(!(ModelType==1)){
    for_R_var <- rep(0,length(Loss.horizon))
    R_var <- rep(0,length(Loss.horizon))
    
    for(t in 1:nrow(Loss)){
      rt2_sim <- object$prevision[t,paste0("rt2p",1:n.ahead)]
      for(k in 1:length(Loss.horizon)) for_R_var[k] <- sum(rt2_sim[1:Loss.horizon[k]])
      names(for_R_var) <- paste0("R_for_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(for_R_var)]<-for_R_var
      for(k in 1:length(Loss.horizon)) for_R_var[k] <- rt2_sim[Loss.horizon[k]]
      names(for_R_var) <- paste0("R_for_m_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(for_R_var)]<-for_R_var
      
      ech_r <- object$prevision[t:nrow(Loss),"rt"]
      for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech_r[1:Loss.horizon[k]])^2)
      names(R_var) <- paste0("R_tru_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(R_var)]<-R_var
      for(k in 1:length(Loss.horizon)) R_var[k] <- (ech_r[Loss.horizon[k]])^2
      names(R_var) <- paste0("R_tru_m_",Loss.horizon)
      Loss[t,colnames(Loss) %in% names(R_var)]<-R_var
    }
  }

  out<-c(object,list(VaR.test     = VaR.test, 
                     Loss.horizon = Loss.horizon, 
                     Loss.window  = Loss.window, 
                     Loss         = Loss,
                     ...          = ...))
  
  class(out) <- "summary.MDSVroll"
  return(out)
}

#' @rdname summary.MDSVroll
#' @export
"print.summary.MDSVroll" <- function(x, ...){
  stopifnot(class(x) == "summary.MDSVroll")
  
  No.refit <- seq(0,x$forecast.length,by=x$refit.every)
  No.refit <- length(No.refit[!(No.refit == x$forecast.length)])
  out      <- x$Loss
  ind      <- (nrow(out)-x$Loss.window + 1):nrow(out)
  
  if(x$ModelType==0) Model_type <- "Univariate log-return"
  if(x$ModelType==1) Model_type <- "Univariate realized variances"
  if(x$ModelType==2) Model_type <- "Joint log-return and realized variances"
  
  cat("=================================================\n")
  cat(paste0("===  MDSV Rolling Estimation and Forecasting ====\n"))
  cat("=================================================\n\n")
  # cat("Conditional Variance Dynamique \n")
  # cat("------------------------------------------------- \n")
  cat(paste0("Model              : MDSV(",x$N,",",x$K,")\n"))
  cat(paste0("Data               : ",Model_type,"\n"))
  cat(paste0("Leverage           : ",x$LEVIER,"\n"))
  cat(paste0("No.refit           : ",No.refit,"\n"))
  cat(paste0("Refit Horizon      : ",x$refit.every,"\n"))
  cat(paste0("No.Forecasts       : ",x$forecast.length,"\n"))
  cat(paste0("n.ahead            : ",x$n.ahead,"\n"))
  cat(paste0("Date (T[0])        : ",x$dates[(nrow(x$data)-x$forecast.length)],"\n\n"))
  
  cat("Forecasting performances \n")
  cat("------------------------------------------------- \n")
  Pred_lik <- x$estimates[ind,grep('predict_loglik', colnames(x$estimates), fixed=TRUE)]
  cat(paste0("Predictive density : ",round(sum(Pred_lik),2),"\n"))
  cat("-------------------- \n\n")
  
  cat("Cummulative Loss Functions : \n")
  cat("---------------------------- \n")
  H_range        <- x$Loss.horizon
  if(!(x$ModelType == 1)){
    cat("Log-returns : \n")
    for_R_var      <- out[ind,grep('R_for', colnames(out), fixed=TRUE)]
    for_R_var      <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
    for_R_err      <- out[ind,grep('R_tru', colnames(out), fixed=TRUE)]
    for_R_err      <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
    if(length(H_range)==1){
      QLIK_R         <- mean(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
      RMSE_R         <- sqrt(mean( (for_R_var - for_R_err)^2, na.rm=TRUE ))/H_range
      MAE_R          <- mean( abs(for_R_var - for_R_err), na.rm=TRUE )/H_range
    }else{
      QLIK_R         <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
      RMSE_R         <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/H_range
      MAE_R          <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/H_range
    }
    
    
    Y              <- matrix(c(QLIK_R,RMSE_R,MAE_R),3,length(QLIK_R),T)
    row.names(Y)   <- c("QLIK","RMSE","MAE")
    colnames(Y)    <- H_range
    
    print(round(Y,3))
    x<-c(x,list(Loss_functions_R=Y))
  }
  if(!(x$ModelType == 0)){
    if(x$ModelType == 2) cat("\n")
    cat("Realized Variances : \n")
    for_RV_var     <- out[ind,grep('RV_for', colnames(out), fixed=TRUE)]
    for_RV_var      <- for_RV_var[,-grep('RV_for_m', colnames(for_RV_var), fixed=TRUE)]
    for_RV_err     <- out[ind,grep('RV_tru', colnames(out), fixed=TRUE)]
    for_RV_err      <- for_RV_err[,-grep('RV_tru_m', colnames(for_RV_err), fixed=TRUE)]
    if(length(H_range)==1){
      QLIK_RV        <- mean(log(for_RV_var) + for_RV_err/for_RV_var, na.rm=TRUE)
      RMSE_RV        <- sqrt(mean( (for_RV_var - for_RV_err)^2, na.rm=TRUE ))/H_range
      MAE_RV         <- mean( abs(for_RV_var - for_RV_err), na.rm=TRUE )/H_range
    }else{
      QLIK_RV        <- colMeans(log(for_RV_var) + for_RV_err/for_RV_var, na.rm=TRUE)
      RMSE_RV        <- sqrt(colMeans( (for_RV_var - for_RV_err)^2, na.rm=TRUE ))/H_range
      MAE_RV         <- colMeans( abs(for_RV_var - for_RV_err), na.rm=TRUE )/H_range
    }
    
    Y              <- matrix(c(QLIK_RV,RMSE_RV,MAE_RV),3,length(QLIK_RV),T)
    row.names(Y)   <- c("QLIK","RMSE","MAE")
    colnames(Y)    <- H_range
    
    print(round(Y,3))
    x<-c(x,list(Loss_functions_RV=Y))
  }
  
  cat(paste0("\n", "Marginal Loss Functions : \n"))
  cat("------------------------- \n")
  H_range        <- x$Loss.horizon
  if(!(x$ModelType == 1)){
    cat("Log-returns : \n")
    for_R_var      <- out[ind,grep('R_for_m', colnames(out), fixed=TRUE)]
    for_R_err      <- out[ind,grep('R_tru_m', colnames(out), fixed=TRUE)]
    if(length(H_range)==1){
      QLIK_R         <- mean(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
      RMSE_R         <- sqrt(mean( (for_R_var - for_R_err)^2, na.rm=TRUE ))/H_range
      MAE_R          <- mean( abs(for_R_var - for_R_err), na.rm=TRUE )/H_range
    }else{
      QLIK_R         <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
      RMSE_R         <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/H_range
      MAE_R          <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/H_range
    }
    
    
    Y              <- matrix(c(QLIK_R,RMSE_R,MAE_R),3,length(QLIK_R),T)
    row.names(Y)   <- c("QLIK","RMSE","MAE")
    colnames(Y)    <- H_range
    
    print(round(Y,3))
    x<-c(x,list(Loss_functions_marg_R=Y))
  }
  if(!(x$ModelType == 0)){
    if(x$ModelType == 2) cat("\n")
    cat("Realized Variances : \n")
    for_RV_var     <- out[ind,grep('RV_for_m', colnames(out), fixed=TRUE)]
    for_RV_err     <- out[ind,grep('RV_tru_m', colnames(out), fixed=TRUE)]
    if(length(H_range)==1){
      QLIK_RV        <- mean(log(for_RV_var) + for_RV_err/for_RV_var, na.rm=TRUE)
      RMSE_RV        <- sqrt(mean( (for_RV_var - for_RV_err)^2, na.rm=TRUE ))/H_range
      MAE_RV         <- mean( abs(for_RV_var - for_RV_err), na.rm=TRUE )/H_range
    }else{
      QLIK_RV        <- colMeans(log(for_RV_var) + for_RV_err/for_RV_var, na.rm=TRUE)
      RMSE_RV        <- sqrt(colMeans( (for_RV_var - for_RV_err)^2, na.rm=TRUE ))/H_range
      MAE_RV         <- colMeans( abs(for_RV_var - for_RV_err), na.rm=TRUE )/H_range
    }
    
    Y              <- matrix(c(QLIK_RV,RMSE_RV,MAE_RV),3,length(QLIK_RV),T)
    row.names(Y)   <- c("QLIK","RMSE","MAE")
    colnames(Y)    <- H_range
    
    print(round(Y,3))
    x<-c(x,list(Loss_functions_marg_RV=Y))
  }
  
  if(x$VaR.test){
    cat(paste0("\n","VaR Tests \n"))
    VaR.alpha <- x$VaR.alpha
    LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
    p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
    for(iter in 1:length(VaR.alpha)){
      cat("------------------------------------------------- \n")
      viol      <- sum(x$estimates[,paste0('I',100*(1-VaR.alpha[iter]))])
      alpha_hat <- sum(x$estimates[,paste0('I',100*(1-VaR.alpha[iter]))])/x$forecast.length
      LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(x$forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(x$forecast.length-viol))))
      decision  <- "No"
      if((1-pchisq(LR.uc,1) < VaR.alpha[iter])) decision  <- "Yes"
      
      cat(paste0("alpha              : ",VaR.alpha[iter],"%\n"))
      cat(paste0("Excepted Exceed    : ",round(VaR.alpha[iter]*x$forecast.length,1),"\n"))
      cat(paste0("Actual VaR Exceed  : ",viol,"\n"))
      cat(paste0("Actual %           : ",round(alpha_hat,2),"%\n\n"))
      
      cat("Unconditionnal Coverage (Kupiec)\n")
      cat("Null-Hypothesis    : Correct exceedances\n")
      cat(paste0("LR.uc Statistic    : ",round(LR.uc,3),"\n"))
      cat(paste0("LR.uc Critical     : ",round(qchisq(1-VaR.alpha[iter],1),3),"\n"))
      cat(paste0("LR.uc p-value      : ",round(1-pchisq(LR.uc,1),3),"\n"))
      cat(paste0("Reject Null        : ",decision,"\n\n"))
      
      viol      <- g(x$estimates[,paste0('I',100*(1-VaR.alpha[iter]))])
      pi        <- (viol[1,2]+viol[2,2])/sum(viol)
      pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
      pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
      LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                             ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
      decision  <- "No"
      if((1-pchisq(LR.ind,1) < VaR.alpha[iter])) decision  <- "Yes"
      
      cat("Independance (Christoffersen)\n")
      cat("Null-Hypothesis    : Independance of failures\n")
      cat(paste0("LR.ind Statistic   : ",round(LR.ind,3),"\n"))
      cat(paste0("LR.ind Critical    : ",round(qchisq(1-VaR.alpha[iter],1),3),"\n"))
      cat(paste0("LR.ind p-value     : ",round(1-pchisq(LR.ind,1),3),"\n"))
      cat(paste0("Reject Null        : ",decision,"\n\n"))
      
      decision  <- "No"
      if((1-pchisq(LR.ind+LR.uc,2) < VaR.alpha[iter])) decision  <- "Yes"
      
      cat("Conditionnal Coverage (Christoffersen)\n")
      cat("Null-Hypothesis    : Correct exceedances and Independance of failures\n")
      cat(paste0("LR.cc Statistic    : ",round(LR.ind+LR.uc,3),"\n"))
      cat(paste0("LR.cc Critical     : ",round(qchisq(1-VaR.alpha[iter],2),3),"\n"))
      cat(paste0("LR.cc p-value      : ",round(1-pchisq(LR.ind+LR.uc,2),3),"\n"))
      cat(paste0("Reject Null        : ",decision,"\n\n"))
      
      LR.uc_vec <- c(LR.uc_vec,LR.uc)
      LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
      LR.ind_vec <- c(LR.ind_vec,LR.ind)
      p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
      p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
      p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
    }
    names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")
    x<-c(x,list(LR.uc  = LR.uc_vec,
                LR.ind = LR.ind_vec,
                LR.cc  = LR.cc_vec,
                p.uc   = p.uc_vec,
                p.ind  = p.ind_vec,
                p.cc   = p.cc_vec))
  }
  
  invisible(x)
}


#' @rdname summary.MDSVroll
#' @export
"print.MDSVroll" <- function(x, ...) {
  stopifnot(class(x) == "MDSVroll")
  print(summary(x, ...))
}


#' @title Plot MDSV Rolling estimates, volatility forecast and backtesting
#' @description Plot methods for the class \link{MDSVroll} as returned by the function \code{\link{MDSVroll}}.
#' @param x An object of class \link{plot.MDSVroll}, output of the function \code{\link{plot.MDSVroll}}
#' or class \link{MDSVroll} of the function \code{\link{MDSVroll}}.
#' @param plot.type The type of plot to be draw.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return One or more of the following plots :
#' \itemize{
#'     \item sigma : Graph of the square of log-returns and/or realized variances with its 1-ahead forecasts.
#'     \item VaR : Graph showing the log-returns and the forecast conditional VaR at 1-ahead.
#'     \item dens : Graph of predicative density of log-retunrs or realized variances at 1-ahead.
#' }
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @export
"plot.MDSVroll" <- function(x, plot.type=c("sigma","VaR","dens"),...) {
  stopifnot(class(x) == "MDSVroll")
  
  if(prod(plot.type %in% c("sigma","VaR","dens"))==0){
    stop("summary.MDSVroll(): input plot.type must be sigma, VaR, or dens!")
  }
  
  if((x$ModelType == 1) & ("VaR" %in% plot.type)){
    print("summary.MDSVroll() WARNING: VaR is compute only for log-returns! remove VaR in plot.type")
    plot.type <- plot.type[!(plot.type == "VaR")]
  }
  
  x              <- c(x, list(plot.type = plot.type,
                              ...       = ...))
  
  class(x)       <- "plot.MDSVroll"
  x
}


#' @rdname plot.MDSVroll
#' @export
"print.plot.MDSVroll" <- function(x, ...) {
  
  ModelType <- x$ModelType
  Prev      <- x$prevision
  Estim     <- x$estimates
  plot.type <- x$plot.type
  
  if("sigma" %in% plot.type){
    if(!(ModelType==1)){
      ylim <- c(min((Prev[,"rt"])^2,Prev[,"rt2p1"])-0.25,
                max(c((Prev[,"rt"])^2,Prev[,"rt2p1"])+0.75))
      tmp           <- c(list(x = Prev[,"date"],
                              y = (Prev[,"rt"])^2,
                              type = "l",
                              main = "Log-returns square : 1.ahead forecast vs realized values",
                              ylab = "",
                              xlab = "Date",
                              col  = 'gray',
                              ylim = ylim,
                              ...  = ...), x[-(1:17)])
      do.call("plot", tmp)
      
      tmp           <- c(list(x = Prev[,"date"],
                              y = Prev[,"rt2p1"],
                              ylab = "",
                              col  = 'blue',
                              ...  = ...), x[-(1:17)])
      do.call("lines", tmp)
      legend("topleft", legend=c("1.ahead forecast", "realized values"), col=c("blue", "gray"), lty=1, cex=0.8)
    }
    
    if(!(ModelType==0)){
      ylim <- c(min(Prev[,"rvt"],Prev[,"rvtp1"])-0.25,
                max(c(Prev[,"rvt"],Prev[,"rvtp1"])+0.75))
      
      tmp           <- c(list(x = Prev[,"date"],
                              y = Prev[,"rvt"],
                              type = "l",
                              main = "Realized Variances : 1.ahead forecast vs realized values",
                              xlab = "Date",
                              ylab = "",
                              col  = 'gray',
                              ylim = ylim,
                              ...  = ...), x[-(1:17)])
      do.call("plot", tmp)
      
      tmp           <- c(list(x = Prev[,"date"],
                              y = Prev[,"rvtp1"],
                              col  = 'blue',
                              ...  = ...), x[-(1:17)])
      do.call("lines", tmp)
      
      legend("topleft", legend=c("1.ahead forecast", "realized values"), col=c("blue", "gray"), lty=1, cex=0.8)
    }
    
    par(mfrow=c(1,1))
  }
  
  if("VaR" %in% plot.type){
    for(iter in 1:length(x$VaR.alpha)){
    
      ylim <- c(min(c(Estim[,"rt"],Estim[,paste0("VaR",100*(1-x$VaR.alpha[iter]))]))-0.25,
                max(c(Estim[,"rt"],Estim[,paste0("VaR",100*(1-x$VaR.alpha[iter]))]))+0.75)
      
      tmp           <- c(list(x = Estim[,"date"],
                              y = Estim[,"rt"],
                              main = paste0("Log-returns and Value-at-Risk Exceedances (alpha = ",x$VaR.alpha[iter],")"),
                              type = "p",
                              xlab = "Date",
                              ylab = "",
                              col  = 'gray',
                              ylim = ylim,
                              pch = 16,
                              ...  = ...), x[-(1:17)])
      do.call("plot", tmp)
      
      tmp           <- c(list(x = Estim[,"date"],
                              y = Estim[,paste0("VaR",100*(1-x$VaR.alpha[iter]))],
                              col  = 'black',
                              ...  = ...), x[-(1:17)])
      do.call("lines", tmp)
      
      tmp           <- c(list(x = Estim[Estim[,paste0("I",100*(1-x$VaR.alpha[iter]))],"date"],
                              y = Estim[Estim[,paste0("I",100*(1-x$VaR.alpha[iter]))],"rt"],
                              col  = 'red', 
                              pch = 16,
                              ...  = ...), x[-(1:17)])
      do.call("points", tmp)
      
      legend('topleft',c('returns','return < VaR', 'VaR'),lty=c(NA,NA,1),pch=c(16,16,NA),col=c('gray','red','black'))
    }
    
  }
  
  if("dens" %in% plot.type){
    ind <- order(Estim[,"rt"])
    lo <- loess(Estim[ind,"predict_loglik"] ~ Estim[ind,"rt"])
    tmp           <- c(list(x = Estim[ind,"rt"],
                            y = predict(lo),
                            main = "Density forecasts",
                            type = "l",
                            xlab = "Log-returns",
                            ylab = "Densities",
                            ...  = ...), x[-(1:17)])
    do.call("plot", tmp)
  }
  
  invisible(x)
}



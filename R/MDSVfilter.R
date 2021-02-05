#' @export
MDSVfilter<-function(N,K,data,para,ModelType=0,LEVIER=FALSE){

  # N is the number of component for the MDSV process
  # K is the number of state of each MDSV process component
  # data is the data to be use for the MDSV estimation. It must be a numeric Txk matrix where
  #        T is the number of observations
  #        k = 1 if ModelType = 0 or (ModelType = 1 and LEVIER = FALSE)
  #        k = 2 if ModelType = 2 or (ModelType = 1 and LEVIER = TRUE).
  # ModelType is the type of model you want to fit.
  # ModelType=0 if you model an "Univariate log-return".
  # ModelType=1 if you model an "Univariate realized variances".
  # ModelType=2 if you model a "Joint log-return and realized variances".
  # LEVIER is a boolean to say if you want to estime a model with leverage or not.
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVfilter(): input N and K must be numeric!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVfilter(): input ModelType must be numeric!")
  }else if((ModelType>2) || (ModelType<0)){
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
      if(para[6]<=0){
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
  
  ### Results
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)       vars <- c(vars,"l","theta")
  
  names(para)<-vars
  
  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  if(N==1) para<-para[(vars[!(vars=='b')])]
  Np<-length(para)
  
  l<-logLik2(ech=data, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N)
  
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
            smoothed_proba = l$smoothed_proba,
            Pred_loglik    = l$Pred_loglik)
  
  if(ModelType==2) out<-c(out,list(Marg_loglik = l$Marg_loglik))
  
  if(!(ModelType==1)){
    pi_0 <- l$w_hat
    sig<- volatilityVector(para=para,N=N,K=K)
    if(LEVIER){
      Levier<-levierVolatility(ech=data[((T-200):T),1],para=para,Model_type=ModelType)$`levier`
      sig<-sig*Levier
    }
    
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

#' @importFrom graphics par
#' @export
"plot.MDSVfilter" <- function(x, ...) {
  stopifnot(class(x) == "MDSVfilter")
  
  if(x$ModelType == "Univariate log-return")                   ModelType <- 0
  if(x$ModelType == "Univariate realized variances")           ModelType <- 1
  if(x$ModelType == "Joint log-return and realized variances") ModelType <- 2
  
  para      <- x$estimates
  sig       <- sqrt(volatilityVector(para,K=x$K,N=x$K))
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
  
  .pardefault <- par()
  do.call("par", list(mar=c(6,6,4,4)))
  
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
  
  par(.pardefault)
  
  invisible(x)
}

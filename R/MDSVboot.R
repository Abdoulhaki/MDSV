#' @export
#' @importFrom mhsmm sim.mc
MDSVboot<-function(fit,n.ahead=100,n.bootpred=500,rtp1=0,rseed=NA){
  stopifnot(class(fit) == "MDSVfit")
  
  if ( (!is.numeric(n.ahead)) || (!is.numeric(n.bootpred)) ) {
    stop("MDSVboot(): inputs must be numeric!")
  }
  if (!is.numeric(rtp1)) {
    print("MDSVboot() WARNING: input rtp1 must be numeric! rtp1 set to 0")
    rtp1 <- 0;
  }
  if ((!is.numeric(rseed)) || is.na(rseed)) {
    print("MDSVboot() WARNING: input rseed must be numeric! rseed set to random")
    rseed <- sample.int(.Machine$integer.max,1)
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
  
  out<-c(ModelType = ModelType, fit[-1], list(dates      = dates,
                                              n.ahead    = n.ahead,
                                              n.bootpred = n.bootpred))
  
  l<-logLik2(ech=data, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N, t=T,r=rtp1)
  
  #simulation
  set.seed(rseed)
  
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
      if(Model_type==0) {
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


#' @export

"summary.MDSVboot" <- function(object, plot.type=NULL, ...){
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
  if(x$LEVIER){
    if(!(x$ModelType == 1)){
      rt_sim <- f.present(X = x$rt_sim, thr = n.ahead)
      rownames(rt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      cat("Log-returns (summary) : \n")
      print(head(round(rt_sim,6), min(n.ahead, 10)))
      if(n.ahead>10)  cat("......................... \n")
      
      if(x$ModelType == 0){
        rt2_sim <- f.present(X = (x$rt_sim)^2, thr = n.ahead)
        rownames(rt2_sim) <- paste("t", 1:n.ahead, sep="+")
        
        cat(paste0("\n","Sigma (summary) : \n"))
        print(head(round(rt2_sim,6), min(n.ahead, 10)))
        if(n.ahead>10)  cat("......................... \n")
      }
    }
    
    if(!(x$ModelType == 0)){
      rvt_sim <- f.present(X = x$rvt_sim, thr = n.ahead)
      rownames(rvt_sim) <- paste("t", 1:n.ahead, sep="+")
      
      if(x$ModelType == 2) cat("\n")
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

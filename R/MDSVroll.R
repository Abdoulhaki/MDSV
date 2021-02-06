#' @export
#' @importFrom mhsmm sim.mc
#' @importFrom Rsolnp solnp
#' @importFrom foreach foreach %dopar%
#' @import doParallel doSNOW
MDSVroll<-function(N, K, data, ModelType=0, LEVIER=FALSE, n.ahead = 1, forecast.length = 500, 
                   refit.every = 25, refit.window = "recursive", window.size = NULL, 
                   calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL, rseed = NA, ...){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVroll(): input N and K must be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVfit(): input N and K must be integer!")
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
  }else {
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

  if ( (!is.numeric(n.ahead)) || (!is.numeric(forecast.length)) || (!is.numeric(refit.every)) ) {
    stop("MDSVroll(): inputs n.ahead, forecast.length and refit.every must all be numeric!")
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
  }
  
  if ((!is.numeric(rseed)) || is.na(rseed)) {
    print("MDSVboot() WARNING: input rseed must be numeric! rseed set to random")
    rseed <- sample.int(.Machine$integer.max,1)
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
  ctrl <- list(TOL=1e-15, trace=0)
  para<-c(0.52,0.99, 2.77,sqrt(var(data[,1])),0.72)
  if(ModelType==1) para <- c(para,2.10)
  if(ModelType==2) para <- c(para,-0.5,	0.93,	-0.09,	0.04,	2.10)
  if(LEVIER)       para <- c(para,0.78,0.87568)
  
  para_tilde <- natWork(para=para,LEVIER=LEVIER,Model_type=ModelType)
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)       vars <- c(vars,"l","theta")
  
  
  ### Estimation
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- (T-forecast.length) - window.size
  opt<-NULL
  if(is.null(cluster)){
    cat("Estimation step : \n")
    
    pb <- txtProgressBar(min=0, max = forecast.length-1, style = 3)
    
    for(t in 0:(forecast.length-1)){
      ech    <- data[strt:(T-forecast.length+t),]
      
      if(t %in% update_date){
        if(!is.null(opt)) para_tilde<-opt$pars
        oldw <- getOption("warn")
        options(warn = -1)
        opt  <- try(solnp(pars=para_tilde,fun=logLik,ech=ech,Model_type=ModelType,K=K,
                       LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
        
        if (class(opt) =='try-error'){
          stop("MDSVroll(): Optimization ERROR")
        }
        options(warn = oldw)
        
        para <- workNat(opt$pars,LEVIER=LEVIER,Model_type=ModelType)
        
        if(refit.window == "moving") strt <- strt + refit.every
      }
      
      model[t+1,vars] <- round(para,5)
      if(!(Model_type == 1)){
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
        
        pi_0 <- l$w_hat
        sig  <- volatilityVector(para=para,N=N,K=K)
        if(LEVIER){
          Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
          sig    <- sig*Levier
        }
        
        for(iter in 1:length(VaR.alpha)){
          model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0)
        }
        
      }else{
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
      }
      
      model[t+1,"loglik"]         <- -as.numeric(opt$values[length(opt$values)])
      model[t+1,"predict_loglik"] <- l$Pred_loglik
      model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
      model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
      
      if(Model_type == 2){
        model[t+1,"Marg_loglik"]  <-l$Marg_loglik
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
    
      Y<- foreach(t=update_date[!(update_date==forecast.length)], .export=c("solnp"), 
                .packages =c("Rcpp","RcppArmadillo","RcppEigen","progressr"), .combine = rbind, .options.snow = opts) %dopar% { 

      if(refit.window == "moving") strt <- (T-forecast.length) - window.size + t*refit.every
      ech    <- data[strt:(T-forecast.length+t),]
      if(!is.null(opt)) para_tilde<-opt$pars
      oldw <- getOption("warn")
      options(warn = -1)
      opt  <- solnp(pars=para_tilde,fun=logLik,ech=ech,Model_type=ModelType,K=K,
                        LEVIER=LEVIER,N=N,Nl=70,control=ctrl)

      options(warn = oldw)
      
      para <- workNat(opt$pars,LEVIER=LEVIER,Model_type=ModelType)
       
     model[t+1,vars] <- round(para,5)
     if(!(Model_type == 1)){
       l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
       
       pi_0 <- l$w_hat
       sig  <- volatilityVector(para=para,N=N,K=K)
       if(LEVIER){
         Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
         sig    <- sig*Levier
       }
       
       for(iter in 1:length(VaR.alpha)){
         model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0)
       }
     }else{
       l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
     }
     
     model[t+1,"loglik"]         <- -as.numeric(opt$values[length(opt$values)])
     model[t+1,"predict_loglik"] <- l$Pred_loglik
     model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
     model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
     
     if(Model_type == 2){
       model[t+1,"Marg_loglik"]  <-l$Marg_loglik
       model[t+1,'AICm']         <- model[(t+1),"Marg_loglik"]-length(para)
       model[t+1,'BICm']         <- model[(t+1),"Marg_loglik"]-length(para)*log(nrow(ech))/2
     }
     model[t+1,]
    }
      close(pb)
      
    model[update_date[!(update_date==forecast.length)]+1,]<-Y
    
    cat("Estimation step 2: \n")
    pb <- txtProgressBar(min=0, max = (forecast.length-1), style = 3)
    for(t in 0:(forecast.length-1)){
      ech    <- data[1:(T-forecast.length+t),]
      
      if(t %in% update_date){
        para <- unlist(model[t+1,vars])
        next
      }
      
      model[t+1,vars] <- round(para,5)
      if(!(Model_type == 1)){
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rt"], Model_type = ModelType)
        
        pi_0 <- l$w_hat
        sig  <- volatilityVector(para=para,N=N,K=K)
        if(LEVIER){
          Levier <- levierVolatility(ech=ech[((nrow(ech)-200):nrow(ech)),1],para=para,Model_type=ModelType)$`levier`
          sig    <- sig*Levier
        }
        
        for(iter in 1:length(VaR.alpha)){
          model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0)
        }
      }else{
        l    <- logLik2(ech=ech,para=para,LEVIER=LEVIER,K=K,N=N,t=nrow(ech),r=model[t+1,"rvt"], Model_type = ModelType)
      }
      
      model[t+1,"loglik"]         <- -l$loglik
      model[t+1,"predict_loglik"] <- l$Pred_loglik
      model[t+1,'AIC']            <- model[(t+1),"loglik"]-length(para)
      model[t+1,'BIC']            <- model[(t+1),"loglik"]-length(para)*log(nrow(ech))/2
      
      if(Model_type == 2){
        model[t+1,"Marg_loglik"]  <- l$Marg_loglik
        model[t+1,'AICm']         <- model[(t+1),"Marg_loglik"]-length(para)
        model[t+1,'BICm']         <- model[(t+1),"Marg_loglik"]-length(para)*log(nrow(ech))/2
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  }
  
  if(!(ModelType == 1)){
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
  n.bootpred            <- 10000
  strt                  <- (T-forecast.length) - window.size
  
  if(is.null(cluster)){
    
    cat("Prevision step : \n")
    pb <- txtProgressBar(min=1, max = forecast.length, style = 3)
    
    for(t in 1:nrow(model)){
      ech    <- data[strt:(T-forecast.length+t-1),]
      
      para <- unlist(model[t, vars])
      l<-logLik2(ech=ech, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N)
      
      set.seed(rseed)
      
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
          if(Model_type==0) {
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
      l<-logLik2(ech=ech, para=para, Model_type=ModelType, LEVIER=LEVIER, K=K, N=N)
      
      set.seed(rseed)
      
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
          if(Model_type==0) {
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

#' @export

"summary.MDSVroll" <- function(object, VaR.test=TRUE, Loss.horizon = c(1,5,10,25,50,75,100), Loss.window = 756, ...){
  stopifnot(class(object) == "MDSVroll")
  
  if(!is.logical(VaR.test)){
    stop("summary.MDSVroll(): input VaR.test must be logical!")
  }else if(object$ModelType == 1){
    print("summary.MDSVroll() WARNING: VaR use only for log-returns! set VaR.test to FALSE")
    VaR.test<-FALSE
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
  if(!(ModelType==1)) vars<-c(vars,paste0("R_for_",Loss.horizon),paste0("R_tru_",Loss.horizon))
  if(!(ModelType==0)) vars<-c(vars,paste0("RV_for_",Loss.horizon),paste0("RV_tru_",Loss.horizon))
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
      
      ech_rv <- object$prevision[t:nrow(Loss),"rvt"]
      for(k in 1:length(Loss.horizon)) RV_var[k] <- sum(ech_rv[1:Loss.horizon[k]])
      names(RV_var) <- paste0("RV_tru_",Loss.horizon)
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
      
      ech_r <- object$prevision[t:nrow(Loss),"rt"]
      for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech_r[1:Loss.horizon[k]])^2)
      names(R_var) <- paste0("R_tru_",Loss.horizon)
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
  
  cat("Loss Functions :\n")
  cat("------------------- \n")
  H_range        <- x$Loss.horizon
  if(!(x$ModelType == 1)){
    cat("Log-returns : \n")
    for_R_var      <- out[ind,grep('R_for', colnames(out), fixed=TRUE)]
    for_R_err      <- out[ind,grep('R_tru', colnames(out), fixed=TRUE)]
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
  }
  if(!(x$ModelType == 0)){
    if(x$ModelType == 2) cat("\n")
    cat("Realized Variances : \n")
    for_RV_var     <- out[ind,grep('RV_for', colnames(out), fixed=TRUE)]
    for_RV_err     <- out[ind,grep('RV_tru', colnames(out), fixed=TRUE)]
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
  }
  
  if(x$VaR.test){
    cat(paste0("\n","VaR Tests \n"))
    VaR.alpha <- x$VaR.alpha
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
      pi        <- (viol[1,2]+viol[2,1])/sum(viol)
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
    }
    
  }
  
  invisible(x)
}


#' @rdname summary.MDSVroll
#' @export
"print.MDSVroll" <- function(x, ...) {
  stopifnot(class(x) == "MDSVroll")
  print(summary(x, ...))
}


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
    if(ModelType==2) par(mfrow=c(2,1))
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
      legend("topleft", legend=c("1.ahead forecast", "realized values"), col=c("red", "gray"), lty=1, cex=0.8)
    }
    
    if(!(ModelType==0)){
      ylim <- c(min((Prev[,"rvt"])^2,Prev[,"rvtp1"])-0.25,
                max(c((Prev[,"rvt"])^2,Prev[,"rvtp1"])+0.75))
      
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
      
      legend("topleft", legend=c("1.ahead forecast", "realized values"), col=c("red", "gray"), lty=1, cex=0.8)
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



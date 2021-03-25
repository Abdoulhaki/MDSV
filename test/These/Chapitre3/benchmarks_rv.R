path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/Returns"
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3"
setwd(path)
library(MDSV)
library(caTools)
library(Rcpp)
if(!require(rugarch)){install.packages("rugarch")}; library(rugarch)
if(!require(Rsolnp)){install.packages("Rsolnp")}; library(Rsolnp)
if(!require(parallel)){install.packages("parallel")}; library(parallel)
if(!require(doSNOW)){install.packages("doSNOW")}; library(doSNOW)
if(!require(KScorrect)){install.packages("KScorrect")}; library(KScorrect)
if(!require(mhsmm)){install.packages("mhsmm")}; library(mhsmm)
if(!require(Rcpp)){install.packages("Rcpp")}; library(Rcpp)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

index_set<-c("sp500", "nasdaq", "ftse", "n225")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

ctrl <- list(TOL=1e-15, trace=0)

source("benchmarks/realized_functions.r") 

#-------------------------------------------------------------------------------
# E S T I M A T I O N
#-------------------------------------------------------------------------------

### HAR

HAR <- expand.grid('index'=index_set,
                   'aggHAR'=c('RVar'),
                   'logHAR'=c(TRUE), 
                   'model'=c('HAR', 'HAR-lev'),
                   'RV'=c('RVar')
)
# HAR <- HAR[,ncol(HAR):1]
HAR <- cbind(HAR, 'v1'='mean', 'dist'='norm')

vars <- c('par', 'cst', 'RVd', 'RVw', 'RVm', 'b', 'Jd', 'RSd', 'Ld', 'Lw', 'Lm',
          'sig', 'loglik', 'AIC', 'BIC', 'n', 'conv')

HAR_add <- matrix(0, nrow=nrow(HAR), ncol=length(vars))
colnames(HAR_add) <- vars
HAR <- cbind(HAR, HAR_add)

for(i in 1:nrow(HAR)){
  index   <- as.character(HAR[i,'index'])
  RV_TYPE <- as.character(HAR[i,'RV'])
  MODEL   <- as.character(HAR[i,'model'])
  LOG     <- as.logical(HAR[i,'logHAR'])
  AGG     <- as.character(HAR[i,'aggHAR'])
  
  RV   <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_ <- rownames(RV)
  n    <- nrow(RV)
  
  x   <- RV[,"r"] - mean(RV[,"r"])
  Ld  <- c(0, x[-n]) #day
  Ld  <- pmin(Ld,0)
  Lw  <- runmean(c(rep(0, 5), x), k=5,
                 endrule='trim', align='right')[-(n+1)]  #week
  Lw  <- pmin(Lw,0)  
  Lm  <- runmean(c(rep(0, 22), x), k=22,
                 endrule='trim', align='right')[-(n+1)] #month
  Lm  <- pmin(Lm,0) 
  
  RV  <- log(RV[,"rv"])
  RVd <- c(mean(RV), RV[-n]) #day
  RVw <- runmean(c(rep(mean(RV), 5), RV), k=5,
                 endrule='trim', align='right')[-(n+1)]  #week
  RVm <- runmean(c(rep(mean(RV), 22), RV), k=22,
                 endrule='trim', align='right')[-(n+1)] #month
  
  if(MODEL=='HAR') {
    estim <- lm(RV ~ RVd + RVw + RVm)
  } else if(MODEL=='HAR-lev') {
    estim <- lm(RV ~ RVd + RVw + RVm + Ld + Lw + Lm)
  }
  params <- c(coef(estim), 'sig'=summary(estim)$sigma)
  names(params)[1] <- 'cst'
  HAR[i,colnames(HAR) %in% names(params)] <- params
  params_n <- length(params)  
  
  #log-likelihood
  HAR[i,'loglik'] <- round(sum( log( dnorm(residuals(estim), sd=params['sig']) ) - RV ), 3)   
  HAR[i,'AIC'] <- round(-length(params) + HAR[i,'loglik'], 3)
  HAR[i,'BIC'] <- round(-0.5*log(length(RV))*length(params) + HAR[i,'loglik'], 3)
  
  HAR[i,'par'] <- params_n
  HAR[i,'b']   <- 0
  HAR[i,'n']   <- length(RV)
  HAR[i,'conv']<- 0
  
}

#### MEM and AMEM

MEM <- function(all_vars,
                Vt=list(model = c('MEM', 'AMEM'),
                        RV = c('RVol', 'RVar'),
                        logHAR = c(FALSE, TRUE),
                        aggHAR = c('-', 'RVar', 'AR'),
                        mult = 100^2 #for realized variance
                ),
                dist = c('gamma'),
                init = list(method=c('mean','estim','unc'), v1=NA),#not used for HAR models
                Data){
  
  params <- working_to_real_MEM(all_vars, Vt, dist, init)
  w <- params$w; ad <- params$ad; aw <- params$aw; am <- params$am;
  b <- params$b; g1 <- params$g1; g2 <- params$g2; 
  Ld <- params$Ld; Lw <- params$Lw; Lm <- params$Lm; 
  o  <- params$o; v1 <- params$v1
  
  mult   <- Vt$mult
  logHAR <- Vt$logHAR
  
  R  <- Data[,'r'] - mean(Data[,'r'])         #returns
  RV <- Data[,'rv'] #realized variance
  n  <- length(RV)
  
  if(Vt$RV=='RVol') RV <- sqrt(RV)
  
  #v1 <- mean(RV)
  e  <- w + ( ad + g1*(R[-n]<0) ) * RV[-n]
  vt <- c( v1, filter(e, b, "r", init=v1) )
  
  loglik <- sum( log( dgamma(RV/vt, shape=o, scale=1/o) / vt ) )
  
  loglik <- -loglik
  attr(loglik, "Vt") <- vt
  return(loglik)
  
}


design <- expand.grid(  'index'=index_set,
                        'v1'=c('mean'),
                        'dist' = c('gamma'),
                        'aggHAR' = c('-'),
                        'model' = c('MEM', 'AMEM'),
                        'logHAR' = c(FALSE), 
                        'RV' = c('RVar'))
design_mat <- design[,ncol(design):1]                          

res_mat <- NULL
for(i in 1:nrow(design_mat)){   
  index   <- as.character(design_mat[i,'index'])
  
  RV   <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_ <- rownames(RV)
  n    <- nrow(RV)
  
  Vt <- list(  'model' = as.vector(design_mat$model[i]),
               'RV' = as.vector(design_mat$RV[i]),
               'logHAR' = as.vector(design_mat$logHAR[i]),
               'aggHAR' = as.vector(design_mat$aggHAR[i])
  )
  dist <- as.vector(design_mat$dist[i])
  
  w  <- 0.05
  ad <- 0.35
  aw <- 0
  am <- 0
  b  <- 0.60
  g1 <- 0.05
  g2 <- 0
  
  Ld <- 0
  Lw <- 0
  Lm <- 0
  
  o  <- 2
  
  v1 <- RV[,"rv"]
  v1 <- mean(v1)
  init <- list('method' = as.vector("mean"),
               'v1' = v1
  )
  
  all_vars <- real_to_working_MEM(w, ad, aw, am, b, g1, g2, Ld, Lw, Lm, o,
                                  Vt, dist, init)
  # working_to_real_MEM(all_vars, Vt, dist, init)
  
  loglik0 <- -MEM(all_vars, Vt, dist, init, RV)
  
  start.time <- proc.time()["elapsed"]
  Optim <- try(solnp(pars=all_vars, fun=MEM, Vt=Vt, dist=dist, init=init,
                     Data=RV, control=ctrl), silent=TRUE)
  end.time <- proc.time()["elapsed"]
  run.time <- as.vector(end.time - start.time)  
  
  if(length(attr(Optim,'class'))==0){#no error
    
    params <- working_to_real_MEM(Optim$pars, Vt, dist, init)
    loglik <- -Optim$values[length(Optim$values)]
    
    params_n <- length(Optim$pars)
    
    res <- c('class'='MEM',
             'RV'=Vt$RV,
             'index'=index,
             'model'=Vt$model,
             'logHAR'=Vt$logHAR,
             'aggHAR'=Vt$aggHAR,
             'v1'=init$method,
             'dist'=dist,
             'par'=params_n,
             unlist(list(cst=params$w, RVd=params$ad, RVw=params$aw, RVm=params$am,
                         b=params$b, Jd=params$g1, RSd=params$g2, 
                         Ld=params$Ld, Lw=params$Lw, Lm=params$Lm, sig=params$o)),
             #'vol_unc'=round(unc,3),
             #'pers'=round(pers, 4),
             'loglik'=round(loglik, 3),
             'AIC'=round(-params_n + loglik, 3),
             'BIC'=round(-0.5*log(n)*params_n + loglik, 3),         
             'n'=n,
             'mult'=Vt$mult,
             'conv'=Optim$convergence
    )
    
  } else if(attr(Optim,'class')=='try-error'){ #ERROR IN OPTIMIZATION
    
    res <- res #we take old value
    res[1:length(res)] <- NA
    res['RV'] <- Vt$RV
    res['model'] <- Vt$model
    res['logHAR'] <- Vt$logHAR
    res['aggHAR'] <- Vt$aggHAR
    res['v1'] <- init$method
    res['dist'] <- dist
  }
  
  res_mat <- rbind(res_mat, res)

}

filename <- paste("Estim_HAR_MEM_all_", start.date, "_", end.date, sep="")
write.csv(rbind(HAR,res_mat[,names(HAR)]), paste(filename,"csv",sep="."), row.names=FALSE)



#### FHMV

#some functions
para_names<-function(LEVIER){
  vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q","shape")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

para_init<-function(LEVIER){
#c(0.4908,4.01718,0.35581,0.99615,8.42876,0.8,0.83108)
  # para<- c(2.43228,	1.73923,	0.99911,	0.99234,	15.89259,	0.93064,	0.92864, 3.25)
  para<- c(1.2146,	2.4756,	0.8975,	0.9857,	3.6969,	0.7180,	0.1229,	6.8726)

  if(LEVIER) para<-c(para,0.44768,	0.8180)
  return(para)
}

ctrl <- list(TOL=1e-15, trace=0)

model<-expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   Model = c("FHMV"), LEVIER=c(FALSE,TRUE))

vars<-c('loglik', 'AIC', 'BIC', c("sigma","c1","theta_c","p","m1","theta_m","q","shape","l","theta"),"time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

#setup parallel backend to use many processors
cluster <- makeCluster(6)
registerDoSNOW(cluster)
pb <- txtProgressBar(min=1, max = nrow(model), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

Y<-foreach(i=1:nrow(model), .combine=rbind, .export=c("solnp"), 
           .packages = c("MDSV","Rcpp","RcppArmadillo","RcppEigen"),
           .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
             start_time <- Sys.time()
             index                 <- as.character(model[i,"index"])
             donne                 <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
             nam_                  <- rownames(donne)
             donne[,"r"]           <- donne[,"r"] - mean(donne[,"r"])
             donne                 <- cbind(donne[,"rv"],donne[,"r"])
             model[i,"start.date"] <- as.Date(nam_[1])
             model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
             model[i,"length"]     <- length(nam_)
             Model                 <- as.character(model[i,"Model"])
             LEVIER                <- model[i,"LEVIER"]
             N                     <- 10
             
             sourceCpp(paste0("benchmarks/FHMV_rv.cpp"))
             
             para<-para_init(LEVIER)
             para_tilde<-natWork(para,LEVIER)
             opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
             vars<-para_names(LEVIER)
             params<-workNat(opt$pars,LEVIER)
             names(params)<-para_names(LEVIER)
             
             model[i,'loglik']    <- -as.numeric(opt$values[length(opt$values)])
             model[i,'AIC']       <- model[i,"loglik"] - length(params)
             model[i,'BIC']       <- model[i,"loglik"]-(length(params))*log(length(nam_))/2
             model[i,vars]        <- params[vars]
             model[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
             
             model[i,]
           }

close(pb)
stopCluster(cluster)

write.csv(Y, "Estim_FHMV_RV.csv", row.names=FALSE)


#-------------------------------------------------------------------------------
# P R E V I S I O N
#-------------------------------------------------------------------------------

### HAR

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050

#--------------------
#HAR regression model
#--------------------
design <- expand.grid(        'v1'='mean',
                              'init.params'='fixed',
                              'dist'='norm',
                              'model'=c('HAR', 'HAR-lev'),
                              'RV'='RVar',
                              'class'='HAR',
                              'index'=index_set
)
design_mat <- design[,ncol(design):1]  

#--------------------
#MEM specification
#--------------------
#MEM
design <- expand.grid(        'v1'='mean',
                              'init.params'='fixed',
                              'dist'=c('gamma'),
                              'model'=c('MEM','AMEM'),
                              'RV'='RVar',
                              'class'='MEM',
                              'index'=index_set
)
design_mat <- rbind(design_mat,design[,ncol(design):1]) 


vars<-c('start.date','end.date','pred_dens',paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon),"time")

design_add <- matrix(0, nrow=nrow(design_mat), ncol=length(vars))
colnames(design_add) <- vars
design_mat <- cbind(design_mat, design_add)

for( i_model in 1:nrow(design_mat) ){   
  start_time <- Sys.time()
  res_mat <- NULL
  
  class_    <- as.vector(design_mat$class[i_model])
  index     <- as.vector(design_mat$index[i_model])
  RV_TYPE   <- as.vector(design_mat$RV[i_model])
  model     <- as.vector(design_mat$model[i_model])
  dist      <- as.vector(design_mat$dist[i_model])
  init.params <- as.vector(design_mat$init.params[i_model])
  v1_ <- as.vector(design_mat$v1[i_model])
  
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- 1
  
  RV   <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_ <- rownames(RV)
  n    <- nrow(RV)
  
  if(class_=='HAR'){
    
    R_   <- RV[,"r"] - mean(RV[,"r"])
    Ld_  <- c(0, R_[-n]) #day
    Ld_  <- pmin(Ld_,0)
    Lw_  <- runmean(c(rep(0, 5), R_), k=5,
                    endrule='trim', align='right')[-(n+1)]  #week
    Lw_  <- pmin(Lw_,0)  
    Lm_  <- runmean(c(rep(0, 22), R_), k=22,
                    endrule='trim', align='right')[-(n+1)] #month
    Lm_  <- pmin(Lm_,0) 
    
    RV_  <- log(RV[,"rv"])
    RVd_ <- c(mean(RV_), RV_[-n]) #day
    RVw_ <- runmean(c(rep(mean(RV_), 5), RV_), k=5,
                    endrule='trim', align='right')[-(n+1)]  #week
    RVm_ <- runmean(c(rep(mean(RV_), 22), RV_), k=22,
                    endrule='trim', align='right')[-(n+1)] #month
    
    for(t in 0:(forecast.length-1)){
      RVout <- exp(RV_[(n-forecast.length+t+1):n])
      Rout <- R_[(n-forecast.length+t+1):n]
      RV   <- RV_[strt:(n-forecast.length+t)]
      R    <- R_[strt:(n-forecast.length+t)]
      RVd  <- RVd_[strt:(n-forecast.length+t)]
      RVw  <- RVw_[strt:(n-forecast.length+t)]
      RVm  <- RVm_[strt:(n-forecast.length+t)]
      Ld   <- Ld_[strt:(n-forecast.length+t)]
      Lw   <- Lw_[strt:(n-forecast.length+t)]
      Lm   <- Lm_[strt:(n-forecast.length+t)]
      
      if(t %in% update_date){
        
        if(model=='HAR') {
          estim <- lm(RV ~ RVd + RVw + RVm)
        } else if(model=='HAR-lev') {
          estim <- lm(RV ~ RVd + RVw + RVm + Ld + Lw + Lm)
        }
        params <- c(coef(estim), 'sig'=summary(estim)$sigma)
        names(params)[1] <- 'cst'
        params_n <- length(params)  
        
        if(refit.window == "moving") strt <- strt + refit.every
        
        resids <- residuals(estim)
      }else{
        if(model=='HAR') {
          resids <- RV-cbind(1,RVd,RVw,RVm)%*%(params[!(names(params)==c('sig'))])
        } else if(model=='HAR-lev') {
          resids <- RV-cbind(1,RVd,RVw,RVm,Ld,Lw,Lm)%*%(params[!(names(params)==c('sig'))])
        }
      }
      
      loglik <- round(sum( log( dnorm(resids, sd=params['sig']) ) - RV ), 3)   
      
      #----------------------------------
      #MULTI-PERIOD VARIANCE FORECASTS
      #----------------------------------
      cst <- as.vector(params['cst'])
      RVd <- as.vector(params['RVd'])
      RVw <- as.vector(params['RVw'])
      RVm <- as.vector(params['RVm'])
      Ld <- as.vector(params['Ld'])
      Lw <- as.vector(params['Lw'])
      Lm <- as.vector(params['Lm'])
      sig <- as.vector(params['sig'])
      
      if(model=='HAR'){
        #AR coefficients of HAR model
        AR <- as.vector(c(RVd+RVw/5+RVm/22, rep(RVw/5+RVm/22,4), rep(RVm/22,17))) 
        RV_st <- RV[length(RV):(length(RV)-21)] #final 22 values
        
        #point forecasts over maximal horizon
        var_point <- as.vector(filter(rep(cst,max(Loss.horizon)), AR, "r", init=RV_st))
      }else if(model=='HAR-lev'){
        
        AR <- as.vector(c(RVd+RVw/5+RVm/22, rep(RVw/5+RVm/22,4), rep(RVm/22,17))) 
        RV_st <- RV[length(RV):(length(RV)-21)] #final 22 values
        r_st <- R #final 22 values
        r_st <- r_st[length(r_st):(length(r_st)-21)] #final 22 values
        rd_ <- r_st[22] #day
        rd_ <- min(rd_,0)
        rw_ <- mean(r_st[18:22])  #week
        rw_ <- pmin(rw_,0)  
        rm_ <- mean(r_st) #month
        rm_ <- pmin(rm_,0) 
        
        #AR coefficients of HAR model
        var_point <- numeric(max(Loss.horizon))
        for(hor_i in 1:max(Loss.horizon)){
          var_point[hor_i] <- cst + sum(RV_st * AR) + Ld*rd_ + Lw*rw_ + Ld*rm_
          RV_st <- c(RV_st[-1],var_point[hor_i])
          r_st <- c(r_st[-1],rnorm(1))
          
          rd_ <- r_st[22] #day
          rd_ <- min(rd_,0)
          rw_ <- mean(r_st[18:22])  #week
          rw_ <- pmin(rw_,0)  
          rm_ <- mean(r_st) #month
          rm_ <- pmin(rm_,0) 
        }
      }
      
      var_point<-exp(var_point)
      
      #aggregated point forecasts
      for_var   <- numeric(length(Loss.horizon))
      for(i in 1:length(Loss.horizon)) for_var[i] <- sum(var_point[1:Loss.horizon[i]])  
      names(for_var) <- paste('var_RV_for',Loss.horizon,sep='_')      
      
      #----------------------------------
      #OUT-OF-SAMPLE DATA
      #----------------------------------
      
      RV <- RVout

      #-----------------------------------------
      #ONE-PERIOD PROBABILITY INTEGRAL TRANSFORM
      #CORSI ET AL. (2008) EQUATION (33)
      #-----------------------------------------
      for_dens <- pnorm(RV[1], mean=var_point[1], sd=sig)
      
      #-----------------------------------------
      #ONE-PERIOD SCORE
      #CORSI ET AL. (2008) EQUATION (34)
      #-----------------------------------------
      for_scr <- dnorm(RV[1], mean=var_point[1], sd=sig, log=TRUE)
      
      #----------------------------------
      #RV
      #----------------------------------
      for_err_sq <- numeric(length(Loss.horizon))
      for(i in 1:length(Loss.horizon)) for_err_sq[i] <- sum(RV[1:Loss.horizon[i]])    
      names(for_err_sq) <- paste('var_RV_tru',Loss.horizon,sep='_')
      
      res <- c('date'=nam_[n-forecast.length+t+1],
               'class'=class_,
               'RV'=RV_TYPE,
               'model'=model,
               'dist'=dist,
               'v1'=v1_,
               'for_dens'=as.vector(for_dens),
               'predict_loglik'=as.vector(for_scr),
               for_var,
               for_err_sq,
               unlist(list(cst=cst, RVd=RVd, RVw=RVw, RVm=RVm, sig=sig)),
               'loglik'=round(loglik, 3),
               'n'=n-forecast.length+t,
               'conv'=0
      )
      
      res_mat <- rbind(res_mat, res)
      
      if(model=="HAR") {LEVIER<-FALSE}else if(model=="HAR-lev"){LEVIER<-TRUE}
      filename <- paste("Forecast_HAR_LEVIER_",LEVIER,"_",index, "_length_",
                        forecast.length, start.date, "_", end.date, sep="")
      
      if( t %in% update_date ) write.csv(res_mat, paste(filename,"csv",sep="."), row.names=FALSE)
    }

  } else if(class_=='MEM'){
    
    R_   <- RV[,"r"] - mean(RV[,"r"])
    RV_  <- RV[,"rv"]
    
    Vt <- list(  'model' = model,
                 'RV' = RV_TYPE,
                 'logHAR' = FALSE,
                 'aggHAR' = '-',
                 'mult' = 100^2
    )
    w  <- 0.05
    ad <- 0.35
    aw <- 0
    am <- 0
    b  <- 0.60
    g1 <- 0.05
    g2 <- 0
    Ld <- 0
    Lw <- 0
    Lm <- 0
    o <- 2 #parameter for error distribution
    
    for(t in 0:(forecast.length-1)){
      RVout <- RV_[(n-forecast.length+t+1):n]
      Rout <- R_[(n-forecast.length+t+1):n]
      RV   <- RV_[strt:(n-forecast.length+t)]
      R    <- R_[strt:(n-forecast.length+t)]

      init <- list('method' = v1_, 'v1' = mean(RV))
      
      all_vars <- real_to_working_MEM(w, ad, aw, am, b, g1, g2, Ld, Lw, Lm, o,
                                      Vt, dist, init)
      # working_to_real_MEM(all_vars, Vt, dist, init)
      
      if(t %in% update_date){
        Optim <- solnp(pars=all_vars, fun=MEM, Vt=Vt, dist=dist, init=init,
                       Data=cbind("rv"=RV,"r"=R), control=ctrl)
      }
      
      loglik <- -MEM(Optim$pars, Vt, dist, init, Data=cbind("rv"=RV,"r"=R))
      params <- working_to_real_MEM(Optim$pars, Vt, dist, init)
      params_n <- length(Optim$pars) #+ (mean_incl==FALSE)
      
      cst <- params$w
      RVd <- params$ad
      RVw <- params$aw
      RVm <- params$am
      b   <- params$b
      Jd  <- params$g1
      o <- sig <- params$o
      
      #----------------------------------
      #MULTI-PERIOD VARIANCE FORECASTS
      #----------------------------------
      
      vt <- attr(MEM(Optim$pars, Vt, dist, init, Data=cbind("rv"=RV,"r"=R)), 'Vt')
      vT <- vt[length(vt)] #T forecast
      vT_1 <- cst + ( RVd + Jd*(R[length(R)]<0) ) *
        RV[length(RV)] + b * vT  #T+1 forecast
      
      #PERSISTENCE        
      pers <- RVd + b + Jd/2
      #UNCONDITIONAL FORECAST
      unc  <- cst / (1-pers)
      
      #point forecasts over maximal horizon
      var_point <- unc + (vT_1 - unc) * pers^(1:max(Loss.horizon)-1)
      
      #--------------------------
      #aggregated point forecasts
      #--------------------------
      for_var   <- numeric(length(Loss.horizon))
      for(i in 1:length(Loss.horizon)) for_var[i] <- sum(var_point[1:Loss.horizon[i]])  
      names(for_var) <- paste('var_RV_for',Loss.horizon,sep='_')      
      
      #----------------------------------
      #OUT-OF-SAMPLE DATA
      #----------------------------------
      
      RV <-RVout
      
      #-----------------------------------------
      #ONE-PERIOD PROBABILITY INTEGRAL TRANSFORM
      #CORSI ET AL. (2008) EQUATION (33)
      #-----------------------------------------
      
      for_dens <- pgamma(RV[1]/var_point[1], shape=o, scale=1/o)
      
      #-----------------------------------------
      #ONE-PERIOD SCORE
      #CORSI ET AL. (2008) EQUATION (34)
      #-----------------------------------------

      for_scr <- log( dgamma(RV[1]/var_point[1], shape=o, scale=1/o) / var_point[1] )
      
      #----------------------------------
      #RV
      #----------------------------------
      for_err_sq <- numeric(length(Loss.horizon))
      for(i in 1:length(Loss.horizon)) for_err_sq[i] <- sum(RV[1:Loss.horizon[i]])    
      names(for_err_sq) <- paste('var_RV_tru',Loss.horizon,sep='_')
      
      res <- c('date'=nam_[n-forecast.length+t+1],
               'class'=class_,
               'RV'=RV_TYPE,
               'model'=model,
               'dist'=dist,
               'v1'=v1_,
               'for_dens'=as.vector(for_dens),
               'predict_loglik'=as.vector(for_scr),
               for_var,
               for_err_sq,
               unlist(list(cst=cst, RVd=RVd, RVw=RVw, RVm=RVm, b=b, Jd=Jd, sig=sig)),
               'unc'=round(sqrt(unc*252), 3),
               'pers'=round(pers, 4),
               'loglik'=round(loglik, 3),
               'n'=n-forecast.length+t,
               'conv'=0
      )
      
      res_mat <- rbind(res_mat, res)
      
      if(model=="MEM") {LEVIER<-FALSE}else if(model=="AMEM"){LEVIER<-TRUE}
      filename <- paste("Forecast_MEM_LEVIER_",LEVIER,"_",index, "_length_",
                        forecast.length, start.date, "_", end.date, sep="")
      
      if( t %in% update_date ) write.csv(res_mat, paste(filename,"csv",sep="."), row.names=FALSE)
    }
  }
  
  write.csv(res_mat, paste(filename,"csv",sep="."), row.names=FALSE)
  
  design_mat[i_model,'pred_dens'] <- sum(as.numeric(res_mat[,"predict_loglik"]))
  
  for_R_var <- matrix(as.numeric(res_mat[,grep('var_RV_for', colnames(res_mat), fixed=TRUE)]),ncol=length(Loss.horizon))
  for_R_err <- matrix(as.numeric(res_mat[,grep('var_RV_tru', colnames(res_mat), fixed=TRUE)]),ncol=length(Loss.horizon))
  design_mat[i_model,paste('QLIKc_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  design_mat[i_model,paste('RMSEc_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  design_mat[i_model,paste('MAEc_RV',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  design_mat[i_model,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
}

write.csv(design_mat, "Forecast_HAR_MEM_All_756.csv", row.names=FALSE)

#### FHMV 

#some functions
para_names<-function(LEVIER){
  vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q","shape")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

para_init<-function(LEVIER){
  para<- c(1.2146,	2.4756,	0.8975,	0.9857,	3.6969,	0.7180,	0.1229,	6.8726)
  if(LEVIER) para<-c(para,0.44768,	0.8180)
  return(para)
}

ctrl <- list(TOL=1e-15, trace=0)

model_extern<-expand.grid(index=index_set, Model = c('FHMV'), LEVIER=c(TRUE,FALSE), n.ahead=c(100),
                          forecast.length=c(756), refit.every=c(63), refit.window=c("recursive"), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
rseed           <- 1050
cluster         <- makeCluster(6)

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

R_var <- rep(0,length(Loss.horizon))

for(i in 1:nrow(model_extern)){
  start_time <- Sys.time()
  index                        <- as.character(model_extern[i,"index"])
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne[,"r"]                  <- donne[,"r"] - mean(donne[,"r"])
  donne                        <- cbind("rv"=donne[,"rv"],"r"=donne[,"r"])
  model_extern[i,"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i,"end.date"]   <- as.character(nam_[length(nam_)])
  Model                        <- as.character(model_extern[i,"Model"])
  LEVIER                       <- model_extern[i,"LEVIER"]
  N                            <- 10
  
  sourceCpp(paste0("benchmarks/FHMV_rv.cpp"))
  set.seed(rseed)
  
  filename <- paste("Forecast_FHMVrv_LEVIER_",LEVIER,"_", index, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  model<-expand.grid(date = nam_[((length(nam_)-forecast.length+1):length(nam_))],rvt=0,Levier=LEVIER)
  model$rvt  <- donne[((length(nam_)-forecast.length+1):length(nam_)),"rv"]
  
  vars<-c(paste0("RV_for",Loss.horizon),paste0("RV_tru",Loss.horizon),'pred_lik','loglik', 'AIC', 'BIC',
          "sigma","c1","theta_c","p","m1","theta_m","q","shape","l","theta")
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  para<-para_init(LEVIER)
  vars<-para_names(LEVIER)
  para_tilde<-natWork(para,LEVIER)
  
  opt<-NULL
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- 1
  
  registerDoSNOW(cluster)
  cat("Estimation step 1: \n")
  
  pb <- txtProgressBar(min= 0, max = length(update_date[!(update_date==forecast.length)]), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  Y<- foreach(t=update_date[!(update_date==forecast.length)], .packages = c("Rcpp","RcppArmadillo","RcppEigen"), 
              .export=c("solnp"), .noexport=c("workNat", "logLik", "logLik2","volatilityVector", "levierVolatility"),
              .combine = rbind, .options.snow = opts) %dopar% {
                # for(t in update_date[!(update_date==forecast.length)]){              
                ech    <- donne[strt:(length(nam_)-forecast.length+t),]
                
                sourceCpp(paste0("benchmarks/FHMV_rv.cpp"))
                
                if(!is.null(opt)) para_tilde<-opt$pars
                oldw <- getOption("warn")
                options(warn = -1)
                opt  <- try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
                options(warn = oldw)
                
                para <- workNat(opt$pars,LEVIER=LEVIER)
                model[t+1,vars] <- round(para,5)
                
                l    <- logLik2(ech,opt$pars,LEVIER,N,t=length(ech[,"r"]),rv=model[(t+1),'rvt'])
                pi_0 <- l$w_hat
                sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
                
                if(LEVIER){
                  Levier<-levierVolatility(ech[(length(ech[,"r"])-100):(length(ech[,"r"])),"r"],Nl=70,para=para)$`levier`
                  sig<- sig*Levier
                }
                
                model[t+1,"loglik"]   <- -as.numeric(opt$values[length(opt$values)])
                model[t+1,"pred_lik"] <- l$Pred_loglik
                model[t+1,'AIC']      <- model[(t+1),"loglik"]-length(para)
                model[t+1,'BIC']      <- model[(t+1),"loglik"]-length(para)*log(length(ech))/2
                
                model[t+1,]
              }
  close(pb)
  
  model[update_date[!(update_date==forecast.length)]+1,]<-Y
  
  cat("Estimation step 2: \n")
  pb <- txtProgressBar(min=0, max = (forecast.length-1), style = 3)
  for(t in 0:(forecast.length-1)){
    ech    <- donne[strt:(length(nam_)-forecast.length+t),]
    
    if(t %in% update_date){
      para <- unlist(model[t+1,vars])
      next
    }
    
    model[t+1,vars] <- round(para,5)
    para_tilde<-natWork(para,LEVIER)
    l    <- logLik2(ech,para_tilde,LEVIER,N,t=length(ech[,"r"]),rv=model[(t+1),'rvt'])
    pi_0 <- l$w_hat
    sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
    
    if(LEVIER){
      Levier<-levierVolatility(ech[(length(ech[,"r"])-100):(length(ech[,"r"])),"r"],Nl=70,para=para)$`levier`
      sig<- sig*Levier
    }
    
    model[t+1,"loglik"]   <- -l$loglik
    model[t+1,"pred_lik"] <- l$Pred_loglik
    model[t+1,'AIC']      <- model[(t+1),"loglik"]-length(para)
    model[t+1,'BIC']      <- model[(t+1),"loglik"]-length(para)*log(length(ech))/2
    setTxtProgressBar(pb, t)
  }
  close(pb)
  
  registerDoSNOW(cluster)
  cat("Prevision step : \n")
  n <- 1000 #number of simulations
  Nl<-70
  H <- max(Loss.horizon)    #nb of years simulated
  pb <- txtProgressBar(min=1, max = forecast.length, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  model <- foreach(t=1:nrow(model), .noexport=c("natWork","workNat", "logLik", "logLik2","volatilityVector", "P",
                                                "levierVolatility", "R_hat","f_sim"),
                   .packages = c("Rcpp","RcppArmadillo","RcppEigen"), .combine = "rbind",
                   .export=c("sim.mc"), .options.snow=opts) %dopar% {
            
                     sourceCpp(paste0("benchmarks/FHMV_rv.cpp"))
                     
                     ech    <- donne[strt:(length(nam_)-forecast.length+t-1),]
                     
                     para <- unlist(model[t, vars])
                     para_tilde<-natWork(para,LEVIER)
                     l<-logLik2(ech,para_tilde,LEVIER,N,t=length(ech[,"r"]),rv=model[t,'rvt'])
                     
                     pi_0 <- l$w_hat
                     if(t %in% update_date+1){
                       if(refit.window == "moving") strt <- strt + refit.every
                       sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
                       matP   <- P(para["p"],N)
                     }
                     
                     MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
                     z_t    <- matrix(rnorm(n*H), nrow = n , ncol = H)
                     
                     rt2_sim<-f_sim(H,sig,pi_0,matP)
                     
                     if(LEVIER){
                       Levier<-rep(1,n)%*%t(levierVolatility(ech[(length(ech[,2])-200):(length(ech[,2])),2],Nl,para)$`Levier`)
                       sim<-R_hat(H,ech[(length(ech[,2])-200):(length(ech[,2])),2],z_t,Levier,para,N,Nl)
                       rt2_sim <-rt2_sim*sim[(ncol(sim)-H+1):ncol(sim)]
                     }
                     
                     
                     for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[1:Loss.horizon[k]])
                     names(R_var) <- paste0("RV_for",Loss.horizon)
                     model[t,colnames(model) %in% names(R_var)]<-R_var
                     
                     ech    <- donne[(length(nam_)-forecast.length+t):length(nam_),"rv"]
                     for(k in 1:length(Loss.horizon)) R_var[k] <- sum(ech[1:Loss.horizon[k]])
                     names(R_var) <- paste0("RV_tru",Loss.horizon)
                     model[t,colnames(model) %in% names(R_var)]<-R_var
                     
                     model[t,]
                   }
  close(pb)
  
  
  model_extern[i,'pred_dens'] <- sum(model[,"pred_lik"])
  
  for_R_var   <- model[,grep('RV_for', colnames(model), fixed=TRUE)]
  for_R_err   <- model[,grep('RV_tru', colnames(model), fixed=TRUE)]
  model_extern[i,paste('QLIKc_RV',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEc_RV',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEc_RV',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  model_extern[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  
  filenam<-"Forecast_FHMV_RV"
  if(file.exists(paste0(filenam,".csv"))){
    tmp <- read.csv(paste0(filenam,".csv"))
    tmp[i,9:ncol(tmp)] <- model_extern[i,9:ncol(model_extern)]
    write.csv(tmp, paste0(filenam,".csv"), row.names=FALSE)
  }else{
    write.csv(model_extern, paste0(filenam,".csv"), row.names=FALSE)
  }
  
}
stopCluster(cluster)
write.csv(model_extern, "Forecast_FHMV_RV.csv", row.names=FALSE)
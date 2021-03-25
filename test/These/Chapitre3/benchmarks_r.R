path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/Returns"
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3"
setwd(path)
library(MDSV)
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

#-------------------------------------------------------------------------------
# Returns : ESTIMATION
#-------------------------------------------------------------------------------

#### GARCH and GJR-GARCH

model<-expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   LEVIER=c(FALSE,TRUE), dist = c("norm","std"))

vars<-c('loglik', 'AIC', 'BIC', c("omega","alpha1","beta1","gamma1","shape"),"conv","time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

for(i in 1:nrow(model)){
  start_time <- Sys.time()
  index                 <- as.character(model[i,"index"])
  donne                 <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                  <- rownames(donne)
  donne                 <- donne[,"r"] - mean(donne[,"r"])
  model[i,"start.date"] <- as.Date(nam_[1])
  model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
  model[i,"length"]     <- length(nam_)
  dist                  <- as.character(model[i,"dist"])
  LEVIER                <- model[i,"LEVIER"]
  
  if(LEVIER){
    spec = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                       mean.model = list(armaOrder =c(0,0), include.mean = FALSE), 
                      distribution.model = dist)
  }else{
    spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder =c(0,0), include.mean = FALSE), 
                      distribution.model = dist)
  }
  
  out_fit = ugarchfit(spec = spec, data = donne)
  
  if(1-convergence(out_fit)){
    model[i,'conv']  <- "Convergence."
  }else{
    model[i,'conv']  <- "No Convergence. Retrun the best result."
  }
  model[i,'loglik']    <- likelihood(out_fit)
  model[i,'AIC']       <- likelihood(out_fit) - length(coef(out_fit))
  model[i,'BIC']       <- likelihood(out_fit) - length(coef(out_fit))*log(length(nam_))/2
  model[i,c("omega","alpha1","beta1","gamma1","shape")] <- coef(out_fit)[c("omega","alpha1","beta1","gamma1","shape")]
  model[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
}

X<-model

write.csv(X, "Estim_ALLGARCH.csv", row.names=FALSE)

#### MSM and FHMV

#some functions
para_names<-function(model,LEVIER){
  if(model=="MSM") vars.names<-c("m0","sigma","b","gamma")
  if(model=="FHMV") vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q")
  if(model=="FHMV_rv") vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q","shape")
  if(model=="DSARV") vars.names<-c("omega","phi","delta","gamma")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

para_init<-function(model,LEVIER){
  if(model=="MSM") para<-c(0.8,1.26,5,0.05)
  if(model=="FHMV") para<-#c(0.4908,4.01718,0.35581,0.99615,8.42876,0.8,0.83108)
  c(2.432275905,	1.739229056,	0.999109205,	0.992337682,	15.89259305,	0.930636033,	0.928642337)
  if(model=="DSARV") para<-c(0.31,  0.96, -2.71,  1.53)
  if(LEVIER) para<-c(para,1.89,0.89)
  return(para)
}

ctrl <- list(TOL=1e-15, trace=0)

model<-expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   Model = c("MSM","FHMV"), LEVIER=c(FALSE,TRUE))

vars<-c('loglik', 'AIC', 'BIC', paste0("param",1:9),"time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

#setup parallel backend to use many processors
cluster <- makeCluster(5)
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
  donne                 <- donne[,"r"] - mean(donne[,"r"])
  model[i,"start.date"] <- as.Date(nam_[1])
  model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
  model[i,"length"]     <- length(nam_)
  Model                 <- as.character(model[i,"Model"])
  LEVIER                <- model[i,"LEVIER"]
  N                     <- 2
  
  sourceCpp(paste0("benchmarks/",Model,".cpp"))
  
  para<-para_init(Model,LEVIER)
  para_tilde<-natWork(para,LEVIER)
  opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
  vars<-para_names(Model,LEVIER)
  params<-workNat(opt$pars,LEVIER)
  names(params)<-para_names(Model,LEVIER)
  
  model[i,'loglik']    <- -as.numeric(opt$values[length(opt$values)])
  model[i,'AIC']       <- model[i,"loglik"] - length(params)
  model[i,'BIC']       <- model[i,"loglik"]-(length(params))*log(length(nam_))/2
  model[i,paste0("param",1:length(params))] <- params
  model[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  model[i,]
}

close(pb)
stopCluster(cluster)

write.csv(Y, "Estim_ALLHMM_returns.csv", row.names=FALSE)


#-------------------------------------------------------------------------------
# Returns : PREVISION
#-------------------------------------------------------------------------------

g<-function(vector){
  Sortie<-matrix(NA,2,2)
  for(i in c(0,1)) for(j in c(0,1)) Sortie[i+1,j+1]<-sum((vector[-1]==j)*(vector[-length(vector)]==i))
  if(vector[1]==0) Sortie[1,1]<-Sortie[1,1]+1
  if(vector[1]==1) Sortie[1,2]<-Sortie[1,2]+1
  colnames(Sortie)<-c(0,1)
  rownames(Sortie)<-c(0,1)
  return(Sortie)
}

#### GARCH and GJR-GARCH

model_extern<-expand.grid(index=index_set, LEVIER=c(FALSE,TRUE), dist = c("norm","std"), n.ahead=c(100),
                          forecast.length=c(756), refit.every=c(63), refit.window=c("recursive"),
                          calculate.VaR=c(TRUE), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

R_var <- rep(0,length(Loss.horizon))

for(i in 1:nrow(model_extern)){
  start_time <- Sys.time()
  index                        <- as.character(model_extern[i,"index"])
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne                        <- donne[,"r"] - mean(donne[,"r"])
  model_extern[i,"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i,"end.date"]   <- as.character(nam_[length(nam_)])
  dist                         <- as.character(model_extern[i,"dist"])
  LEVIER                       <- model_extern[i,"LEVIER"]
  calculate.VaR                <- as.logical(model_extern[i,"calculate.VaR"])
  
  set.seed(rseed)
  
  if(LEVIER){
    spec = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder =c(0,0), include.mean = FALSE), 
                      distribution.model = dist)
  }else{
    spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder =c(0,0), include.mean = FALSE), 
                      distribution.model = dist)
  }
  
  filename <- paste("Forecast_GARCH_LEVIER_",LEVIER,"_", index, "_", dist, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  model<-expand.grid(date = nam_[((length(nam_)-forecast.length+1):length(nam_))],rt=0,dist=dist,Levier=LEVIER)
  
  model$rt  <- donne[((length(nam_)-forecast.length+1):length(nam_))]
  
  vars<-c(paste0("R_for",Loss.horizon),paste0("R_tru",Loss.horizon),
          paste0("R_for_m",Loss.horizon),paste0("R_tru_m",Loss.horizon),'pred_lik','loglik', 'AIC', 'BIC',
          c("omega","alpha1","beta1","gamma1","shape"), paste0("VaR",100*(1-VaR.alpha)), 
          paste0("I",100*(1-VaR.alpha)))
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- 1
  
  for(t in 0:(forecast.length-1)){
    ech    <- donne[strt:(length(nam_)-forecast.length+t)]
    
    if(t %in% update_date){
      out_fit <- ugarchfit(spec, data=ech)
      spec2<-spec
      setfixed(spec2) <- as.list(coef(out_fit))
      
      if(refit.window == "moving") strt <- strt + refit.every
    }
    
    opt     <- ugarchboot(spec2,data=ech, method = c("Partial", "Full")[1], n.ahead = n.ahead, 
                          n.bootpred = 10000)
    
    tmp <- colMeans((as.data.frame(opt, which = "series"))^2)
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(tmp[1:Loss.horizon[k]])
    names(R_var) <- paste0("R_for",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(tmp[Loss.horizon[k]])
    names(R_var) <- paste0("R_for_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t+1):length(nam_)]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[1:Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var

    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru_m",Loss.horizon)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    model[(t+1),c("omega","alpha1","beta1","gamma1","shape")] <- coef(out_fit)[c("omega","alpha1","beta1",
                                                                                 "gamma1","shape")]
    model[(t+1),'loglik']    <- likelihood(out_fit)
    model[(t+1),'AIC']       <- likelihood(out_fit) - length(coef(out_fit))
    model[(t+1),'BIC']       <- likelihood(out_fit) - length(coef(out_fit))*log(length(nam_))/2
    if(dist == 'norm') {
      model[(t+1),'pred_lik']  <- dnorm(model[(t+1),'rt'],0,as.data.frame(opt, which = "sigma")[1,1],log=TRUE)
    }else if(dist == 'std'){
      sigma<-as.data.frame(opt, which = "sigma")[1,1]
      shape<-coef(out_fit)["shape"]
      x <- model[(t+1),'rt']/(sigma*sqrt((shape-2)/shape))
      model[(t+1),'pred_lik']  <- log(dt(x,shape)/(sigma*sqrt((shape-2)/shape)))
    }
    
    if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
        model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- quantile(as.data.frame(opt, which = "series")[,1],VaR.alpha[iter])
      }
      
  }

  if(calculate.VaR){
    for(iter in 1:length(VaR.alpha)){
      model[,paste0('I',100*(1-VaR.alpha[iter]))] <- (model[,'rt'] < model[,paste0('VaR',100*(1-VaR.alpha[iter]))])
    }
  }
  
  model_extern[i,'pred_dens'] <- sum(model[,"pred_lik"])
  
  for_R_var   <- model[,grep('R_for', colnames(model), fixed=TRUE)]
  for_R_var   <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru', colnames(model), fixed=TRUE)]
  for_R_err   <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
  model_extern[i,paste('QLIKc_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEc_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEc_R',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  for_R_var   <- model[,grep('R_for_m', colnames(model), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru_m', colnames(model), fixed=TRUE)]
  model_extern[i,paste('QLIKm_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEm_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEm_R',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
  p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
  
  for(iter in 1:length(VaR.alpha)){
    viol      <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    alpha_hat <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])/forecast.length
    LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(forecast.length-viol))))
    
    viol      <- g(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    pi        <- (viol[1,2]+viol[2,2])/sum(viol)
    pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
    pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
    LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                           ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
    
    LR.uc_vec <- c(LR.uc_vec,LR.uc)
    LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
    LR.ind_vec <- c(LR.ind_vec,LR.ind)
    p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
    p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
    p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
  }
  
  names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
  names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")

  model_extern[i,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.uc_vec
  model_extern[i,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.ind_vec
  model_extern[i,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.cc_vec
  
  model_extern[i,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.uc_vec
  model_extern[i,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.ind_vec
  model_extern[i,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.cc_vec
  
  model_extern[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
}

write.csv(model_extern, "Forecast_ALLGARCH_756.csv", row.names=FALSE)



#### MSM and FHMV

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

model_extern<-expand.grid(index=c("n225","sp500","ftse"), Model = c('FHMV'), LEVIER=c(FALSE), n.ahead=c(100),
                          forecast.length=c(756), refit.every=c(63), refit.window=c("recursive"),
                          calculate.VaR=c(TRUE), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050
cluster         <- makeCluster(detectCores()[1]-1)

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

R_var <- rep(0,length(Loss.horizon))

for(i in 1:nrow(model_extern)){
  start_time <- Sys.time()
  index                        <- as.character(model_extern[i,"index"])
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne                        <- donne[,"r"] - mean(donne[,"r"])
  model_extern[i,"start.date"] <- as.character(as.Date(nam_[length(nam_)])-forecast.length)
  model_extern[i,"end.date"]   <- as.character(nam_[length(nam_)])
  Model                        <- as.character(model_extern[i,"Model"])
  LEVIER                       <- model_extern[i,"LEVIER"]
  calculate.VaR                <- as.logical(model_extern[i,"calculate.VaR"])
  N                            <- 10
  
  sourceCpp(paste0("benchmarks/",Model,".cpp"))
  set.seed(rseed)
  
  filename <- paste("Forecast_",Model,"_LEVIER_",LEVIER,"_", index, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  model<-expand.grid(date = nam_[((length(nam_)-forecast.length+1):length(nam_))],rt=0,Levier=LEVIER)
  model$rt  <- donne[((length(nam_)-forecast.length+1):length(nam_))]
  
  vars<-c(paste0("R_for",Loss.horizon),paste0("R_tru",Loss.horizon),
          paste0("R_for_m",Loss.horizon),paste0("R_tru_m",Loss.horizon),'pred_lik','loglik', 'AIC', 'BIC',
          paste0("param",1:9), paste0("VaR",100*(1-VaR.alpha)), 
          paste0("I",100*(1-VaR.alpha)))
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  para<-para_init(Model,LEVIER)
  vars<-para_names(Model,LEVIER)
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
                
                ech    <- donne[strt:(length(nam_)-forecast.length+t)]
                
                sourceCpp(paste0("benchmarks/",Model,".cpp"))
                
                if(!is.null(opt)) para_tilde<-opt$pars
                oldw <- getOption("warn")
                options(warn = -1)
                opt  <- try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
                options(warn = oldw)

                para <- workNat(opt$pars,LEVIER=LEVIER)

                model[t+1,paste0("param",1:length(para))] <- round(para,5)
                l    <- logLik2(ech,opt$pars,LEVIER,N,t=length(ech),r=model[(t+1),'rt'])

                pi_0 <- l$w_hat
                if(Model=="FHMV"){
                  sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
                }else{
                  sig <- volatilityVector(para,N)
                }
                if(LEVIER){
                  Levier<-levierVolatility(ech[(length(ech)-100):(length(ech))],Nl=70,para=para)$`levier`
                  sig<- sig*Levier
                }

                if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
                  va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
                  if(class(va) =='try-error') {
                    va <- try(qmixnorm((1-VaR.alpha[iter]), rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
                    if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
                  }else{
                    model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
                  }
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
    ech    <- donne[strt:(length(nam_)-forecast.length+t)]
    
    if(t %in% update_date){
      para <- unlist(model[t+1,paste0("param",1:length(para))])
      next
    }
    
    model[t+1,paste0("param",1:length(para))] <- round(para,5)
    para_tilde<-natWork(para,LEVIER)
    l    <- logLik2(ech,para_tilde,LEVIER,N,t=length(ech),r=model[(t+1),'rt'])
    
    pi_0 <- l$w_hat
    if(Model=="FHMV"){
      sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
    }else{
      sig <- volatilityVector(para,N)
    }
    if(LEVIER){
      Levier<-levierVolatility(ech[(length(ech)-100):(length(ech))],Nl=70,para=para)$`levier`
      sig<- sig*Levier
    }
    
    if(calculate.VaR) for(iter in 1:length(VaR.alpha)){
      va <- try(qmist2n(VaR.alpha[iter], sigma=sig, p=pi_0),silent=T)
      if(class(va) =='try-error') {
        va <- try(qmixnorm((1-VaR.alpha[iter]), rep(0,length(sig)), sqrt(sig), pi_0),silent=T)
        if(!(class(va) =='try-error')) model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
      }else{
        model[t+1,paste0('VaR',100*(1-VaR.alpha[iter]))] <- va
      }
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
    
    sourceCpp(paste0("benchmarks/",Model,".cpp"))

    ech    <- donne[strt:(length(nam_)-forecast.length+t-1)]
    
    para <- unlist(model[t, paste0("param",1:length(para))])
    para_tilde<-natWork(para,LEVIER)
    l<-logLik2(ech,para_tilde,LEVIER,N,t=length(ech),r=model[t,'rt'])
    
    pi_0 <- l$w_hat
    if(t %in% update_date+1){
      if(refit.window == "moving") strt <- strt + refit.every
      if(Model=="FHMV"){
        sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
      }else{
        sig <- volatilityVector(para,N)
      }
      matP   <- P(para,N)
    }
    
    MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
    z_t    <- matrix(rnorm(n*H), nrow = n , ncol = H)
    
    if(LEVIER){
      Levier<-rep(1,n)%*%t(levierVolatility(ech[(length(ech)-100):(length(ech))],Nl,para)$`Levier`)
      sim<-R_hat(H,ech[(length(ech)-100):(length(ech))],MC_sim,z_t,Levier,sig,para,N,Nl)
      rt2_sim<-sim$`rt2`
      rt2_sim <- rt2_sim[,(ncol(rt2_sim)-H+1):ncol(rt2_sim)]
    } else {
      rt2_sim<-f_sim(H,sig,pi_0,matP)
    }
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[1:Loss.horizon[k]])
    names(R_var) <- paste0("R_for",Loss.horizon)
    model[t,colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[Loss.horizon[k]])
    names(R_var) <- paste0("R_for_m",Loss.horizon)
    model[t,colnames(model) %in% names(R_var)]<-R_var
    
    ech    <- donne[(length(nam_)-forecast.length+t):length(nam_)]
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[1:Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru",Loss.horizon)
    model[t,colnames(model) %in% names(R_var)]<-R_var
    
    for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[Loss.horizon[k]])^2)
    names(R_var) <- paste0("R_tru_m",Loss.horizon)
    model[t,colnames(model) %in% names(R_var)]<-R_var
    
    model[t,]
    #print(paste("===",round(100*t/nrow(model),2) , "%" ,"===="))
  }
  close(pb)
  
  if(calculate.VaR){
    for(iter in 1:length(VaR.alpha)){
      model[,paste0('I',100*(1-VaR.alpha[iter]))] <- (model[,'rt'] < model[,paste0('VaR',100*(1-VaR.alpha[iter]))])
    }
  }
  
  model_extern[i,'pred_dens'] <- sum(model[,"pred_lik"])
  
  for_R_var   <- model[,grep('R_for', colnames(model), fixed=TRUE)]
  for_R_var   <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru', colnames(model), fixed=TRUE)]
  for_R_err   <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
  model_extern[i,paste('QLIKc_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEc_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEc_R',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  for_R_var   <- model[,grep('R_for_m', colnames(model), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru_m', colnames(model), fixed=TRUE)]
  model_extern[i,paste('QLIKm_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEm_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEm_R',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
  p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
  
  for(iter in 1:length(VaR.alpha)){
    viol      <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    alpha_hat <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])/forecast.length
    LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(forecast.length-viol))))
    
    viol      <- g(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    pi        <- (viol[1,2]+viol[2,2])/sum(viol)
    pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
    pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
    LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                           ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
    
    LR.uc_vec <- c(LR.uc_vec,LR.uc)
    LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
    LR.ind_vec <- c(LR.ind_vec,LR.ind)
    p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
    p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
    p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
  }
  
  names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
  names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")
  
  model_extern[i,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.uc_vec
  model_extern[i,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.ind_vec
  model_extern[i,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.cc_vec
  
  model_extern[i,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.uc_vec
  model_extern[i,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.ind_vec
  model_extern[i,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.cc_vec
  
  model_extern[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
}

stopCluster(cluster)
write.csv(model_extern, "Forecast_ALLHMM_756_FALSE.csv", row.names=FALSE)












para_init<-function(model,LEVIER){
  if(model=="MSM") para<-c(0.8,1.26,5,0.05)
  if(model=="FHMV") para<-#c(0.4908,4.01718,0.35581,0.99615,8.42876,0.8,0.83108)
      c(2.432275905,	1.739229056,	0.999109205,	0.992337682,	15.89259305,	0.930636033,	0.928642337)
  if(model=="DSARV") para<-c(0.31,  0.96, -2.71,  1.53)
  if(LEVIER) para<-c(para,1.89,0.89)
  return(para)
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


############### correction

path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess/Forecasts/Returns"
setwd(path)
library(MDSV)
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



model_extern<-expand.grid(index=index_set, Model = c('MSM','FHMV'), LEVIER=c(FALSE,TRUE), n.ahead=c(100),
                          forecast.length=c(756), refit.every=c(63), refit.window=c("recursive"),
                          calculate.VaR=c(TRUE), rseed=1050)

Loss.horizon    <- c(1,5,10,25,50,75,100)
n.ahead         <- 100
forecast.length <- 756
refit.every     <- 63
refit.window    <- "recursive"
VaR.alpha       <- c(0.01,0.05,0.10)
rseed           <- 1050
R_var <- rep(0,length(Loss.horizon))

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_R',Loss.horizon), paste('RMSEc_R',Loss.horizon), 
        paste('MAEc_R',Loss.horizon), paste('QLIKm_R',Loss.horizon), paste('RMSEm_R',Loss.horizon),
        paste('MAEm_R',Loss.horizon), paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10))), paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10))), 
        paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10))),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

files.all <- c(Sys.glob("Forecast_FHMV*.csv"),
               Sys.glob("Forecast_MSM*.csv"))

for(filename in files.all){
  model<-read.csv(filename)
  X<-model
  
  Model <- gsub("Forecast_","",filename)
  Model <- as.character(gsub("\\_.*", "", Model))
  temp<-paste0("Forecast_",Model)
  
  LEVIER <- gsub(paste0(temp,"_LEVIER_"),"",filename)
  LEVIER <- as.logical(gsub("\\_.*", "", LEVIER))
  para<-para_init(Model,LEVIER)
  temp<-paste0(temp,"_LEVIER_",LEVIER)
  
  index <- gsub(paste0(temp,"_"),"",filename)
  index <- gsub("\\_.*", "", index)
  temp<-paste0(temp,"_",index,"_length_756")
  
  donne                        <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                         <- rownames(donne)
  donne                        <- donne[,"r"] - mean(donne[,"r"])
  
  i<-which((model_extern$index==index)&(model_extern$LEVIER==LEVIER)&(model_extern$Model==Model))
  model_extern[i,"start.date"] <- gsub(paste0(temp,"_"),"",filename)
  model_extern[i,"start.date"] <- gsub("\\_.*", "", model_extern[i,"start.date"])
  temp <- paste0(temp,"_",model_extern[i,"start.date"])
  
  model_extern[i,"end.date"]   <- gsub(paste0(temp,"_"),"",filename)
  model_extern[i,"end.date"]   <- gsub("\\..*", "", model_extern[i,"end.date"])
  
  calculate.VaR                <- TRUE
  N                            <- 10
  
  sourceCpp(paste0("C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3/benchmarks/",Model,".cpp"))
  set.seed(rseed)
  
  filename <- paste("Forecast_",Model,"_LEVIER_",LEVIER,"_", index, "_length_",forecast.length,
                    "_", model_extern[i,"start.date"], "_", model_extern[i,"end.date"], sep="")
  
  opt<-NULL
  update_date<-seq(0,forecast.length,by=refit.every)
  strt <- 1
  
  R_for <- model[,grep('R_for', colnames(model), fixed=TRUE)]
  
  n <- 1000 #number of simulations
  Nl<-70
  H <- max(Loss.horizon)    #nb of years simulated
  
  x<-model[1,]
  
  ech    <- donne[strt:(length(nam_)-forecast.length)]
  
  para <- unlist(model[1, paste0("param",1:length(para))])
  para_tilde<-natWork(para,LEVIER)
  l<-logLik2(ech,para_tilde,LEVIER,N,t=length(ech),r=model[1,'rt'])
  
  pi_0 <- l$w_hat
  if(Model=="FHMV"){
    sig <- volatilityVector(para,N)%*%c((1/(N-1))*rep(para[7],N-1),1-para[7])
  }else{
    sig <- volatilityVector(para,N)
  }
  matP   <- P(para,N)
  
  MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
  z_t    <- matrix(rnorm(n*H), nrow = n , ncol = H)
  
  if(LEVIER){
    Levier<-rep(1,n)%*%t(levierVolatility(ech[(length(ech)-100):(length(ech))],Nl,para)$`Levier`)
    sim<-R_hat(H,ech[(length(ech)-100):(length(ech))],MC_sim,z_t,Levier,sig,para,N,Nl)
    rt2_sim<-sim$`rt2`
    rt2_sim <- rt2_sim[,(ncol(rt2_sim)-H+1):ncol(rt2_sim)]
  } else {
    rt2_sim<-f_sim(H,sig,pi_0,matP)
  }
  
  for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[1:Loss.horizon[k]])
  names(R_var) <- paste0("R_for",Loss.horizon)
  x[colnames(model) %in% names(R_var)]<-R_var
  
  for(k in 1:length(Loss.horizon)) R_var[k] <- sum(rt2_sim[Loss.horizon[k]])
  names(R_var) <- paste0("R_for_m",Loss.horizon)
  x[colnames(model) %in% names(R_var)]<-R_var
  
  ech    <- donne[(length(nam_)-forecast.length+1):length(nam_)]
  for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[1:Loss.horizon[k]])^2)
  names(R_var) <- paste0("R_tru",Loss.horizon)
  x[colnames(model) %in% names(R_var)]<-R_var
  
  for(k in 1:length(Loss.horizon)) R_var[k] <- sum((ech[Loss.horizon[k]])^2)
  names(R_var) <- paste0("R_tru_m",Loss.horizon)
  x[colnames(model) %in% names(R_var)]<-R_var
  
  R_for <- rbind(x[,grep('R_for', colnames(model), fixed=TRUE)],R_for[-nrow(R_for)])
  
  model[,grep('R_for', colnames(model), fixed=TRUE)] <-R_for
  
  model_extern[i,'pred_dens'] <- sum(model[,"pred_lik"])
  
  for_R_var   <- model[,grep('R_for', colnames(model), fixed=TRUE)]
  for_R_var   <- for_R_var[,-grep('R_for_m', colnames(for_R_var), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru', colnames(model), fixed=TRUE)]
  for_R_err   <- for_R_err[,-grep('R_tru_m', colnames(for_R_err), fixed=TRUE)]
  model_extern[i,paste('QLIKc_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEc_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEc_R',Loss.horizon)] <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  for_R_var   <- model[,grep('R_for_m', colnames(model), fixed=TRUE)]
  for_R_err   <- model[,grep('R_tru_m', colnames(model), fixed=TRUE)]
  model_extern[i,paste('QLIKm_R',Loss.horizon)] <- colMeans(log(for_R_var) + for_R_err/for_R_var, na.rm=TRUE)
  model_extern[i,paste('RMSEm_R',Loss.horizon)] <- sqrt(colMeans( (for_R_var - for_R_err)^2, na.rm=TRUE ))/Loss.horizon
  model_extern[i,paste('MAEm_R',Loss.horizon)]  <- colMeans( abs(for_R_var - for_R_err), na.rm=TRUE )/Loss.horizon
  
  LR.uc_vec <- LR.cc_vec <- LR.ind_vec <- NULL
  p.uc_vec  <- p.cc_vec  <- p.ind_vec  <- NULL
  
  for(iter in 1:length(VaR.alpha)){
    viol      <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    alpha_hat <- sum(model[,paste0('I',100*(1-VaR.alpha[iter]))])/forecast.length
    LR.uc     <- 2*log(((alpha_hat^viol)*((1-alpha_hat)^(forecast.length-viol)))/((VaR.alpha[iter]^viol)*((1-VaR.alpha[iter])^(forecast.length-viol))))
    
    viol      <- g(model[,paste0('I',100*(1-VaR.alpha[iter]))])
    pi        <- (viol[1,2]+viol[2,2])/sum(viol)
    pi0       <- (viol[1,2])/(viol[1,1]+viol[1,2])
    pi1       <- (viol[2,2])/(viol[2,2]+viol[2,1])
    LR.ind    <- - 2*log((((1-pi)^(viol[1,1]+viol[2,1]))*(pi^(viol[1,2]+viol[2,2])))/
                           ((((1-pi0)^viol[1,1])*(pi0^viol[1,2]))*(((1-pi1)^viol[2,1])*(pi1^viol[2,2]))))
    
    LR.uc_vec <- c(LR.uc_vec,LR.uc)
    LR.cc_vec <- c(LR.cc_vec,LR.uc+LR.ind)
    LR.ind_vec <- c(LR.ind_vec,LR.ind)
    p.uc_vec <- c(p.uc_vec,1-pchisq(LR.uc,1))
    p.cc_vec <- c(p.cc_vec,1-pchisq(LR.uc+LR.ind,2))
    p.ind_vec <- c(p.ind_vec,1-pchisq(LR.ind,1))
  }
  
  names(LR.uc_vec) <- names(LR.cc_vec) <- names(LR.ind_vec) <- paste(100*(1-VaR.alpha), "%")
  names(p.uc_vec)  <- names(p.cc_vec)  <- names(p.ind_vec) <- paste(100*(1-VaR.alpha), "%")
  
  model_extern[i,paste('LR.uc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.uc_vec
  model_extern[i,paste('LR.ind_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.ind_vec
  model_extern[i,paste('LR.cc_stat',100*(1-c(0.01,0.05,0.10)))] <- LR.cc_vec
  
  model_extern[i,paste('LR.uc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.uc_vec
  model_extern[i,paste('LR.ind_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.ind_vec
  model_extern[i,paste('LR.cc_pvalue',100*(1-c(0.01,0.05,0.10)))] <- p.cc_vec
  
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  
}
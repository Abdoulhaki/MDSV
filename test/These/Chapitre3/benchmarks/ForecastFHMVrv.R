path<-"/home/maoudek/Rsim/Article1/MDSV/Returns"
# path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/MDSV_package/MDSV/test/These/Chapitre3/benchmarks"
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
cluster         <- makeCluster(3)

vars<-c('start.date','end.date','pred_dens',paste('QLIKc_RV',Loss.horizon), paste('RMSEc_RV',Loss.horizon), 
        paste('MAEc_RV',Loss.horizon),"time")

model_add <- matrix(0, nrow=nrow(model_extern), ncol=length(vars))
colnames(model_add) <- vars
model_extern <- cbind(model_extern, model_add)

R_var <- rep(0,length(Loss.horizon))

# for(i in 1:nrow(model_extern)){
for(i in 3:3){
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
  
  sourceCpp(paste0("FHMV_rv.cpp"))
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
                
                sourceCpp(paste0("FHMV_rv.cpp"))
                
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
  # for(t in 1:nrow(model)){                   
                     sourceCpp(paste0("FHMV_rv.cpp"))
                     
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
  
  filenam<-"Forecast_FHMV_RV2"
  if(file.exists(paste0(filenam,".csv"))){
    tmp <- read.csv(paste0(filenam,".csv"))
    tmp[i,9:ncol(tmp)] <- model_extern[i,9:ncol(model_extern)]
    write.csv(tmp, paste0(filenam,".csv"), row.names=FALSE)
  }else{
    write.csv(model_extern, paste0(filenam,".csv"), row.names=FALSE)
  }
  
}
stopCluster(cluster)

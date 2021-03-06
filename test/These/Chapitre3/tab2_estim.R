path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
setwd(path)
library(MDSV)
if(!require(timeDate)){install.packages("timeDate")}; library(timeDate)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

index_set<-c("sp500", "nasdaq", "ftse", "n225")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

#-------------------------------------------------------------------------------
# Mise en place
#-------------------------------------------------------------------------------

model<-expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   K=c(1:10,15,20,25,30,32,40,50,75,100,250,500,750,1000,1024), N=c(1:10), 
                   ModelType=c(0:2), LEVIER=c(FALSE,TRUE))

model$nstate <- (model$K)^(model$N)
model        <- model[!(model$nstate>1024),]
model        <- model[!(model$nstate==1),]
model        <- model[order(model$nstate),]

vars<-c('loglikm','loglik', 'AIC', 'BIC', c("omega","a","b","sigma","v0","xi",
                                            "varphi","delta1","delta2","shape","l","theta"),"conv","time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)


#-------------------------------------------------------------------------------
# Calcul des maxima de vraisemblance
#-------------------------------------------------------------------------------

for(i in 1060:nrow(model)){
  start_time <- Sys.time()
  index                 <- as.character(model[i,"index"])
  donne                 <- get(index)[rownames(get(index))>=start.date & rownames(get(index))<=end.date,]
  nam_                  <- rownames(donne)
  donne[,"r"]           <- donne[,"r"] - mean(donne[,"r"])
  model[i,"start.date"] <- as.Date(nam_[1])
  model[i,"end.date"]   <- as.Date(nam_[length(nam_)])
  model[i,"length"]     <- length(nam_)
  ModelType             <- as.numeric(model[i,"ModelType"])
  LEVIER                <- model[i,"LEVIER"]
  N                     <- model[i,"N"]
  K                     <- model[i,"K"]
  
  #optimization
  out_fit     <- MDSVfit(N = N, K = K, data = donne, ModelType = ModelType, LEVIER = LEVIER)
  para        <- out_fit$estimates # parameter
  model[i,names(para)] <- round(para,5)
  if(1-out_fit$convergence){
    model[i,'conv']  <- "Convergence."
  }else{
    model[i,'conv']  <- "No Convergence. Retrun the best result."
  }
  model[i,'loglik']    <- out_fit$LogLikelihood
  model[i,'AIC']       <- out_fit$AIC
  model[i,'BIC']       <- out_fit$BIC
  
  if(ModelType==2){
    if(N==1) para       <-c(para[1:2],b=2,para[3:length(para)])
    out_filter          <- MDSVfilter(K = K, N = N, data = donne, para = para, ModelType = ModelType, LEVIER = LEVIER)
    model[i,'loglikm']  <- out_filter$Marg_loglik
  }
  
  if(ModelType==0) model[i,"ModelType"] <- "log-return"
  if(ModelType==1) model[i,"ModelType"] <- "realized variances"
  if(ModelType==2) model[i,"ModelType"] <- "Joint"
  model[i,"time"] <- difftime(Sys.time(), start_time,units = "secs")[[1]]
  
  write.csv(model, paste0("All_estimations_suite.csv"), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model),2) , "%" ,"===="))
}

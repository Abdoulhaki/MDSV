path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/Returns"
setwd(path)
library(MDSV)
if(!require(parallel)){install.packages("parallel")}; library(parallel)
if(!require(doSNOW)){install.packages("doSNOW")}; library(doSNOW)

index_set<-c("aex", "aord", "bell", "bsesn", "bvsp", "cac40", "dax", "dji", "ftmib", "ftse",
             "hsi", "ibex35", "kospi", "kse", "mxx", "n225", "nasdaq", "nifty50", "omxc20", "omxhpi",
             "omxspi", "oseax", "psi", "rut", "smsi", "sp500", "ssec", "ssmi", "sti", "stoxx50", "tsx")

index_set<-c("sp500", "nasdaq", "ftse", "n225")

start.date <- as.Date("2000-01-01") 
end.date   <- as.Date("2019-12-31")

#-------------------------------------------------------------------------------
# Mise en place
#-------------------------------------------------------------------------------

# model<-rbind(expand.grid(index=c("sp500"), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
#                          K=c(2:1024), N=c(1), ModelType=c(2), LEVIER=c(FALSE)),
#              expand.grid(index=c("sp500"), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
#                          K=c(2), N=c(2:10), ModelType=c(2), LEVIER=c(FALSE)))

model<-rbind(expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                   K=c(2), N=c(10), ModelType=c(0:2), LEVIER=c(FALSE,TRUE)),
             expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(4), N=c(5), ModelType=c(0:2), LEVIER=c(FALSE,TRUE)),
             expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(10), N=c(3), ModelType=c(0:2), LEVIER=c(FALSE,TRUE)),
             expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(32), N=c(2), ModelType=c(0:2), LEVIER=c(FALSE,TRUE)),
             expand.grid(index=index_set, start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(1024), N=c(1), ModelType=c(0:2), LEVIER=c(FALSE,TRUE)))

model<-rbind(expand.grid(index=c("ftse"), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(1024), N=c(1), ModelType=c(0), LEVIER=c(FALSE)),
             expand.grid(index=("nikkei"), start.date=as.Date("2000-01-01"), end.date=as.Date("2019-12-31"), length=0,
                         K=c(1024), N=c(1), ModelType=c(0), LEVIER=c(TRUE)))

model$nstate <- (model$K)^(model$N)

vars<-c('loglikm','loglik', 'AIC', 'BIC', c("omega","a","b","sigma","v0","xi",
                                            "varphi","delta1","delta2","shape","l","theta"),"conv","time")

model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

#-------------------------------------------------------------------------------
# Calcul des maxima de vraisemblance
#-------------------------------------------------------------------------------
cluster <- makeCluster(5)
registerDoSNOW(cluster)
pb <- txtProgressBar(min=61, max = nrow(model), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

Y<- foreach(i=61:nrow(model), .export=c("MDSVfit","MDSVfilter"),.packages=c("MDSV"), .combine = "rbind",.options.snow=opts) %dopar% {
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
  
  model[i,]
}

close(pb)
stopCluster(cluster)

write.csv(Y, paste0("Estimations_MDSV_ALL_2.csv"), row.names=FALSE)

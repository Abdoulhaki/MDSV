#-------------------------------------------------------------------------------
#Data input
#-------------------------------------------------------------------------------
library(MDSV)
if(!require(parallel)){install.packages("parallel")}; library(parallel)
if(!require(doSNOW)){install.packages("doSNOW")}; library(doSNOW)

# Loading of functions code
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/Code_These/Ess"
path<-"/home/maoudek/Rsim/Article1/MDSV/MonteCarlo"
setwd(path)

panel<-"D" # "A", "B", "C", "D"

#QUELQUES FONCTIONS 

para_names<-function(LEVIER){
  vars.names<-c("omega","a","b","sigma","v0")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

colSdColMeans <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

colFSSE <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

colRMSE <- function(x,m, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  
  colVar <- (colFSSE(x))^2 + (colMeans(x, na.rm=na.rm) - m)^2
  return(sqrt(colVar))
}


# Setting of parameters
ctrl <- list(TOL=1e-15, trace=0) #optimization parameter

sigma <- 1
a <- 0.95
b <- 3

#setup parallel backend to use many processors
cores=detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
cl <- makeCluster(17) #not to overload your computer

if(panel == "A"){
  #### PANEL A : K = 10, N = 2
  
  model_extern <- expand.grid(K=c(10),N=c(2),sigma=sigma,a=a,b=b,omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000),Time=0)
  
  for(ij in 1:nrow(model_extern)){
    
    Time1<-Sys.time()
    para      <- c(model_extern[ij,"omega"],model_extern[ij,"a"],model_extern[ij,"b"],model_extern[ij,"sigma"],model_extern[ij,"v0"])
    N         <- model_extern[ij,"N"]
    K         <- model_extern[ij,"K"]
    n.sim     <- model_extern[ij,"T"]
    m.sim     <- 1
    n.start   <- 0
    ModelType <- 0
    LEVIER    <- FALSE
    
    set.seed(1005)
    
    if(LEVIER) {
      para<-c(para,model_extern[ij,"l"],model_extern[ij,"theta_l"])
    }
    
    filename <- paste0("MonteCarlo1_",ij)
    
    nordi <- 1
    n     <- 1000
    seed  <- sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = n, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    ordi <- 1 
    model<-data.frame(LEVIER=FALSE,K=K,N=N,N_T=n.sim,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
    
    Y<-foreach(i=1:n, .combine=rbind, .export=c("MDSVsim", "MDSVfit"), .packages = c("MDSV"), .options.snow = opts) %dopar% {
      X   <- as.matrix(MDSVsim(N = N, K = N, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
                           n.start = n.start, m.sim = m.sim, rseed = seed[i+n*(ordi-1)])$sim.1$r_t)
      
      opt <- MDSVfit(N = N, K = N, data = X, ModelType = ModelType, LEVIER = LEVIER, ctrl=ctrl)
      
      model[1,colnames(model) %in% names(opt$estimates)] <- round(opt$estimates,5)
      model[1,"loglik"]<-opt$LogLikelihood
      model[1,'N_T']<-length(X)
      model[1,'Np']<-length(opt$estimates)
      model[1,'AIC'] <- opt$AIC
      model[1,'BIC'] <- opt$BIC
      
      model[1,]
    }
    
    close(pb)
    # stopCluster(cl)
    on.exit(stopCluster(cl))
    
    model_extern$Time[ij]<-Sys.time()-Time1
    
    write.csv(Y, paste0(filename,"_K_",K,"_N_",N,".csv"), row.names=FALSE)
  }
  
  
  ##### #Traitement des bases de donnees : PANEL A
  
  files.all<-Sys.glob("MonteCarlo1_*_K_10_N_2.csv")
  
  S<-expand.grid(omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000))
  S$Name<-apply(S,1,FUN=function(x) paste(c("w","v0","T"),x,collapse = '_'))
  X<-matrix(0,16,18)
  X<-as.data.frame(X)
  rownames(X)<-c("length","v0","FSSEv0","RMSEv0","sigma","FSSEsig","RMSEsig","a","FSSEa","RMSEa","b","FSSEb","RMSEb","w","FSSEw","RMSEw")
  names(X)<-S$Name
  filename<-files.all[1]
  K<-10
  N<-2
  
  for(filename in files.all){
    out<-read.csv(filename)
    
    index<-(out[,"a"]>0)&(out[,"a"]<1)&(out[,"b"]>1)&(out[,"v0"]>0)&(out[,"v0"]<1)&(out[,"omega"]>0)&(out[,"omega"]<1)&(out[,"omega"]>0)&(out[,"b"]<10)
    out<-out[index,]
    k<-as.numeric(sub("\\_.*", "",sub(".*MonteCarlo1_", "", filename)))
    X["length",which(filename == files.all)]<-nrow(out)
    X[c("v0","sigma","a","b","w"),S$Name[k]]<-round(colMeans(out[,c("v0","sigma","a","b","omega")]),5)
    X[paste0("FSSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colFSSE(out[,c("v0","sigma","a","b","omega")]),5)
    X[paste0("RMSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colRMSE(out[,c("v0","sigma","a","b","omega")],c(S$v0[k],1.00,0.95,3,S$omega[k])),5)
  }
  
  index<-kronecker(kronecker(c(1:3),c(0,3),"+"),c(0,6,12),"+")
  
  View(X[,index])
  
  write.csv(X[,index], paste0("MONTECARLO_K_",K,"_N_",N,".csv"))
}else if(panel == "B"){
  #### PANEL B
  
  model_extern <- expand.grid(K=c(2),N=c(8),sigma=sigma,a=a,b=b,omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000),Time=0)
  
  for(ij in 1:(nrow(model_extern))){
    
    Time1<-Sys.time()
    para      <- c(model_extern[ij,"omega"],model_extern[ij,"a"],model_extern[ij,"b"],model_extern[ij,"sigma"],model_extern[ij,"v0"])
    N         <- model_extern[ij,"N"]
    K         <- model_extern[ij,"K"]
    n.sim     <- model_extern[ij,"T"]
    m.sim     <- 1
    n.start   <- 0
    rseed     <- 1005
    ModelType <- 0
    LEVIER    <- FALSE
    
    set.seed(1005)
    
    if(LEVIER) {
      para<-c(para,model_extern[ij,"l"],model_extern[ij,"theta_l"])
    }
    
    filename <- paste0("MonteCarlo1_",ij)
    
    nordi <- 1
    n     <- 1000
    seed  <- sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = n, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    ordi <- 1
    model<-data.frame(LEVIER=FALSE,K=K,N=N,N_T=n.sim,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
    
    Y<-foreach(i=1:n, .combine=rbind, .export=c("MDSVsim", "MDSVfit"), .packages = c("MDSV"), .options.snow = opts) %dopar% {
      
      X<-as.matrix(MDSVsim(N = N, K = K, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
                           n.start = n.start, m.sim = m.sim, rseed = seed[i+n*(ordi-1)])$sim.1$r_t)
      
      opt<-MDSVfit(N = N, K = K, data = X, ModelType = ModelType, LEVIER = LEVIER, ctrl=ctrl)
      
      model[1,colnames(model) %in% names(opt$estimates)] <- round(opt$estimates,5)
      model[1,"loglik"]<-opt$LogLikelihood
      model[1,'N_T']<-length(X)
      model[1,'Np']<-length(opt$estimates)
      model[1,'AIC'] <- opt$AIC
      model[1,'BIC'] <- opt$BIC
      
      model[1,]
    }
    
    close(pb)
    # stopCluster(cl)
    on.exit(stopCluster(cl))
    
    model_extern$Time[ij]<-Sys.time()-Time1
    
    write.csv(Y, paste0(filename,"_K_",K,"_N_",N,".csv"), row.names=FALSE)
  }
  
  ##### #Traitement des bases de donnees : PANEL B

  files.all<-Sys.glob("MonteCarlo1_*_K_2_N_8.csv")

  S<-expand.grid(omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000))
  S$Name<-apply(S,1,FUN=function(x) paste(c("w","v0","T"),x,collapse = '_'))
  X<-matrix(0,16,18)
  X<-as.data.frame(X)
  rownames(X)<-c("length","v0","FSSEv0","RMSEv0","sigma","FSSEsig","RMSEsig","a","FSSEa","RMSEa","b","FSSEb","RMSEb","w","FSSEw","RMSEw")
  names(X)<-S$Name
  filename<-files.all[1]
  K<-2
  N<-8

  for(filename in files.all){
    out<-read.csv(filename)

    index<-(out[,"a"]>0)&(out[,"a"]<1)&(out[,"b"]>1)&(out[,"v0"]>0)&(out[,"v0"]<1)&(out[,"omega"]>0)&(out[,"omega"]<1)&(out[,"omega"]>0)&(out[,"b"]<10)
    out<-out[index,]
    X["length",which(filename == files.all)]<-nrow(out)
    k<-as.numeric(sub("\\_.*", "",sub(".*MonteCarlo1_", "", filename)))
    X[c("v0","sigma","a","b","w"),S$Name[k]]<-round(colMeans(out[,c("v0","sigma","a","b","omega")]),5)
    X[paste0("FSSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colFSSE(out[,c("v0","sigma","a","b","omega")]),5)
    X[paste0("RMSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colRMSE(out[,c("v0","sigma","a","b","omega")],c(S$v0[k],1.00,0.95,3,S$omega[k])),5)
  }

  index<-kronecker(kronecker(c(1:3),c(0,3),"+"),c(0,6,12),"+")

  View(X[,index])

  write.csv(X[,index], paste0("MONTECARLO_K_",K,"_N_",N,".csv"))
}else if(panel == "C"){
  ###### PANEL C : MONTE CARLO SUR N #####
  
  # Setting of parameters
  sigma     <- 0.8
  a         <- 0.95
  b         <- 3
  omega     <- 0.8
  v0        <- 0.5
  para      <- c(omega, a, b, sigma, v0)
  n.sim     <- 10000
  m.sim     <- 1
  n.start   <- 0
  ModelType <- 0
  K<-2
  N<-5
  N_max<-8
  K_max<-8
  LEVIER<-FALSE
  
  # Data simulation
  set.seed(1005)
  
  nordi<-2
  n <-500
  seed<-sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ordi<-2
  model<-expand.grid(ech=1,LEVIER=FALSE,K=K,N=c(1:N_max),N_T=n.sim,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
  
  Y<-foreach(i=1:n, .combine=rbind, .export=c("MDSVsim", "MDSVfit"), .packages = c("MDSV"), .options.snow = opts) %dopar% {
    X<-as.matrix(MDSVsim(N = N, K = K, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
                         n.start = n.start, m.sim = m.sim, rseed = seed[i+n*(ordi-1)])$sim.1$r_t)
    for(N in 1:N_max){
      opt<-MDSVfit(N = N, K = K, data = X, ModelType = ModelType, LEVIER = LEVIER, ctrl=ctrl)
      
      model[N,"ech"] <- i+n*(ordi-1)
      model[N,colnames(model) %in% names(opt$estimates)] <- round(opt$estimates,5)
      model[N,"loglik"]<-opt$LogLikelihood
      model[N,'N_T']<-length(X)
      model[N,'Np']<-length(opt$estimates)
      model[N,'AIC'] <- opt$AIC
      model[N,'BIC'] <- opt$BIC
    }
    model[1:N_max,]
  }
  
  close(pb)
  # stopCluster(cl)
  on.exit(stopCluster(cl))
  
  write.csv(Y, paste0("____MonteCarlo2_",ordi,"_K_",2,"_N_",5,".csv"), row.names=FALSE)
  
  ##### #Traitement des bases de donnees : PANEL C

  files.all<-Sys.glob("____MonteCarlo2_*.csv")

  files.all<-files.all[grepl("_N_5", files.all)]
  Y<-NULL
  for(filename in files.all){
    Y<-rbind(Y,read.csv(filename))
  }
  
  Z<-Y[1:8,"loglik"]
  lig<-8
  for(i in 2:1000){
    Z<-rbind(Z,Y[(lig+1):(lig+8),"loglik"])
    lig<-lig+8
  }
  rownames(Z)<-1:1000
  colnames(Z)<-paste0("N=",1:8)

  U2<-round(colMeans(Z-Z[,"N=5"]),3)
  U3<-round(colSdColMeans(Z-Z[,"N=5"])/sqrt(1000),3)

  Z<-cbind(Z,apply(Z,1,which.max))
  U1<-table(Z[,ncol(Z)])

  U1
  U2
  U3
  
}else if(panel == "D"){
  ###### PANEL D : MONTE CARLO SUR K #####
  
  # Setting of parameters
  sigma     <- 0.8
  a         <- 0.95
  b         <- 3
  omega     <- 0.8
  v0        <- 0.5
  para      <- c(omega, a, b, sigma, v0)
  n.sim     <- 10000
  m.sim     <- 1
  n.start   <- 0
  ModelType <- 0
  K<-2
  N<-5
  N_max<-8
  K_max<-8
  LEVIER<-FALSE
  
  # Data simulation
  set.seed(1005)
  
  nordi<-2
  n <-500
  seed<-sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ordi<-2
  model<-expand.grid(ech=1,LEVIER=FALSE,K=c(1:K_max),N=N,N_T=n.sim,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
  
  Y<-foreach(i=1:n, .combine=rbind, .export=c("MDSVsim", "MDSVfit"), .packages = c("MDSV"), .options.snow = opts) %dopar% {
    X<-as.matrix(MDSVsim(N = N, K = K, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
                         n.start = n.start, m.sim = m.sim, rseed = seed[i+n*(ordi-1)])$sim.1$r_t)
    for(K in 1:K_max){
      opt<-MDSVfit(N = N, K = K, data = X, ModelType = ModelType, LEVIER = LEVIER, ctrl=ctrl)
      
      model[K,"ech"] <- i+n*(ordi-1)
      model[K,colnames(model) %in% names(opt$estimates)] <- round(opt$estimates,5)
      model[K,"loglik"]<-opt$LogLikelihood
      model[K,'N_T']<-length(X)
      model[K,'Np']<-length(opt$estimates)
      model[K,'AIC'] <- opt$AIC
      model[K,'BIC'] <- opt$BIC
    }
    model[1:K_max,]
  }
  
  close(pb)
  # stopCluster(cl)
  on.exit(stopCluster(cl))
  
  write.csv(Y, paste0("MonteCarlo2_",ordi,"_K_",5,"_N_",2,".csv"), row.names=FALSE)
  
  
  # ##### #Traitement des bases de donnees : PANEL D
  # 
  # ## K
  # 
  # files.all<-Sys.glob("MonteCarlo1_*.csv")
  # 
  # files.all<-files.all[grepl("_K_5", files.all)]
  # 
  # Y<-read.csv(filename)
  # 
  # Z<-Y[1:8,"loglik"]
  # lig<-8
  # for(i in 2:1000){
  #   Z<-rbind(Z,Y[(lig+1):(lig+8),"loglik"])
  #   lig<-lig+8
  # }
  # rownames(Z)<-1:1000
  # colnames(Z)<-paste0("K=",2:9)
  # 
  # U2<-round(colMeans(Z-Z[,"K=5"]),3)
  # U3<-round(colSdColMeans(Z-Z[,"K=5"])/sqrt(1000),3)
  # 
  # Z<-cbind(Z,apply(Z,1,which.max)+1)
  # U1<-table(Z[,ncol(Z)])
  # 
  # U1
  # U2
  # U3
  
}

#-------------------------------------------------------------------------------
#Data input
#-------------------------------------------------------------------------------
if(!require(MDSV)){install.packages("MDSV")}; library(MDSV)
if(!require(mhsmm)){install.packages("mhsmm")}; library(mhsmm)
if(!require(Rsolnp)){install.packages("Rsolnp")}; library(Rsolnp)
if(!require(parallel)){install.packages("parallel")}; library(parallel)
if(!require(doSNOW)){install.packages("doSNOW")}; library(doSNOW)

# Loading of functions code
path<-"C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Code/MDSV/MonteCarlo"
setwd(path)


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
cl <- makeCluster(cores[1]-1) #not to overload your computer

#### PANEL A

model_extern <- expand.grid(K=c(10),N=c(2),sigma=sigma,a=a,b=b,omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000),Time=0)

for(ij in 1:nrow(model_extern)){

  Time1<-Sys.time()
  para <- c(model_extern[ij,"omega"],model_extern[ij,"a"],model_extern[ij,"b"],model_extern[ij,"sigma"],model_extern[ij,"v0"])
  N <- model_extern[ij,"N"]
  K <- model_extern[ij,"K"]
  v<-volatilityVector(para,K,N)
  MatrixP<-P(para,K,N)
  stationnaryDist<-probapi(model_extern[ij,"omega"],K,N)
  T1<-model_extern[ij,"T"]
  
  LEVIER<-FALSE
  if(LEVIER) {
    para<-c(para,model_extern[ij,"l"],model_extern[ij,"theta_l"])
  }
  para_tilde<-natWork(para+0.01*sample(c(-1,1),length(para),TRUE,c(0.5,0.5)),LEVIER)
  
  # Data simulation
  
  set.seed(1005)
  
  filename <- paste0("MonteCarlo1_",ij)
  
  n <-1000
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  model<-data.frame(N=1,LEVIER=FALSE,K=K,N=N,N_T=T1,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
  
  Y<-foreach(i=1:n, .combine=rbind, .export=c("sim.mc", "solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
             .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
               
       sourceCpp("MDSV_r.cpp")
       X<-rnorm(T1,0,sqrt(v[sim.mc(stationnaryDist,MatrixP,T1)]))
       
       opt<-try(solnp(pars=para_tilde,fun=logLik,ech=X,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
       
       if(!is(opt,"try-error")){
         vars<-c("omega","a","b","sigma","v0")
         if(LEVIER) {
           vars<-c(vars,"l","theta_l")
         }
         
         params<-workNat(opt$pars,LEVIER)
         names(params)<-para_names(LEVIER)
         
         model[1,colnames(model) %in% names(params)] <- round(params[vars],5)
         model[1,"loglik"]<--as.numeric(opt$values[length(opt$values)])
         model[1,'N_T']<-length(X)
         model[1,'Np']<-length(params)
         model[1,'AIC'] <- model[1,"loglik"]-model[1,'Np']
         model[1,'BIC'] <- model[1,"loglik"]-(model[1,'Np'])*log(model[1,'N_T'])/2
         model[1,'BIC'] <- model[1,"loglik"]-(model[1,'Np'])*log(model[1,'N_T'])/2
       }
       model[1,]

  }
  
  close(pb)
  # stopCluster(cl)
  on.exit(stopCluster(cl))
  
  model_extern$Time[ij]<-Sys.time()-Time1
  
  write.csv(Y, paste0(filename,"_N_D_",K,"_N_",N,".csv"), row.names=FALSE)
}

#### PANEL B

model_extern <- expand.grid(K=c(2),N=c(8),sigma=sigma,a=a,b=b,omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000),Time=0)

for(ij in 1:nrow(model_extern)){
  
  Time1<-Sys.time()
  para <- c(model_extern[ij,"omega"],model_extern[ij,"a"],model_extern[ij,"b"],model_extern[ij,"sigma"],model_extern[ij,"v0"])
  N <- model_extern[ij,"N"]
  K <- model_extern[ij,"K"]
  v<-volatilityVector(para,K,N)
  MatrixP<-P(para,K,N)
  stationnaryDist<-probapi(model_extern[ij,"omega"],K,N)
  T1<-model_extern[ij,"T"]
  
  LEVIER<-FALSE
  if(LEVIER) {
    para<-c(para,model_extern[ij,"l"],model_extern[ij,"theta_l"])
  }
  para_tilde<-natWork(para+0.01*sample(c(-1,1),length(para),TRUE,c(0.5,0.5)),LEVIER)
  
  # Data simulation
  
  set.seed(1005)
  
  filename <- paste0("MonteCarlo1_",ij)
  
  n <-1000
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  model<-data.frame(N=1,LEVIER=FALSE,K=K,N=N,N_T=T1,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)
  
  Y<-foreach(i=1:n, .combine=rbind, .export=c("sim.mc", "solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
             .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
               
               sourceCpp("MDSV_r.cpp")
               X<-rnorm(T1,0,sqrt(v[sim.mc(stationnaryDist,MatrixP,T1)]))
               
               opt<-try(solnp(pars=para_tilde,fun=logLik,ech=X,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
               
               if(!is(opt,"try-error")){
                 vars<-c("omega","a","b","sigma","v0")
                 if(LEVIER) {
                   vars<-c(vars,"l","theta_l")
                 }
                 
                 params<-workNat(opt$pars,LEVIER)
                 names(params)<-para_names(LEVIER)
                 
                 model[1,colnames(model) %in% names(params)] <- round(params[vars],5)
                 model[1,"loglik"]<--as.numeric(opt$values[length(opt$values)])
                 model[1,'N_T']<-length(X)
                 model[1,'Np']<-length(params)
                 model[1,'AIC'] <- model[1,"loglik"]-model[1,'Np']
                 model[1,'BIC'] <- model[1,"loglik"]-(model[1,'Np'])*log(model[1,'N_T'])/2
                 model[1,'BIC'] <- model[1,"loglik"]-(model[1,'Np'])*log(model[1,'N_T'])/2
               }
               model[1,]
               
             }
  
  close(pb)
  # stopCluster(cl)
  on.exit(stopCluster(cl))
  
  model_extern$Time[ij]<-Sys.time()-Time1
  
  write.csv(Y, paste0(filename,"_N_D_",K,"_N_",N,".csv"), row.names=FALSE)
}


##### #Traitement des bases de donnees : PANEL A ET PANEL B

files.all<-Sys.glob("MonteCarlo1_*_2_N_8.csv")

S<-expand.grid(omega=c(0.2,0.5,0.8),v0=c(0.5,0.8),T=c(2500,5000,10000))
S$Name<-apply(S,1,FUN=function(x) paste(c("w","v0","T"),x,collapse = '_'))
X<-matrix(0,15,18)
X<-as.data.frame(X)
rownames(X)<-c("v0","FSSEv0","RMSEv0","sigma","FSSEsig","RMSEsig","a","FSSEa","RMSEa","b","FSSEb","RMSEb","w","FSSEw","RMSEw")
names(X)<-S$Name
filename<-files.all[1]

K<-as.numeric(sub("\\_.*", "",sub(".*MonteCarlo1_10_N_D_", "", filename)))
N<-as.numeric(sub("\\..*", "",sub(paste0(".*MonteCarlo1_10_N_D_",K,"_N_"), "", filename)))

for(filename in files.all){
  out<-read.csv(filename)
  
  index<-(out[,"a"]>0)&(out[,"a"]<1)&(out[,"b"]>1)&(out[,"v0"]>0)&(out[,"v0"]<1)&(out[,"omega"]>0)&(out[,"omega"]<1)&(out[,"omega"]>0)&(out[,"b"]<10)
  out<-out[index,]
  k<-as.numeric(sub("\\_.*", "",sub(".*MonteCarlo1_", "", filename)))
  X[c("v0","sigma","a","b","w"),S$Name[k]]<-round(colMeans(out[,c("v0","sigma","a","b","omega")]),5)
  X[paste0("FSSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colFSSE(out[,c("v0","sigma","a","b","omega")]),5)
  X[paste0("RMSE",c("v0","sig","a","b","w")),S$Name[k]]<-round(colRMSE(out[,c("v0","sigma","a","b","omega")],c(S$v0[k],1.00,0.95,3,S$omega[k])),5)
}

index<-kronecker(kronecker(c(1:3),c(0,3),"+"),c(0,6,12),"+")

View(X[,index])

abs(X[c(1,4,7,10,13),16]-c(0.8,1,0.95,3,0.2))
abs(X[c(1,4,7,10,13),10]-c(0.8,1,0.95,3,0.2))

write.table(X[,index], paste0("MONTECARLO_N_D_",K,"_N_",N,".csv"), sep = ";")


# Setting of parameters
sigma <- 1
a <- 0.99
b <- 3
omega <- 0.5
v0 <- 0.8
T1 <- 10000
K<-2
N<-5
N_max<-8
LEVIER<-FALSE

para <- c(omega,a,b,sigma,v0)
v<-volatilityVector(para,K,N)
MatrixP<-P(para,K,N)
stationnaryDist<-probapi(omega,K,N)

para_tilde<-natWork(para+0.01*sample(c(-1,1),length(para),TRUE,c(0.5,0.5)),FALSE)

para_tilde[2]<--6.9067548

###### PANEL C : MONTE CARLO SUR N #####

# Data simulation
set.seed(1005)

nordi<-1
n <-1000
seed<-sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
registerDoSNOW(cl)
pb <- txtProgressBar(max = n, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

ordi<-1
model<-expand.grid(ech=1,LEVIER=FALSE,K=K,N=c(1:N_max),N_T=T1,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)

Y<-foreach(i=1:n, .combine=rbind, .export=c("sim.mc", "solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
           .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
             sourceCpp("MDSV_r.cpp")
             set.seed(seed[i+n*(ordi-1)])
             X<-rnorm(T1,0,sqrt(v[sim.mc(stationnaryDist,MatrixP,T1)]))
             
             for(N in 1:N_max){
               opt<-try(solnp(pars=para_tilde,fun=logLik,ech=X,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
               
               if(!is(opt,"try-error")){
                 vars<-c("omega","a","b","sigma","v0")
                 
                 params<-workNat(opt$pars,LEVIER)
                 names(params)<-para_names(LEVIER)
                 
                 model[N,"ech"] <- i+n*(ordi-1)
                 model[N,colnames(model) %in% names(params)] <- round(params[vars],5)
                 model[N,"loglik"]<--as.numeric(opt$values[length(opt$values)])
                 model[N,'N_T']<-length(X)
                 model[N,'Np']<-length(params)
                 model[N,'AIC'] <- model[N,"loglik"]-model[N,'Np']
                 model[N,'BIC'] <- model[N,"loglik"]-(model[N,'Np'])*log(model[N,'N_T'])/2
                 model[N,'BIC'] <- model[N,"loglik"]-(model[N,'Np'])*log(model[N,'N_T'])/2
                 
               }
             }
             
             model[1:N_max,]
             
           }

close(pb)
# stopCluster(cl)
on.exit(stopCluster(cl))

write.csv(Y, paste0("MonteCarlo2_",ordi,"_N_D_",2,"_N_",5,".csv"), row.names=FALSE)


###### PANEL D : MONTE CARLO SUR K #####

# Data simulation
set.seed(1005)

nordi<-1
n <-1000
seed<-sample.int(min(100*(n*nordi),.Machine$integer.max),(n*nordi))
registerDoSNOW(cl)
pb <- txtProgressBar(max = n, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

ordi<-1
model<-expand.grid(ech=1,LEVIER=FALSE,K=c(1:N_Dmax),N=N,N_T=T1,Np=0,loglik=0,AIC=0,BIC=0,omega=0,a=0,b=0,sigma=0,v0=0,l=0,theta_l=0)

Y<-foreach(i=1:n, .combine=rbind, .export=c("sim.mc", "solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
           .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
             sourceCpp("MDSV_r.cpp")
             set.seed(seed[i+n*(ordi-1)])
             X<-rnorm(T1,0,sqrt(v[sim.mc(stationnaryDist,MatrixP,T1)]))
             
             for(K in 1:N_Dmax){
               opt<-try(solnp(pars=para_tilde,fun=logLik,ech=X,K=K,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
               
               if(!is(opt,"try-error")){
                 vars<-c("omega","a","b","sigma","v0")
                 
                 params<-workNat(opt$pars,LEVIER)
                 names(params)<-para_names(LEVIER)
                 
                 model[K,"ech"] <- i+n*(ordi-1)
                 model[K,colnames(model) %in% names(params)] <- round(params[vars],5)
                 model[K,"loglik"]<--as.numeric(opt$values[length(opt$values)])
                 model[K,'N_T']<-length(X)
                 model[K,'Np']<-length(params)
                 model[K,'AIC'] <- model[K,"loglik"]-model[K,'Np']
                 model[K,'BIC'] <- model[K,"loglik"]-(model[K,'Np'])*log(model[K,'N_T'])/2
                 model[K,'BIC'] <- model[K,"loglik"]-(model[K,'Np'])*log(model[K,'N_T'])/2
                 
               }
             }
             
             model[1:N_Dmax,]
             
           }

close(pb)
# stopCluster(cl)
on.exit(stopCluster(cl))

write.csv(Y, paste0("MonteCarlo2_",ordi,"_N_D_",5,"_N_",2,".csv"), row.names=FALSE)


##### #Traitement des bases de donnees : PANEL C ET PANEL D

## N

files.all<-Sys.glob("_MonteCarlo2_*.csv")

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


## K

files.all<-Sys.glob("_MonteCarlo1_*.csv")

Y<-read.csv(filename)

Z<-Y[1:8,"loglik"]
lig<-8
for(i in 2:1000){
  Z<-rbind(Z,Y[(lig+1):(lig+8),"loglik"])
  lig<-lig+8
}
rownames(Z)<-1:1000
colnames(Z)<-paste0("K=",2:9)

U2<-round(colMeans(Z-Z[,"K=5"]),3)
U3<-round(colSdColMeans(Z-Z[,"K=5"])/sqrt(1000),3)

Z<-cbind(Z,apply(Z,1,which.max)+1)
U1<-table(Z[,ncol(Z)])

U1
U2
U3


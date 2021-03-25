#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::VectorXd workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,
                        const bool &LEVIER){
  Eigen::VectorXd para=para_tilde;
  
  para[0]=exp(para_tilde[0]);        //sigma
  para[1]=1+exp(para_tilde[1]);      //c1
  para[2]=1/(1+exp(para_tilde[2]));  //thetac
  para[3]=1/(1+exp(para_tilde[3]));  //p
  para[4]=1+exp(para_tilde[4]);      //m1
  para[5]=1/(1+exp(para_tilde[5]));  //thetam
  para[6]=1/(1+exp(para_tilde[6]));  //q
  para(7)=exp(para_tilde(7));//shape
  if(LEVIER){ 
    para[8]=exp(para_tilde[8]);      //l
    para[9]=1/(1+exp(para_tilde[9]));//thetal
  }
  return para;
}

// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para,
                        const bool &LEVIER){
  Eigen::VectorXd para_tilde=para;
  
  para_tilde[0]=log(para[0]);        //sigma
  para_tilde[1]=log(para[1]-1);      //c1
  para_tilde[2]=log((1/para[2])-1);  //theta_c
  para_tilde[3]=log((1/para[3])-1);  //p
  para_tilde[4]=log(para[4]-1);      //m1
  para_tilde[5]=log((1/para[5])-1);  //theta_m
  para_tilde[6]=log((1/para[6])-1);  //q
  para_tilde[7]=log(para[7]);//shape
  if(LEVIER){ 
    para_tilde[8]=log(para[8]);      //l
    para_tilde[9]=log((1/para[9])-1);//thetal
  }
  return para_tilde;
}

Eigen::VectorXd volatilityVector_C(const double &c1, const double &thetac, const int &N){
  Eigen::VectorXd C(N+1);
  C[1]=c1;
  C[0]=1/(1+(c1-1)/2);
  for(int i(2); i<(N+1);i++) {
    C[i]=1+pow(thetac,i-1)*(c1-1);
    C[0]=C[0]/(1+pow(thetac,i-1)*(c1-1)/2);
  }
  Eigen::VectorXd x=as<Eigen::VectorXd>(NumericVector::create(c1,1));
  for(int i=2; i<(N+1); i++) x=kroneckerProduct(x,as<Eigen::VectorXd>(NumericVector::create(C[i],1))).eval();
  x=x*C[0];
  return x;
}

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector_M(const double &m1, const double &thetam, const double &q, const int &N){
  NumericVector M(N);
  double M0;
  M[0]=m1;
  M0=1/(1+q*((m1-1)*(1-pow(thetam,N-1)))/((N-1)*(1-thetam)));
  M[N-1]=1;
  for(int i(1); i<(N-1);i++) {M[i]=1+pow(thetam,i)*(m1-1);}
  M=M0*M;
  Eigen::VectorXd x=as<Eigen::VectorXd>(M);
  return x;
}

// [[Rcpp::export]]
Eigen::MatrixXd volatilityVector(const Eigen::VectorXd &para, const int &N){
  Eigen::VectorXd volatility_C=volatilityVector_C(para[1], para[2], N);
  volatility_C=para[0]*volatility_C;
  Eigen::VectorXd volatility_M=volatilityVector_M(para[4],para[5],para[6],N);
  Eigen::MatrixXd volatility=kroneckerProduct(volatility_C,volatility_M.transpose()).eval();
  return volatility;
}

// [[Rcpp::export]]
Eigen::MatrixXd transitionMatrix(const double &p){
  return (2*p-1)*(Eigen::MatrixXd::Identity(2,2))+(1-p)*(Eigen::VectorXd::Ones(2))*(Eigen::VectorXd::Ones(2)).transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd P(const double &p,const int &N){
  Eigen::MatrixXd P=transitionMatrix(p);
  if(N==1) return(P);
  for(int k=1; k<N; k++) P=kroneckerProduct(P,transitionMatrix(p)).eval();
  return(P);
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[8]*pow(para[9],i);
  for(int t=Nl;t<ech.size();t++){
    levier=1;
    for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
    Levier[t]=levier;
  }
  levier=1;
  for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("Levier") = Levier,
      Rcpp::Named("levier") = levier
    );
  return output;
} 

// [[Rcpp::export]]
Eigen::MatrixXd rv_dens(const Eigen::MatrixXd &x,const double &shape){
  int n = x.rows();
  int m = x.cols();
  Eigen::MatrixXd res(n,m);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      res(i,j) = R::dgamma(x(i,j),shape,1.0/shape,false);
    }
  }
  return(res);
}


// [[Rcpp::export]]
double logLik(const Eigen::MatrixXd &ech,
                    const Eigen::Map<Eigen::VectorXd> &para_tilde,
                    const bool &LEVIER,
                    const int &N,
                    const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  // para(0)=0.5;
  int n=ech.rows();
  double shape=para[7];
  Eigen::VectorXd aj;
  Eigen::MatrixXd sigma = volatilityVector(para,N);
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  Eigen::VectorXd p0M=(para[6]/(N-1))*(Eigen::VectorXd::Ones(N));
  p0M[N-1]=1-para[6];
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = P(para[3],N);
  double y=ech(0,0);
  Eigen::MatrixXd bj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
  aj=p0.array()*(bj*p0M).array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  if(!LEVIER){
    for(int i(1); i<n;i++){
      y=ech(i,0);
      bj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
      aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::MatrixXd sigma1;
    for(int i(1); i<n;i++){
      y=ech(i,0);
      sigma1=sigma*Levier[i];
      bj=(rv_dens(y*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
      aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const Eigen::MatrixXd &ech,
                   const Eigen::Map<Eigen::VectorXd> &para_tilde,
                   const bool &LEVIER,
                   const int &N,
                   const double &rv=0,
                   const int &t=2,
                   const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  // para(0)=0.5;
  int n=ech.rows();

  double shape=para[7];
  Eigen::VectorXd aj;
  Eigen::MatrixXd sigma = volatilityVector(para,N);
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  Eigen::VectorXd p0M=(para[6]/(N-1))*(Eigen::VectorXd::Ones(N));
  p0M[N-1]=1-para[6];
  double a(0);
  double lik(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  double y=ech(0,0);
  Eigen::MatrixXd matP = P(para[3],N);
  Eigen::MatrixXd bj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
  aj=p0.array()*(bj*p0M).array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  if(!LEVIER){
    for(int i(1); i<n;i++){
      y=ech(i,0);
      bj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
      aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;

      if(i==(t-1)) w_hat=w;
    }
    bj=(rv_dens(rv*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
    aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::MatrixXd sigma1;
    for(int i(1); i<n;i++){
      y=ech(i,0);
      sigma1=sigma*Levier[i];
      bj=(rv_dens(y*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
      aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;

      if(i==(t-1)) w_hat=w;
    }
    double levier=L["levier"];
    sigma1=sigma*levier;
    bj=(rv_dens(rv*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
    aj=(((w.transpose())*matP).array())*((bj*p0M).transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("loglik") = -lik,
      Rcpp::Named("Pred_loglik") = pred_lik,
      Rcpp::Named("w_hat") = w_hat
    );
  return output;
}

// [[Rcpp::export]]
NumericMatrix levierVolatilityMat(const NumericMatrix &ech, const NumericMatrix Levier, const int &Nl, const NumericVector &para){
  NumericVector li(Nl);
  colvec levier=as<colvec>(wrap(Eigen::VectorXd::Ones(ech.rows())));
  for(int i=0;i<Nl;i++) li[i]=para[8]*pow(para[9],i);
  NumericVector y1;
  NumericVector y2;
  LogicalVector x;
  int t=ech.cols();
  for(int i=0; i<Nl;i++) {
    y1 = -ech(_,t-i-1);
    x=(y1>0);
    y2=1/sqrt(Levier(_,t-i-1));
    levier=levier%(1+(li[i]*as<colvec>(y1)%as<colvec>(x)%as<colvec>(y2)));
  }
  NumericVector l=wrap(levier);
  NumericMatrix Levier2=cbind(Levier,l);
  return Levier2;
}

// [[Rcpp::export]]
arma::rowvec colMeansRcpp(NumericMatrix x){
  arma::mat X = arma::mat(x.begin(), x.nrow(), x.ncol(), false); 
  return arma::mean(X, 0); 
}

// [[Rcpp::export]]
NumericMatrix Pow(const NumericMatrix &x,const double &n){
  Eigen::MatrixXd y = as<Eigen::MatrixXd>(x); 
  return wrap(y.array().pow(n)); 
}

// [[Rcpp::export]]
arma::rowvec R_hat(const int &H,
                   const Eigen::Map<Eigen::VectorXd> &ech,
                   const Eigen::Map<Eigen::MatrixXd> &z_t,
                   NumericMatrix Levier,
                   const NumericVector &para,
                   const int &N=3,
                   const int &Nl=70){
  NumericVector temp;
  NumericMatrix rt_sim(z_t.rows(),ech.size());
  rt_sim=wrap(kroneckerProduct(Eigen::VectorXd::Ones(z_t.rows()),ech.transpose()));
  for(int h=0; h<H;h++){
    Levier=levierVolatilityMat(rt_sim,Levier,Nl,para);
    temp = wrap(z_t.col(h));
    rt_sim=cbind(rt_sim,temp);
  }
  return colMeansRcpp(as<NumericMatrix>(Levier));
}

// [[Rcpp::export]]
NumericVector f_sim(const int &H,
                    const Eigen::Map<Eigen::VectorXd> &sig,
                    const Eigen::Map<Eigen::VectorXd> &pi_0,
                    const Eigen::Map<Eigen::MatrixXd> &matP){
  Eigen::MatrixXd temp = matP;
  NumericVector temp2(H);
  NumericVector temp3 =  wrap((pi_0.transpose())*sig);
  temp2[0] = temp3[0];
  for(int h=1;h<H;h++){
    temp3 = wrap((pi_0.transpose() * temp)*sig);
    temp2[h] = temp3[0];
    temp = temp * matP;
  }
  return temp2;
}


/*** R

# para<-c(0.67044,  0.99481,  5.58986,  2.33798,  0.80352, -0.31636,  0.90798, -0.03042,  0.07726,  0.14792)
# vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
# ech<-c(0.7,0.5,-0.3,-1.75,-0.7,0.5,-0.5,0.8,0.3,-0.7,0.8)
# ech<-cbind(rv=c(0.15,1.5,0.7,0.8,1.9,0.8,0.3,0.54,0.76,0.98,0.01),r=ech)
# ech1<-donnees[[1]]
# ech<-ech1[1:5012,]
# N_D<-10
# N<-3
# LEVIER<-FALSE
# if(LEVIER) {
#   para<-c(para,1.75,0.89)
#   vars<-c(vars,"l","theta")
# }
# names(para)<-vars
# 
# para_tilde<-natWork(para,LEVIER)
# 
# logLik(ech,para_tilde,LEVIER,N_D,N,Nl)
# 
# LEVIER<-TRUE
# Nl<-70
# if(LEVIER) {
#   para<-c(para,1.75,0.89)
#   vars<-c(vars,"l","theta")
# }
# names(para)<-vars
# 
# levierVolatility(ech,Nl,para)
# para_tilde<-natWork(para,LEVIER)
# workNat(para_tilde,LEVIER)
# 
# logLik(ech,para_tilde,LEVIER,N_D,N,Nl=Nl)
# 
# n<-100
# H <- 100    #nb of years simulated
# 
# Levier<-rep(1,n)%*%t(levierVolatility(ech[,2],Nl,para)$`Levier`)
# echMat<-rep(1,n)%*%t(ech[,2])
# 
# l<-logLik2(ech,para_tilde,LEVIER,N_D,N,Nl=Nl,r=0)
# 
# pi_0 <- l$w_hat
# sig <- volatilityVector(para,N_D,N)
# matP<-P(para,N_D,N)
# 
# #simulation
# set.seed(1984)
# 
# library(mhsmm)
# MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
# z_t<-matrix(rnorm(n*H),nrow=n,ncol=H)
# r_t2<-R_hat(H,ech[,2],MC_sim,z_t,Levier,sig,para,N,Nl)
# rvt_sim <- r_t2$rvt
# rt2_sim  <- r_t2$rt2


*/

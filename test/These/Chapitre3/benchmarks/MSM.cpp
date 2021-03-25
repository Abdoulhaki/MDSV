#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
Eigen::VectorXd workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,
                        const bool &LEVIER){
  Eigen::VectorXd para=para_tilde;
  
  para(0)=1/(1+exp(para_tilde(0)));   //m0
  para(1)=exp(para_tilde(1));         //sigma
  para(2)=1+exp(para_tilde(2));       //b
  para(3)=1/(1+exp(para_tilde(3)));   //gamma1
  if(LEVIER){
    para(4)=exp(para_tilde(4)); //l
    para(5)=1/(1+exp(para_tilde(5))); //theta
  }
  return para;
}

// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para,
                        const bool &LEVIER){
  Eigen::VectorXd para_tilde=para;
  
  para_tilde[0]=log((1/para[0])-1);  //m0
  para_tilde[1]=log((para[1]));      //sigma
  para_tilde[2]=log(para[2]-1);      //b
  para_tilde[3]=log((1/para[3])-1);  //gamma
  if(LEVIER){
    para_tilde[4]=log(para[4]); //l
    para_tilde[5]=log((1/para[5])-1); //theta
  }
  return para_tilde;
}

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector(const Eigen::VectorXd &para, const int &N){
  Eigen::VectorXd sigma_i=Eigen::VectorXd::Ones(2);
  sigma_i[0]=para[0];
  sigma_i[1]=2-para[0];
  Eigen::VectorXd sigma=sigma_i;
  for(int i(1);i<N;i++ ) sigma=kroneckerProduct(sigma,sigma_i).eval();
  sigma=para[1]*sigma;
  return sigma;
}

// [[Rcpp::export]]
Eigen::MatrixXd transitionMatrix(const double &gamma){
  return (1-gamma)*(Eigen::MatrixXd::Identity(2,2))+(gamma/2)*(Eigen::VectorXd::Ones(2))*(Eigen::VectorXd::Ones(2)).transpose();
}

// [[Rcpp::export]]
mat P(const Eigen::VectorXd &para,const int &N){
  colvec gamma(N);
  gamma[0]=para[3];
  mat P=as<mat>(wrap(transitionMatrix(gamma[0])));
  if(N==1) return(P);
  for(int k=1; k<N; k++) {
    gamma[k]=1-pow(1-para[3],(pow(para[2],k)));
    P=kron(P,as<mat>(wrap(transitionMatrix(gamma[k])))).eval();
  }
  return(P);
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[4]*pow(para[5],i);
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
double logLik(const NumericVector &ech,
              const Eigen::Map<Eigen::VectorXd> &para_tilde,
              const bool &LEVIER,
              const int &N,
              const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n(ech.size());
  Eigen::VectorXd sigma = volatilityVector(para,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = as<Eigen::MatrixXd>(wrap(P(para,N)));
  Eigen::VectorXd Res;
  aj=r_dens(ech(0),sigma).array();
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }

  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const NumericVector &ech,
              const Eigen::Map<Eigen::VectorXd> &para_tilde,
              const bool &LEVIER,
              const int &N,
              const double &r=0,
              const int &t=2,
              const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n(ech.size());
  Eigen::VectorXd sigma = volatilityVector(para,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  double a(0);
  double lik(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  Eigen::MatrixXd matP = as<Eigen::MatrixXd>(wrap(P(para,N)));
  Eigen::VectorXd Res;
  aj=r_dens(ech(0),sigma).array();
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
      
      if(i==(t-1)) w_hat=w;
    }
    Res=r_dens(r,sigma).array();
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  
  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
      
      if(i==(t-1)) w_hat=w;
    }
    double levier=L["levier"];
    sigma1=sigma*levier;
    Res=r_dens(r,sigma1).array();
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
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
  for(int i=0;i<Nl;i++) li[i]=para[4]*pow(para[5],i);
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
Rcpp::List R_hat(const int &H,
                 const Eigen::Map<Eigen::VectorXd> &ech,
                 const Eigen::Map<Eigen::MatrixXi> &MC_sim,
                 const Eigen::Map<Eigen::MatrixXd> &z_t,
                 NumericMatrix Levier,
                 const NumericVector &sig,
                 const NumericVector &para,
                 const int &N=3,
                 const int &Nl=70){
  NumericVector temp;
  NumericVector temp1;
  NumericVector temp2;
  NumericMatrix rt_sim(z_t.rows(),ech.size());
  rt_sim=wrap(kroneckerProduct(Eigen::VectorXd::Ones(z_t.rows()),ech.transpose()));
  for(int h=0; h<H;h++){
    Levier=levierVolatilityMat(rt_sim,Levier,Nl,para);
    IntegerVector idx=wrap(MC_sim.col(h).array()-1);
    temp1=sig[idx];
    temp2=wrap(z_t.col(h));
    temp=sqrt(Levier(_,Levier.ncol()-1)*temp1)*temp2;
    rt_sim=cbind(rt_sim,temp);
  }
  NumericMatrix rt2 = Pow(rt_sim,2);
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt_sim") = rt_sim,
      Rcpp::Named("rt2") = colMeansRcpp(as<NumericMatrix>(rt2))
    );
  return output;
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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

# MSM

# logR<-c(1.2,1.5,-2.3,-0.5,0.7,-0.8,0.8,0.8,0.7,-0.5,-1.3)
# NMSM <- 5
# 
# para<-c(0.8,1.26,5,0.05)
# para_tilde<-natWork(para,FALSE)
# volatilityVector(para,NMSM)
# logLik(logR,para_tilde,FALSE,NMSM)
# nlm(logLik,ech=logR,para_tilde,LEVIER=FALSE,N=NMSM)
# 
# para<-c(0.8,1.26,5,0.05,1.89,0.89)
# para_tilde<-natWork(para,TRUE)
# volatilityVector(para,NMSM)
# levierVolatility(logR,3,para)
# logLik(logR,para_tilde,TRUE,NMSM,Nl=10)
# s<-nlm(logLik,ech=logR,para_tilde,LEVIER=TRUE,N=NMSM,Nl=3)
# 
# workNat(s$estimate,TRUE)[1:4]

***/
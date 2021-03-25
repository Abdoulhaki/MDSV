#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::MatrixXd r_dens(const double &x, const Eigen::MatrixXd &sd){
  int n = sd.rows();
  int m = sd.cols();
  Eigen::MatrixXd res(n,m);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      res(i,j) = R::dnorm(x, 0, sqrt(sd(i,j)), FALSE);
    }
  }
  return(res);
}

// Function that transform the working parameters to the natural parameters of FHMV
// [[Rcpp::export]]
NumericVector workNatFHMV(const NumericVector &para_tilde,const bool &LEVIER){
  int taillePara(para_tilde.size());
  NumericVector para(taillePara);
    para[0]=exp(para_tilde[0]);//sigma
    para[1]=1+exp(para_tilde[1]);//c1
    para[2]=1/(1+exp(para_tilde[2]));//thetac
    para[3]=1/(1+exp(para_tilde[3]));//p
    para[4]=1+exp(para_tilde[4]);//m1
    para[5]=1/(1+exp(para_tilde[5]));//thetam
    para[6]=1/(1+exp(para_tilde[6]));//q
  if(LEVIER){ 
    para[7]=exp(para_tilde[7]);//l
    para[8]=1/(1+exp(para_tilde[8]));//thetal
  }  
  return para;
}

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
  if(LEVIER){ 
    para[7]=exp(para_tilde[7]);      //l
    para[8]=1/(1+exp(para_tilde[8]));//thetal
  }
  return para;
}

// Function that transform the natural parameters to the working parameters of FHMV
// [[Rcpp::export]]
NumericVector natWorkFHMV(const NumericVector &para,const bool &LEVIER){
  int taillePara(para.size());
  NumericVector para_tilde(taillePara);
    para_tilde[0]=log(para[0]);//sigma
    para_tilde[1]=log(para[1]-1);//c1
    para_tilde[2]=log((1/para[2])-1);//theta_c
    para_tilde[3]=log((1/para[3])-1);//p
    para_tilde[4]=log(para[4]-1);//m1
    para_tilde[5]=log((1/para[5])-1);//theta_m
    para_tilde[6]=log((1/para[6])-1);//q
  if(LEVIER){ 
    para_tilde[7]=log(para[7]);//l
    para_tilde[8]=log((1/para[8])-1);//thetal
  }
  return para_tilde;
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
  if(LEVIER){ 
    para_tilde[7]=log(para[7]);      //l
    para_tilde[8]=log((1/para[8])-1);//thetal
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
Eigen::MatrixXd P(const Eigen::VectorXd &para,const int &N){
  Eigen::MatrixXd P=transitionMatrix(para[3]);
  if(N==1) return(P);
  for(int k=1; k<N; k++) P=kroneckerProduct(P,transitionMatrix(para[3])).eval();
  return(P);
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){
  int T=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(T));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[7]*pow(para[8],i);
  for(int t=Nl;t<ech.size();t++){
    levier=1;
    for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
    Levier[t]=levier;
  }
  levier=1;
  for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[T-i-1]))*(ech[T-i-1]<0)/(sqrt(Levier[T-i-1]))));
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
  Eigen::MatrixXd sigma = volatilityVector(para,N);
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  Eigen::VectorXd p0M=(para[6]/(N-1))*(Eigen::VectorXd::Ones(N));
  p0M[N-1]=1-para[6];
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = P(para,N);
  Eigen::VectorXd Res;
  Eigen::VectorXd aj=r_dens(ech(0),sigma)*p0M;
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma)*p0M;
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }

  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::MatrixXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1)*p0M;
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
  Eigen::MatrixXd sigma = volatilityVector(para,N);
  Eigen::VectorXd p0=(1.0/pow(2,N))*(Eigen::VectorXd::Ones(pow(2,N)));
  Eigen::VectorXd p0M=(para[6]/(N-1))*(Eigen::VectorXd::Ones(N));
  p0M[N-1]=1-para[6];
  double a(0);
  double lik(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  Eigen::MatrixXd matP = P(para,N);
  Eigen::VectorXd Res;
  Eigen::VectorXd aj=r_dens(ech(0),sigma)*p0M;
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma)*p0M;
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
      
      if(i==(t-1)) w_hat=w;
    }
    Res=r_dens(r,sigma)*p0M;
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  
  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::MatrixXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1)*p0M;
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
      
      if(i==(t-1)) w_hat=w;
    }
    double levier=L["levier"];
    sigma1=sigma*levier;
    Res=r_dens(r,sigma1)*p0M;
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
  for(int i=0;i<Nl;i++) li[i]=para[7]*pow(para[8],i);
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


# # FHMV
# 
# N <- 3
# logR<-c(1.2,1.5,-2.3,-0.5,0.7,-0.8,0.8,0.8,0.7,-0.5,-1.3)
# 
# para<-c(sqrt(0.22),1.99,0.51,0.9986,23.55,0.87,0.93) # Essai FHMV simple
# para_tilde<-natWork(para,FALSE)
# # volatilityVector(para,3)
# # P(para[4],N)
# logLik(logR,para_tilde,LEVIER=FALSE,N=N)
# nlm(logLik,ech=logR,para_tilde,N=N,LEVIER=FALSE,Nl=3)
# para<-c(sqrt(0.22),1.99,0.51,0.9986,23.55,0.87,0.93,1,0.92) # Essai FHMV simple
# para_tilde<-natWork(para,TRUE)
# logLik(logR,para_tilde,LEVIER=TRUE,N,Nl=3)
# levierVolatility(ech=logR,Nl=3,para)
# logLik(logR,para_tilde,LEVIER=TRUE,N=N,Nl=3)
# nlm(logLik,ech=logR,para_tilde,N=N,LEVIER=TRUE,Nl=3)
#  solnp(para0_tildeFHMV,logLikFHMV,ech=logR,density="normal",modelType="levier",N=NFHMV)

***/
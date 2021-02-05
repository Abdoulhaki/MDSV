#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,
                   const bool &LEVIER=false, const int &Model_type=0){
  Eigen::VectorXd para=para_tilde;

  para(0)=1/(1+exp(para_tilde(0)));      //omega
  para(1)=1/(1+exp(para_tilde(1)));      //a
  para(2)=1+exp(para_tilde(2));          //b
  para(3)=exp(para_tilde(3));            //sigma
  para(4)=1/(1+exp(para_tilde(4)));      //v_0
  
  int j=0;
  if(Model_type==1){
    para(5)=exp(para_tilde(5));            //shape
    j=1;
  }
  if(Model_type==2){
    para(5)=para_tilde(5);                 //xi
    para(6)=para_tilde(6);                 //varphi
    para(7)=para_tilde(7);                 //delta_1
    para(8)=para_tilde(8);                 //delta_2
    para(9)=exp(para_tilde(9));            //shape
    j=5;
  }
  if(LEVIER){
    para(5+j)=exp(para_tilde(5+j));        //l
    para(6+j)=1/(1+exp(para_tilde(6+j)));  //theta
  }
  return para;
}

// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para,
                              const bool &LEVIER=false, const int &Model_type=0){
  Eigen::VectorXd para_tilde=para;

  para_tilde[0]=log((1/para[0])-1);      //omega
  para_tilde[1]=log((1/para[1])-1);      //a
  para_tilde[2]=log(para[2]-1);      //b
  para_tilde[3]=log(para[3]);            //sigma
  para_tilde[4]=log((1/para[4])-1);    //v_0

  int j=0;
  if(Model_type==1){
    para_tilde[5]=log(para[5]);//shape
    j=1;
  }
  if(Model_type==2){
    para_tilde[5]=para[5];    //xi
    para_tilde[6]=para[6];    //varphi
    para_tilde[7]=para[7];    //delta_1
    para_tilde[8]=para[8];    //delta_2
    para_tilde[9]=log(para[9]);//shape
    j=5;
  }
  if(LEVIER){
    para_tilde[5+j]=log(para[5+j]); //l
    para_tilde[6+j]=log((1/para[6+j])-1); //theta
  }
  return para_tilde;
}


NumericVector proba_pi(const double &omega, const int &K){
  NumericVector probaPi=dbinom(seq(0,K-1),K-1,omega);
  return probaPi;
}

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector(const Eigen::VectorXd &para, const int &K, const int &N){
  Eigen::VectorXd sigma=Eigen::VectorXd::Ones(1);
  Eigen::VectorXd sigma_i(K);
  Eigen::VectorXd probaPi=as<Eigen::VectorXd>(proba_pi(para[0],K));
  double e_i;
  for(int k=0; k<K;k++) sigma_i[k]=para[4]*pow(((2-para[4])/para[4]),k);
  e_i=(probaPi.transpose())*sigma_i;
  sigma_i=sigma_i/e_i;
  for(int i(0);i<N;i++ ) sigma=kroneckerProduct(sigma,sigma_i).eval();
  sigma=para[3]*sigma;
  return sigma;
}

// [[Rcpp::export]]
Eigen::VectorXd probapi(const double &omega, const int &K, const int &N){
  Eigen::VectorXd probaPi=as<Eigen::VectorXd>(proba_pi(omega,K));
  Eigen::VectorXd proba=as<Eigen::VectorXd>(proba_pi(omega,K));
  if(N==1) return(probaPi);
  for(int i(1);i<N;i++) probaPi=kroneckerProduct(probaPi,proba).eval();
  return(probaPi);
}

Eigen::MatrixXd transitionMatrix(const double &phi, const double &omega, const int &K){
  return phi*(Eigen::MatrixXd::Identity(K,K))+(1-phi)*(Eigen::VectorXd::Ones(K))*(as<Eigen::VectorXd>(proba_pi(omega,K))).transpose();
}

// [[Rcpp::export]]
Eigen::MatrixXd P(const Eigen::VectorXd &para, const int &K,const int &N){
  Eigen::VectorXd phi(N);
  phi[0]=para[1];
  for(int k=1; k<N; k++) phi[k]=pow(para[1],(pow(para[2],k)));
  Eigen::MatrixXd P=transitionMatrix(phi[0],para[0],K);
  if(N==1) return(P);
  for(int i(1);i<N;i++) P=kroneckerProduct(P,transitionMatrix(phi[i],para[0],K)).eval();
  return(P);
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const Eigen::VectorXd &para, const int &Nl=70,
                            const int &Model_type=0){
  int j=0; if(Model_type==1)  j=1; if(Model_type==2)  j=5;
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[5+j]*pow(para[6+j],i);
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

Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

Eigen::VectorXd rv_dens(const Eigen::VectorXd &x,const double &shape){
  int n = x.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i]=R::dgamma(x[i],shape,1.0/shape,false);
  }
  return(res);
}

Eigen::VectorXd jrv_dens(const double &x,const double &shape,const Eigen::VectorXd &mu_rv){
  int n = mu_rv.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i]=R::dlnorm(x,mu_rv[i],sqrt(shape),FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
double logLik(const Eigen::MatrixXd &ech,
              const Eigen::Map<Eigen::VectorXd> &para_tilde, 
              const int &Model_type=0,
              const bool &LEVIER=false,
              const int &K=2,
              const int &N=2,
              const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER,Model_type);
  int n=ech.rows();
  int k=ech.cols();
  Eigen::VectorXd sigma = volatilityVector(para,K,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=probapi(para[0],K,N);
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  double x=ech(0,0);
  double y=0;
  if(k>1) y=ech(0,1);
  Eigen::MatrixXd matP = P(para,K,N);
  Eigen::VectorXd Res;

  Eigen::VectorXd temp;
  Eigen::VectorXd mu_rv;
  double xi(0);
  double varphi(0);
  double delta1(0);
  double delta2(0);
  double shape(0);
  if(Model_type==2){
    xi=para[5];
    varphi=para[6];
    delta1=para[7];
    delta2=para[8];
    shape=para[9];
    temp=x*(sigma.array().sqrt().inverse());
    mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
    aj=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
  }else if(Model_type==1){
    shape=para[5];
    aj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
  }else{
    aj=r_dens(x,sigma);
  }
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  if(!LEVIER){
    if(Model_type==2){
      for(int i(1); i<n;i++){
        x=ech(i,0);
        y=ech(i,1);
        temp=x*(sigma.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else if(Model_type==1){
      for(int i(1); i<n;i++){
        y=ech(i,1);
        Res=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else {
      for(int i(1); i<n;i++){
        x=ech(i,0);
        Res=r_dens(x,sigma);
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(0));
    Rcpp::List L=levierVolatility(ech1,para,Nl,Model_type);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    if(Model_type==2){
      for(int i(1); i<n;i++){
        x=ech(i,0);
        y=ech(i,1);
        sigma1=sigma*Levier[i];
        temp=x*(sigma1.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma1).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else if(Model_type==1){
      for(int i(1); i<n;i++){
        y=ech(i,1);
        sigma1=sigma*Levier[i];
        Res=(rv_dens(y*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else {
      for(int i(1); i<n;i++){
        x=ech(i,0);
        sigma1=sigma*Levier[i];
        Res=r_dens(x,sigma1);
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }
  }
  
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const Eigen::MatrixXd &ech, 
              const Eigen::VectorXd &para,
              const int &Model_type=0,
              const bool &LEVIER=false,
              const int &K=2,
              const int &N=2,
              const double &r=0,
              const int &t=2,
              const int &Nl=70){
  int n=ech.rows();
  int k=ech.cols();
  Eigen::VectorXd sigma = volatilityVector(para,K,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=probapi(para[0],K,N);
  double a(0);
  double lik(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  double x=ech(0,0);
  double y=0;
  if(k>1) y=ech(0,1);
  Eigen::MatrixXd matP = P(para,K,N);
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(n));
  Eigen::VectorXd Res;
  
  Eigen::VectorXd ajm;
  double likm(0);
  Eigen::VectorXd Resm;
  double xi(0);
  double varphi(0);
  double delta1(0);
  double delta2(0);
  double shape(0);
  
  int tmp =pow(K,N);
  Eigen::MatrixXd proba1(tmp,n);
  Eigen::MatrixXd proba2(tmp,n);
  
  Eigen::VectorXd temp;
  Eigen::VectorXd mu_rv;
  if(Model_type==2){
    xi=para[5];
    varphi=para[6];
    delta1=para[7];
    delta2=para[8];
    shape=para[9];
    temp=x*(sigma.array().sqrt().inverse());
    mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
    aj=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
    
    ajm=r_dens(x,sigma);
    ajm=p0.array()*ajm.array();
    likm=log(ajm.sum());
  }else if(Model_type==1){
    shape=para[5];
    aj=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
  }else{
    aj=r_dens(x,sigma);
  }
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  proba1.col(0) = w;
  
  if(!LEVIER){
    if(Model_type==2){
      for(int i(1); i<n;i++){
        x=ech(i,0);
        y=ech(i,1);
        temp=x*(sigma.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
        
        Resm=r_dens(x,sigma);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
      Res=r_dens(r,sigma);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }else if(Model_type==1){
      for(int i(1); i<n;i++){
        y=ech(i,1);
        Res=(rv_dens(y*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
      }
      Res=(rv_dens(r*sigma.array().inverse(),shape).array())*(sigma.array().inverse());
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }else {
      for(int i(1); i<n;i++){
        x=ech(i,0);
        Res=r_dens(x,sigma);
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
      }
      Res=r_dens(r,sigma);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(0));
    Rcpp::List L=levierVolatility(ech1,para,Nl,Model_type);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    if(Model_type==2){
      for(int i(1); i<n;i++){
        x=ech(i,0);
        y=ech(i,1);
        sigma1=sigma*Levier[i];
        temp=x*(sigma1.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma1).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
        
        Resm=r_dens(x,sigma1);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
      double levier=L["levier"];
      sigma1=sigma*levier;
      Res=r_dens(r,sigma1);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }else if(Model_type==1){
      for(int i(1); i<n;i++){
        y=ech(i,1);
        sigma1=sigma*Levier[i];
        Res=(rv_dens(y*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
      }
      double levier=L["levier"];
      sigma1=sigma*levier;
      Res=(rv_dens(r*sigma1.array().inverse(),shape).array())*(sigma1.array().inverse());
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }else {
      for(int i(1); i<n;i++){
        x=ech(i,0);
        sigma1=sigma*Levier[i];
        Res=r_dens(x,sigma1);
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        proba1.col(i) = w;
      }
      double levier=L["levier"];
      sigma1=sigma*levier;
      Res=r_dens(r,sigma1);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      pred_lik=log(a);
    }
  }
  w_hat = proba1.col(t-1);
  
  proba2.col(n-1) = w;
  
  for(int i=(n-2); i>=0;--i){
    temp = (proba1.col(i).transpose()*matP).array().inverse();
    temp = (proba2.col(i+1).array())*(temp.array());
    temp = matP*temp;
    proba2.col(i) = proba1.col(i).array()*temp.array();
  }

  Rcpp::List output =
        Rcpp::List::create(
          Rcpp::Named("loglik") = -lik,
          Rcpp::Named("Marg_loglik") = likm,
          Rcpp::Named("Pred_loglik") = pred_lik,
          Rcpp::Named("w_hat") = w_hat,
          Rcpp::Named("Levier") = Levier,
          Rcpp::Named("filtred_proba") = proba1,
          Rcpp::Named("smoothed_proba") = proba2
        );
  return output;
}

// [[Rcpp::export]]
NumericMatrix levierVolatilityMat(const NumericMatrix &ech, const NumericMatrix Levier, 
                                  const NumericVector &para, const int &Nl=70, const int &Model_type=0){
  int j=0; if(Model_type==1)  j=1; if(Model_type==2)  j=5;
  arma::colvec levier=as<arma::colvec>(wrap(Eigen::VectorXd::Ones(ech.rows())));
  NumericVector li(Nl);
  for(int i=0;i<Nl;i++) li[i]=para[5+j]*pow(para[6+j],i);
  NumericVector y1;
  NumericVector y2;
  LogicalVector x;
  int t=ech.cols();
  for(int i=0; i<Nl;i++) {
    y1 = -ech(_,t-i-1);
    x=(y1>0);
    y2=1/sqrt(Levier(_,t-i-1));
    levier=levier%(1+(li[i]*as<arma::colvec>(y1)%as<arma::colvec>(x)%as<arma::colvec>(y2)));
  }
  NumericVector l=wrap(levier);
  NumericMatrix Levier2=cbind(Levier,l);
  return Levier2;
}

arma::rowvec colMeansRcpp(NumericMatrix x){
  arma::mat X = arma::mat(x.begin(), x.nrow(), x.ncol(), false); 
  return arma::mean(X, 0); 
}

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
                    const int &Model_type=0,
                    const int &N=3,
                    const int &Nl=70){
  NumericVector temp;
  NumericVector temp1;
  NumericVector temp2;
  NumericMatrix rt_sim(z_t.rows(),ech.size());
  NumericMatrix rvt_sim(z_t.rows(),ech.size());
  rt_sim=wrap(kroneckerProduct(Eigen::VectorXd::Ones(z_t.rows()),ech.transpose()));
  NumericMatrix rvt;
  Rcpp::List output;
  if(Model_type==1){
    for(int h=0; h<H;h++){
      Levier=levierVolatilityMat(rt_sim,Levier,para,Nl,Model_type);
      temp = wrap(z_t.col(h));
      rt_sim=cbind(rt_sim,temp);
    }
  }else{
    for(int h=0; h<H;h++){
      Levier=levierVolatilityMat(rt_sim,Levier,para,Nl,Model_type);
      IntegerVector idx=wrap(MC_sim.col(h).array()-1);
      temp1=sig[idx];
      temp2=wrap(z_t.col(h));
      temp=sqrt(Levier(_,Levier.ncol()-1)*temp1)*temp2;
      rt_sim=cbind(rt_sim,temp);
    }
  }
  
  NumericMatrix rt2 = Pow(rt_sim,2);
  if(Model_type==2){
    rvt     = Pow(rt2,para[6])*exp(para[5]+0.5*para[9]+0.5*pow(para[7],2)-para[8])/sqrt(1-2*para[8]);
  }
  
  if(Model_type==1){
    output = Rcpp::List::create(
        Rcpp::Named("LevierMat") = Levier,
        Rcpp::Named("Levier") = colMeansRcpp(as<NumericMatrix>(Levier))
      );
  }else if(Model_type==0){
    output = Rcpp::List::create(
        Rcpp::Named("rt_sim") = rt_sim,
        Rcpp::Named("rt2") = colMeansRcpp(as<NumericMatrix>(rt2))
      );
  }else if(Model_type==2){
    output = Rcpp::List::create(
        Rcpp::Named("rt_sim") = rt_sim,
        Rcpp::Named("rvt_sim") = rvt,
        Rcpp::Named("rt2") = colMeansRcpp(as<NumericMatrix>(rt2)),
        Rcpp::Named("rvt") = colMeansRcpp(as<NumericMatrix>(rvt))
      );
  }
  
  return output;
}

// [[Rcpp::export]]
Rcpp::List f_sim(const int &H,
                      const Eigen::Map<Eigen::VectorXd> &sig,
                      const Eigen::Map<Eigen::VectorXd> &pi_0,
                      const Eigen::Map<Eigen::MatrixXd> &matP,
                      const double &varphi=0,
                      const double &xi=0,
                      const double &shape=0,
                      const double &delta1=0,
                      const double &delta2=0){
  Eigen::MatrixXd temp = matP;
  NumericVector temp2_r(H);
  NumericVector temp2_rv(H);
  NumericVector temp3_r =  wrap((pi_0.transpose())*sig);
  Eigen::VectorXd sig2=sig.array().pow(varphi);
  NumericVector temp3_rv =  wrap((pi_0.transpose())*sig2);
  temp2_r[0] = temp3_r[0];
  temp2_rv[0] = temp3_rv[0]*(exp(xi+0.5*shape+0.5*pow(delta1,2)-delta2)/sqrt(1-2*delta2));
  for(int h=1;h<H;h++){
    temp3_r = wrap((pi_0.transpose() * temp)*sig);
    temp3_rv = wrap((pi_0.transpose() * temp)*sig2);
    temp2_r[h] = temp3_r[0];
    temp2_rv[h] = temp3_rv[0]*(exp(xi+0.5*shape+0.5*pow(delta1,2)-delta2)/sqrt(1-2*delta2));
    temp = temp * matP;
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt2") = temp2_r,
      Rcpp::Named("rvt") = temp2_rv
    );
  return output;

}



/*** R

*/

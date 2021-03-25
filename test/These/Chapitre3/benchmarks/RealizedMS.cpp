#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,const bool &LEVIER, 
                   const int &K, const std::string &Model_type){
  
  int logN=0;
  if(Model_type=="logN") logN=4;
  int l=0;
  if(LEVIER) l=2;
  Eigen::VectorXd para(1+l+logN);
  Eigen::MatrixXd P(K,K);
  Eigen::VectorXd Vol(K);
  Eigen::VectorXd probapi(K);
  double temp=0;

  if(Model_type=="logN"){
    para(0)=para_tilde(0);                 //xi
    para(1)=para_tilde(1);                 //varphi
    para(2)=para_tilde(2);                 //delta_1
    para(3)=para_tilde(3);                 //delta_2
  }
  para(logN)=exp(para_tilde(logN));        //shape
  
  if(LEVIER){
    para(logN+1)=exp(para_tilde(logN+1)); //l
    para(logN+2)=1/(1+exp(para_tilde(logN+2))); //theta
  }

  for(int k=0; k<K; k++){
    Vol(k)=exp(para_tilde(1+l+logN+k));   //volatility
  }

  for(int i=0; i<K; i++){
    for(int j=0; j<(K-1); j++){
    temp=temp+exp(para_tilde(1+l+logN+K+i*(K-1)+j));
    }
    for(int j=0; j<(K-1); j++){
      P(i,j)=exp(para_tilde(1+l+logN+K+i*(K-1)+j))/(1+temp);      //pij
    }
    P(i,(K-1))=1/(1+temp);
    temp=0;
  }

  Eigen::MatrixXd A = (Eigen::MatrixXd::Identity(K,K)-P+Eigen::MatrixXd::Ones(K,K)).transpose();
  probapi=A.colPivHouseholderQr().solve(Eigen::VectorXd::Ones(K));
  
    Rcpp::List output =
      Rcpp::List::create(
        Rcpp::Named("para") = para,
        Rcpp::Named("P") = P,
        Rcpp::Named("Vol") = Vol,
        Rcpp::Named("probapi") = probapi
      );
    return output;
}


// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para,const Eigen::Map<Eigen::VectorXd> &Vol,
                        const Eigen::Map<Eigen::MatrixXd> &P,const bool &LEVIER, 
                        const int &K, const std::string &Model_type){
  int logN=0;
  if(Model_type=="logN") logN=4;
  int l=0;
  if(LEVIER) l=2;
  Eigen::VectorXd para_tilde(1+l+logN+K*K);

  if(Model_type=="logN"){
    para_tilde[0]=para[0];    //xi
    para_tilde[1]=para[1];    //varphi
    para_tilde[2]=para[2];    //delta_1
    para_tilde[3]=para[3];    //delta_2
  }
  para_tilde[logN]=log(para[logN]);//shape
  if(LEVIER){
    para_tilde[1+logN]=log(para[1+logN]); //l
    para_tilde[2+logN]=log((1/para[2+logN])-1); //theta
  }

  for(int k=0; k<K; k++){
    para_tilde[1+l+logN+k]=log(Vol(k));            //volatility
  }

  for(int i=0; i<K; i++){
    for(int j=0; j<(K-1); j++){
      para_tilde[1+l+logN+K+i*(K-1)+j]=log((P(i,j))/(P(i,K-1)));      //pij
    }
  }
  return para_tilde;
}

// // [[Rcpp::export]]
// NumericVector proba_pi(const double &omega, const int &N_D){
//   NumericVector probaPi=dbinom(seq(0,N_D-1),N_D-1,omega);
//   return probaPi;
// }

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, 
                            const Eigen::VectorXd &para, const std::string &Model_type){
  int logN=0;
  if(Model_type=="logN") logN=4;
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[logN+1]*pow(para[logN+2],i);
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
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
Eigen::VectorXd rv_dens(const double &x,const double &shape, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i]=R::dgamma((1.0/x),shape+1,1.0/(shape*sd[i]),false)/pow(x,2);
  }
  return(res);
}

// [[Rcpp::export]]
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
                    const bool &LEVIER,
                    const int &K,
                    const std::string &Model_type="InvG",
                    const int &Nl=70){
  Rcpp::List para=workNat(para_tilde,LEVIER,K,Model_type);
  Eigen::VectorXd opara=para["para"];
  int n=ech.rows();
  // double shape=opara[0];
  Eigen::VectorXd sigma = para["Vol"];
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=para["probapi"];
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  double x=ech(0,1);
  double y=ech(0,0);
  Eigen::MatrixXd matP = para["P"];
  Eigen::VectorXd Res;
  
  Eigen::VectorXd temp;
  Eigen::VectorXd mu_rv;
  double xi(0);
  double varphi(0);
  double delta1(0);
  double delta2(0);
  double shape(0);
  if(Model_type=="logN"){
    xi=opara[0];
    varphi=opara[1];
    delta1=opara[2];
    delta2=opara[3];
    shape=opara[4];
    temp=x*(sigma.array().sqrt().inverse());
    mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
    aj=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
  }else{
    shape=opara[0];
    aj=r_dens(x,sigma).array()*rv_dens(y,shape,sigma).array();
  }
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    if(Model_type=="logN"){
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        temp=x*(sigma.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else{
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        Res=r_dens(x,sigma).array()*rv_dens(y,shape,sigma).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        // if(a==0) {
        //   a=sigma.maxCoeff();
        //   aj=(sigma.array()==a).cast<double>();
        //   a=aj.sum();
        // }
        lik+=log(a);
        w=aj/a;
      }
    }
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,opara,Model_type);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    
    if(Model_type=="logN"){
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        sigma1=sigma*Levier[i];
        temp=x*(sigma1.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma1).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
      }
    }else{
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        sigma1=sigma*Levier[i];
        Res=r_dens(x,sigma1).array()*rv_dens(y,shape,sigma1).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        // if(a==0) {
        //   a=sigma.maxCoeff();
        //   aj=(sigma.array()==a).cast<double>();
        //   a=aj.sum();
        // }
        lik+=log(a);
        w=aj/a;
      }
    }
  }
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const Eigen::MatrixXd &ech,
                   const Eigen::Map<Eigen::VectorXd> &para_tilde,
                   const bool &LEVIER,
                   const int &K,
                   const std::string &Model_type="InvG",
                   const double &r=0,
                   const int &t=2,
                   const int &Nl=70){
  Rcpp::List para=workNat(para_tilde,LEVIER,K,Model_type);
  Eigen::VectorXd opara=para["para"];
  int n=ech.rows();
  // double shape=opara[0];
  Eigen::VectorXd sigma = para["Vol"];
  Eigen::VectorXd aj;
  Eigen::VectorXd ajm;
  Eigen::VectorXd p0=para["probapi"];
  double a(0);
  double lik(0);
  double likm(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  double x=ech(0,1);
  double y=ech(0,0);
  Eigen::MatrixXd matP = para["P"];
  Eigen::VectorXd Res;
  Eigen::VectorXd Resm;
  
  Eigen::VectorXd temp;
  Eigen::VectorXd mu_rv;
  double xi(0);
  double varphi(0);
  double delta1(0);
  double delta2(0);
  double shape(0);
  if(Model_type=="logN"){
    xi=opara[0];
    varphi=opara[1];
    delta1=opara[2];
    delta2=opara[3];
    shape=opara[4];
    temp=x*(sigma.array().sqrt().inverse());
    mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
    aj=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
  }else{
    shape=opara[0];
    aj=r_dens(x,sigma).array()*rv_dens(y,shape,sigma).array();
  }
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  ajm=r_dens(x,sigma);
  ajm=p0.array()*ajm.array();
  likm=log(ajm.sum());

  if(!LEVIER){
    
    if(Model_type=="logN"){
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        temp=x*(sigma.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        
        if(i==(t-1)) w_hat=w;
        
        Resm=r_dens(x,sigma);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
    }else{
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        Res=r_dens(x,sigma).array()*rv_dens(y,shape,sigma).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        // if(a==0) {
        //   a=sigma.maxCoeff();
        //   aj=(sigma.array()==a).cast<double>();
        //   a=aj.sum();
        // }
        lik+=log(a);
        w=aj/a;
        
        if(i==(t-1)) w_hat=w;
        
        Resm=r_dens(x,sigma);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
    }
    
    Res=r_dens(r,sigma);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,opara,Model_type);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    
    if(Model_type=="logN"){
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        sigma1=sigma*Levier[i];
        temp=x*(sigma1.array().sqrt().inverse());
        mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
        Res=r_dens(x,sigma1).array()*jrv_dens(y,shape,mu_rv).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        lik+=log(a);
        w=aj/a;
        
        if(i==(t-1)) w_hat=w;
        
        Resm=r_dens(x,sigma1);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
    }else{
      for(int i(1); i<n;i++){
        x=ech(i,1);
        y=ech(i,0);
        sigma1=sigma*Levier[i];
        Res=r_dens(x,sigma1).array()*rv_dens(y,shape,sigma1).array();
        aj=(((w.transpose())*matP).array())*(Res.transpose().array());
        a=aj.sum();
        // if(a==0) {
        //   a=sigma.maxCoeff();
        //   aj=(sigma.array()==a).cast<double>();
        //   a=aj.sum();
        // }
        lik+=log(a);
        w=aj/a;
        
        if(i==(t-1)) w_hat=w;
        
        Resm=r_dens(x,sigma1);
        ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
        likm+=log(ajm.sum());
      }
    }
    double levier=L["levier"];
    sigma1=sigma*levier;
    Res=r_dens(r,sigma1);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("loglik") = -lik,
      Rcpp::Named("Marg_loglik") = likm,
      Rcpp::Named("Pred_loglik") = pred_lik,
      Rcpp::Named("w_hat") = w_hat
    );
  return output;
}

// [[Rcpp::export]]
NumericMatrix levierVolatilityMat(const NumericMatrix &ech, const NumericMatrix Levier, const int &Nl,
                                  const NumericVector &para, const std::string &Model_type="InvG"){
  int logN=0;
  if(Model_type=="logN") logN=4;
  NumericVector li(Nl);
  colvec levier=as<colvec>(wrap(Eigen::VectorXd::Ones(ech.rows())));
  for(int i=0;i<Nl;i++) li[i]=para[logN+1]*pow(para[logN+2],i);
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
                    const Eigen::Map<Eigen::MatrixXd> &n_t,
                    NumericMatrix Levier,
                    const NumericVector &sig,
                    const NumericVector &para,
                    const int &K,
                    const std::string &Model_type="InvG",
                    const int &Nl=70){
  NumericVector temp;
  NumericVector temp1;
  NumericVector temp2;
  NumericMatrix rt_sim(z_t.rows(),ech.size());
  NumericMatrix rvt(n_t.rows(),H);
  rt_sim=wrap(kroneckerProduct(Eigen::VectorXd::Ones(z_t.rows()),ech.transpose()));
  for(int h=0; h<H;h++){
    Levier=levierVolatilityMat(rt_sim,Levier,Nl,para,Model_type);
    IntegerVector idx=wrap(MC_sim.col(h).array()-1);
    temp1=sig[idx];
    temp2=wrap(z_t.col(h));
    temp=sqrt(Levier(_,Levier.ncol()-1)*temp1)*temp2;
    rt_sim=cbind(rt_sim,temp);
    temp2=wrap(n_t.col(h));
    rvt(_,h)=(Levier(_,Levier.ncol()-1)*temp1)*temp2;
  }
  NumericMatrix rt2 = Pow(rt_sim,2);
  if(Model_type=="logN"){
    rvt=Pow(rt2,para[1])*exp(para[0]+0.5*para[4]+0.5*pow(para[2],2)-para[3])/sqrt(1-2*para[3]);
  }
  
  Rcpp::List output =
      Rcpp::List::create(
        Rcpp::Named("rt_sim") = rt_sim,
        Rcpp::Named("rt2") = colMeansRcpp(as<NumericMatrix>(rt2)),
        Rcpp::Named("rvt") = colMeansRcpp(as<NumericMatrix>(rvt))
      );
  
  return output;
}

// [[Rcpp::export]]
Rcpp::List f_sim_InvG(const int &H,
                      const Eigen::Map<Eigen::VectorXd> &sig,
                      const Eigen::Map<Eigen::VectorXd> &pi_0,
                      const Eigen::Map<Eigen::MatrixXd> &matP
                      ){
  Eigen::MatrixXd temp = matP;
  NumericVector temp2_r(H);
  NumericVector temp3_r =  wrap((pi_0.transpose())*sig);
  temp2_r[0] = temp3_r[0];
  for(int h=1;h<H;h++){
    temp3_r = wrap((pi_0.transpose() * temp)*sig);
    temp2_r[h] = temp3_r[0];
    temp = temp * matP;
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt2") = temp2_r,
      Rcpp::Named("rvt") = temp2_r
    );
  return output;
}

// [[Rcpp::export]]
Rcpp::List f_sim_logN(const int &H,
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
# K<-2
# Nl<-3
# LEVIER<-FALSE
# para<-c(-0.31636,  0.90798, -0.03042,  0.07726,  0.14792,0.8,0.7)
# Vol<-c(0.2,0.5)
# P<-matrix(c(0.2,0.5,0.8,0.5),2,2,F)
# solve(t(diag(2)-P+1),rep(1,2))
# para_tilde<-natWork(para,Vol,P,LEVIER,2)
# workNat(para_tilde,LEVIER,2)
# ech<-c(0.7,0.5,-0.3,-1.75,-0.7,0.5,-0.5,0.8,0.3,-0.7,0.8)
# ech<-cbind(rv=c(0.15,1.5,0.7,0.8,1.9,0.8,0.3,0.54,0.76,0.98,0.01),r=ech)

# para_tilde<-c(-0.7729722,-0.9088886,1.1451953,2.9137813,-3.8600031)

# l<-logLik(ech,para_tilde,LEVIER,K)
# 
# LEVIER<-TRUE


*/

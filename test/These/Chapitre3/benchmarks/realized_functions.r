#-------------------------------------------------------------------------------
#Load functions and variables
#-------------------------------------------------------------------------------

library(Rsolnp)
library(caTools)#runmean
library(actuar) #inverse gamma and other distributions
#library(highfrequency)
library(mhsmm) #allows to simulate MC quickly

#-------------------------------------------------------------------------------
#MEM and HAR models
#-------------------------------------------------------------------------------

real_to_working_MEM <- function(w, ad, aw, am, b, g1, g2, 
                                  Ld, Lw, Lm, o, Vt, dist, init){

  w.tr <- ad.tr <- aw.tr <- am.tr <- b.tr <- g1.tr <- g2.tr <- numeric(0)
  Ld.tr <- Lw.tr <- Lm.tr <- o.tr <- v1.tr <- numeric(0)

  if(Vt$model %in% c('MEM','AMEM')){

    w.tr  <- sqrt(w)         #>0
    ad.tr <- sqrt(1/ad - 1)  #(0,1)
    b.tr  <- sqrt(1/b  - 1)  #(0,1)
    if(Vt$model=='AMEM') g1.tr <- sqrt(g1) #>0

  } else if(Vt$model %in% c('HAR','HAR-J','HAR-RS')){

    if(Vt$logHAR) w.tr  <- w  else w.tr  <- sqrt(w)        #>0
    if(Vt$logHAR) ad.tr <- ad else ad.tr <- sqrt(1/ad - 1) #(0,1)
    if(Vt$logHAR) aw.tr <- aw else aw.tr <- sqrt(1/aw - 1) #(0,1)
    if(Vt$logHAR) am.tr <- am else am.tr <- sqrt(1/am - 1) #(0,1)

    if(Vt$model=='HAR-J') {
      if(Vt$logHAR) g1.tr <- g1 else g1.tr <- tan(g1*pi/2) #[-1,1]
    } else if(Vt$model=='HAR-RS') {
      if(Vt$logHAR) g2.tr <- g2 else g2.tr <- sqrt(g2) #>0
    }
  } else if(Vt$model %in% c('HAR-lev','HAR-J-lev','HAR-RS-lev')){

    if(Vt$logHAR) w.tr  <- w  else w.tr  <- sqrt(w)        #>0
    if(Vt$logHAR) ad.tr <- ad else ad.tr <- sqrt(1/ad - 1) #(0,1)
    if(Vt$logHAR) aw.tr <- aw else aw.tr <- sqrt(1/aw - 1) #(0,1)
    if(Vt$logHAR) am.tr <- am else am.tr <- sqrt(1/am - 1) #(0,1)

    if(Vt$logHAR) Ld.tr <- Ld else Ld.tr <- sqrt(-Ld) #<0
    if(Vt$logHAR) Lw.tr <- Lw else Lw.tr <- sqrt(-Lw) #<0
    if(Vt$logHAR) Lm.tr <- Lm else Lm.tr <- sqrt(-Lm) #<0

    if(Vt$model=='HAR-J-lev') {
      if(Vt$logHAR) g1.tr <- g1 else g1.tr <- tan(g1*pi/2) #[-1,1]
    } else if(Vt$model=='HAR-RS-lev') {
      if(Vt$logHAR) g2.tr <- g2 else g2.tr <- sqrt(g2) #>0
    }
  }

  #parameter of error distribution
  #must be at least positive
  o <- o - 1e-5 #we bound it away from zero
  if(dist %in% c('igamma', 'iweibull', 'loglogistic', 'paralogistic',
                 'iparalogistic', 'pareto')){
    o <- o - 1
  }

  o.tr <- sqrt(o)

  if(Vt$model %in% c('MEM','AMEM')){
    if(init$method=='estim') v1.tr <- sqrt(init$v1)
  }

  all_vars <- c(w.tr, ad.tr, aw.tr, am.tr, b.tr, g1.tr, g2.tr, 
                  Ld.tr, Lw.tr, Lm.tr, o.tr, v1.tr)

  return(as.vector(all_vars))
}


working_to_real_MEM <- function(all_vars, Vt, dist, init){

  if(Vt$model %in% c('MEM','AMEM')){
    aw <- am <- g2 <- Ld <- Lw <- Lm <- 0

    w  <- all_vars[1]^2        #>0
    ad <- 1/(1+all_vars[2]^2)  #(0,1)
    b  <- 1/(1+all_vars[3]^2)  #(0,1)
    all_vars <- all_vars[-c(1:3)]

    if(Vt$model=='MEM') {
      g1 <- 0
    } else {
      g1 <- all_vars[1]^2 #>0
      all_vars <- all_vars[-1]
    }
  } else if(Vt$model %in% c('HAR','HAR-J','HAR-RS')){
    b <- Ld <- Lw <- Lm <- 0

    if(Vt$logHAR) w  <- all_vars[1] else w  <- all_vars[1]^2       #>0
    if(Vt$logHAR) ad <- all_vars[2] else ad <- 1/(1+all_vars[2]^2) #(0,1)
    if(Vt$logHAR) aw <- all_vars[3] else aw <- 1/(1+all_vars[3]^2) #(0,1)
    if(Vt$logHAR) am <- all_vars[4] else am <- 1/(1+all_vars[4]^2) #(0,1)
    all_vars <- all_vars[-c(1:4)]

    g1 <- g2 <- 0
    if(Vt$model=='HAR-J') {
      if(Vt$logHAR) g1 <- all_vars[1] else g1 <- atan(all_vars[1])*2/pi #[-1,1]
      all_vars <- all_vars[-1]
    } else if(Vt$model=='HAR-RS') {
      if(Vt$logHAR) g2 <- all_vars[1] else g2 <- all_vars[1]^2 #>0
      all_vars <- all_vars[-1]
    }
  } else if(Vt$model %in% c('HAR-lev','HAR-J-lev','HAR-RS-lev')){
    b <- 0

    if(Vt$logHAR) w  <- all_vars[1] else w  <- all_vars[1]^2       #>0
    if(Vt$logHAR) ad <- all_vars[2] else ad <- 1/(1+all_vars[2]^2) #(0,1)
    if(Vt$logHAR) aw <- all_vars[3] else aw <- 1/(1+all_vars[3]^2) #(0,1)
    if(Vt$logHAR) am <- all_vars[4] else am <- 1/(1+all_vars[4]^2) #(0,1)
    all_vars <- all_vars[-c(1:4)]

    g1 <- g2 <- 0
    if(Vt$model=='HAR-J-lev') {
      if(Vt$logHAR) g1 <- all_vars[1] else g1 <- atan(all_vars[1])*2/pi #[-1,1]
      all_vars <- all_vars[-1]
    } else if(Vt$model=='HAR-RS-lev') {
      if(Vt$logHAR) g2 <- all_vars[1] else g2 <- all_vars[1]^2 #>0
      all_vars <- all_vars[-1]
    }
    
    if(Vt$logHAR) Ld <- all_vars[1] else Ld <- -all_vars[1]^2 #<0
    if(Vt$logHAR) Lw <- all_vars[2] else Lw <- -all_vars[2]^2 #<0
    if(Vt$logHAR) Lm <- all_vars[3] else Lm <- -all_vars[3]^2 #<0
    all_vars <- all_vars[-c(1:3)]
  }

  #parameter of error distribution
  #must be at least positive
  o <- all_vars[1]^2 + 1e-5 #we bound it away from zero
  all_vars <- all_vars[-1]

  if(dist %in% c('igamma', 'iweibull', 'loglogistic', 'paralogistic',
                 'iparalogistic', 'pareto')){
    o <- 1 + o
  }

  if(Vt$model %in% c('MEM','AMEM')){
    if(init$method=='estim') {
      v1 <- all_vars[1]^2
    } else if(init$method=='unc'){
      v1 <- w / (1-ad-g1-b)
    } else if(init$method=='mean'){
      v1 <- init$v1
    }
  } else v1 <- NA

  return(list(w=w, ad=ad, aw=aw, am=am, b=b, g1=g1, g2=g2, 
                Ld=Ld, Lw=Lw, Lm=Lm, o=o, v1=v1))
}


MEM <- function(all_vars,
                Vt=list(model = c('MEM', 'AMEM', 'HAR', 'HAR-J', 'HAR-RS',
                                  'HAR-lev', 'HAR-J-lev', 'HAR-RS-lev'),
                           RV = c('RVol', 'RVar'),
                       logHAR = c(FALSE, TRUE),
                       aggHAR = c('-', 'RVar', 'AR'),
                         mult = 100^2 #for realized variance
                         ),
                dist = c('gamma', 'igamma', 'logN', 'weibull', 'iweibull',
                         'loglogistic', 'paralogistic', 'iparalogistic',
                         'pareto'),
                init = list(method=c('mean','estim','unc'), v1=NA),#not used for HAR models
                Data){

      params <- working_to_real_MEM(all_vars, Vt, dist, init)
      w <- params$w; ad <- params$ad; aw <- params$aw; am <- params$am;
      b <- params$b; g1 <- params$g1; g2 <- params$g2; 
      Ld <- params$Ld; Lw <- params$Lw; Lm <- params$Lm; 
      o  <- params$o; v1 <- params$v1

      mult   <- Vt$mult
      logHAR <- Vt$logHAR

      R  <- Data$RV$r         #returns
      RV <- Data$RV$rv * mult #realized variance
      n  <- length(RV)

      if(Vt$RV=='RVol') RV <- sqrt(RV)

      if(Vt$model %in% c('MEM','AMEM')){

        #v1 <- mean(RV)
        e  <- w + ( ad + g1*(R[-n]<0) ) * RV[-n]
        vt <- c( v1, filter(e, b, "r", init=v1) )

      } else if(Vt$model %in% c('HAR', 'HAR-lev', 'HAR-J', 'HAR-J-lev', 'HAR-RS', 'HAR-RS-lev')){

        if(Vt$RV=='RVol'){
          if(Vt$aggHAR=='RVar'){
            reg <- as.data.frame(Data$reg$RVar)
            RVd <- sqrt( reg$rvd * mult )
            RVw <- sqrt( reg$rvw * mult )
            RVm <- sqrt( reg$rvm * mult )
            Jd  <- sqrt( reg$jd  * mult )
            RSd <- sqrt( reg$rsd * mult )

            if(logHAR==TRUE){
              RVd <- log(RVd)
              RVw <- log(RVw)
              RVm <- log(RVm)
              Jd  <- log( 1 + sqrt(reg$jd * mult) )
              if(min(RSd)>0) RSd <- log(RSd)
            }
          } else if(Vt$aggHAR=='AR'){
            if(logHAR==TRUE){
              reg <- as.data.frame(Data$reg$logRVol)
              RVd <- reg$rvd + sqrt( mult )
              RVw <- reg$rvw + sqrt( mult )
              RVm <- reg$rvm + sqrt( mult )
              Jd  <- log( 1 + exp(reg$jd)*sqrt(mult) )
              RSd <- reg$rsd + sqrt( mult )
            } else if(logHAR==FALSE){
              reg <- as.data.frame(Data$reg$RVol)
              RVd <- reg$rvd * sqrt( mult )
              RVw <- reg$rvw * sqrt( mult )
              RVm <- reg$rvm * sqrt( mult )
              Jd  <- reg$jd  * sqrt( mult )
              RSd <- reg$rsd * sqrt( mult )
            }
          }
        } else if(Vt$RV=='RVar'){

            if(logHAR==TRUE & Vt$aggHAR=='AR'){
              reg <- as.data.frame(Data$reg$logRVar)
              RVd <- reg$rvd + log( mult )
              RVw <- reg$rvw + log( mult )
              RVm <- reg$rvm + log( mult )
              Jd  <- log( 1 + exp(reg$jd)*mult )
              RSd <- reg$rsd + log( mult )
            } else {#without log, Vt$aggHAR='RVar' and Vt$aggHAR='AR' are equivalent
              reg <- as.data.frame(Data$reg$RVar)
              RVd <- reg$rvd * mult
              RVw <- reg$rvw * mult
              RVm <- reg$rvm * mult
              Jd  <- (1 + reg$jd * mult)
              RSd <- reg$rsd * mult

              if(logHAR==TRUE){#Vt$aggHAR='RVar'
                RVd <- log(RVd)
                RVw <- log(RVw)
                RVm <- log(RVm)
                Jd  <- log(Jd)
                if(min(RSd)>0) RSd <- log(RSd)
              }
            }
        }

        #regressors for leverage effect
        reg <- as.data.frame(Data$reg$lev)
        regLd <- reg$rd * 100
        regLw <- reg$rw * 100
        regLm <- reg$rm * 100

        vt <- w + ad*RVd + aw*RVw + am*RVm + g1*Jd + g2*RSd + 
                Ld*regLd + Lw*regLw + Lm*regLm

        if(logHAR==TRUE) vt <- exp(vt)
      }

      if(dist=='gamma'){
        #shape, scale > 0 -> o > 0
        loglik <- sum( log( dgamma(RV/vt, shape=o, scale=1/o) / vt ) )
      } else if(dist=='igamma'){
        #shape, scale > 0, but mean exists only for shape > 1 -> o > 1
        #https://cran.r-project.org/web/packages/actuar/vignettes/lossdist.pdf
        #f(x) = (theta/x)^shape * e^(-theta/x) / ( x * gamma(shape )
        #mean = theta / (shape - 1) for shape > 1
        #moments: minvgamma(order=1, shape=o, scale=o-1)
        loglik <- sum( log( dinvgamma(RV/vt, shape=o, scale=o-1) / vt ) )
      } else if(dist=='logN'){
        #sdlog > 0 -> o > 0
        loglik <- sum( log( dlnorm(RV/vt, meanlog=-o^2/2, sdlog=o) / vt ) )
      } else if(dist=='weibull'){
        #shape (a), scale (b) > 0 -> o > 0
        #f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a), x > 0
        #F(x) = 1 - exp(- (x/b)^a),  x > 0
        #E(X) = b G(1 + 1/a)
        #Var(X) = b^2 * (G(1 + 2/a) - (G(1 + 1/a))^2)
        loglik <- sum( log( dweibull(RV/vt, shape=o, scale=1/gamma(1+1/o)) / vt ) )
      } else if(dist=='iweibull'){
        #The Inverse Weibull distribution with parameters shape = a and scale = s
        #has density: f(x) = a (s/x)^a exp(-(s/x)^a)/x, for x > 0, a > 0 and s > 0.
        #The special case shape == 1 is an Inverse Exponential distribution.
        #E(X^k) = b^k G(1 - k/a)
        #minvweibull(1, shape=2, scale = 1/gamma(1-1/2))
        #we must have o > 1
        loglik <- sum( log( dinvweibull(RV/vt, shape=o, scale=1/gamma(1-1/o)) / vt ) )
      } else if(dist=='loglogistic'){
        #The Loglogistic distribution with parameters shape = a and scale = s
        #has density: f(x) = a (x/s)^a / (x [1 + (x/s)^a]^2), for x > 0, a > 0 and b > 0.
        #E(X) = a * pi/s / sin(pi/s), if a > 1 -> o > 1
        #mllogis(1, shape=2, scale = sin(pi/2)/(pi/2))
        loglik <- sum( log( dllogis(RV/vt, shape=o, scale=sin(pi/o)/(pi/o)) / vt ) )
      } else if(dist=='paralogistic'){
        #The Paralogistic distribution with parameters shape = a and scale = s
        #has density: f(x) = a^2 (x/s)^a / (x [1 + (x/s)^a]^(a + 1)), for x > 0, a > 0 and b > 0.
        #E(X) = s * gamma(1+1/a)*gamma(a-1/a)/gamma(a) provided a > 1 -> o > 1
        #mparalogis(1, shape=2, scale = gamma(2)/(gamma(1+1/2)*gamma(2-1/2)))
        loglik <- sum( log( dparalogis(RV/vt, shape=o,
                              scale=gamma(o)/(gamma(1+1/o)*gamma(o-1/o)) ) / vt ) )
      } else if(dist=='iparalogistic'){
        #The Inverse Paralogistic distribution with parameters shape = a and scale = s
        #has density: f(x) = a^2 (x/s)^(a^2)/(x [1 + (x/s)^a]^(a + 1)), for x > 0, a > 0 and b > 0.
        #E(X) = s * gamma(1+1/a)*gamma(1-1/a)/gamma(a) provided a > 1 -> o > 1
        #minvparalogis(1, shape=2, scale = gamma(2)/(gamma(2+1/2)*gamma(1-1/2)))
        loglik <- sum( log( dinvparalogis(RV/vt, shape=o,
                              scale=gamma(o)/(gamma(o+1/o)*gamma(1-1/o)) ) / vt ) )
      } else if(dist=='pareto'){
        #also known as lomax
        #The Pareto distribution with parameters shape = a and scale = s has density:
        #f(x) = a s^a / (x + s)^(a + 1), for x > 0, a > 0 and s > 0.
        #E(X) = s / (a-1) provided a > 1 -> o > 1
        #mpareto(1, shape=2, scale=2-1)
        loglik <- sum( log( dpareto(RV/vt, shape=o, scale=o-1) / vt ) )
      } #for inverse lomax, the mean does not exist


      loglik <- -loglik
      attr(loglik, "Vt") <- vt
      return(loglik)

}


#-------------------------------------------------------------------------------
#Regime-Switching MEM based on Haas et al.
#-------------------------------------------------------------------------------

real_to_working_RSMEM <- function(a0, a1, b, g1, o, P, v1,
                                    n_regime, equal_by_reg, Vt, dist, init){

  a0.tr <- a1.tr <- b.tr <- g1.tr <- o.tr <- P.tr <- v1.tr <- numeric(0)
 
  len_par <- (!equal_by_reg)*(n_regime-1)+1

  a0.tr <- sqrt(a0[1:len_par[1]])         #>0
  a1.tr <- sqrt(a1[1:len_par[2]])         #>0
  #a1.tr <- sqrt(1/a1[1:len_par[2]] - 1)  #(0,1)
  b.tr  <- sqrt(1/b[1:len_par[3]]  - 1)   #(0,1)
  if(Vt$model=='AMEM') g1.tr <- sqrt(g1[1:len_par[4]]) #>0

  #parameter of error distribution
  #must be at least positive
  o <- o[1:len_par[5]] - 1e-5 #we bound it away from zero
  if(dist %in% c('igamma', 'iweibull', 'loglogistic', 'paralogistic',
                 'iparalogistic', 'pareto')){
    o <- o - 1
  }

  o.tr <- sqrt(o)

  P.tr <- sqrt(P / P[,n_regime])[,-n_regime]
  P.tr <- as.vector(t(P.tr))

  if(init$v$method=='estim') v1.tr <- sqrt(init$v$value)

  all_vars <- c(a0.tr, a1.tr, b.tr, g1.tr, o.tr, P.tr, v1.tr)

  return(as.vector(all_vars))
}


working_to_real_RSMEM <- function(all_vars, n_regime, equal_by_reg, Vt, dist, init){

  len_par <- (!equal_by_reg)*(n_regime-1)+1
  #there are 5 parameters that can change regimes
  #(a0, a1, b, g1, o)
  #equal_by_reg specifies whether each one should be held
  #constant in all regimes

  a0 <- all_vars[1:len_par[1]]^2 #>0
  all_vars <- all_vars[-c(1:len_par[1])]
  
  a1 <- all_vars[1:len_par[2]]^2 #>0
  #a1 <- 1/(1+all_vars[1:len_par[2]]^2) #(0,1)
  all_vars <- all_vars[-c(1:len_par[2])]  
  
  b  <- 1/(1+all_vars[1:len_par[3]]^2) #(0,1)
  all_vars <- all_vars[-c(1:len_par[3])]

  if(Vt$model=='MEM') {
    g1 <- 0
  } else {
    g1 <- all_vars[1:len_par[4]]^2 #>0
    all_vars <- all_vars[-c(1:len_par[4])]
  }

  #parameters of error distribution (must be at least positive)
  o <- all_vars[1:len_par[5]]^2 + 1e-5 #we bound it away from zero
  all_vars <- all_vars[-c(1:len_par[5])]

  #for some distributions, the parameter must be greater than 1
  if(dist %in% c('igamma', 'iweibull', 'loglogistic', 'paralogistic',
                 'iparalogistic', 'pareto')){
    o <- 1 + o
  }

  #parameters of the transition matrix
  P.tr <- all_vars[1:((n_regime-1)*n_regime)]
  all_vars <- all_vars[-c(1:((n_regime-1)*n_regime))]
  P <- cbind(matrix(P.tr, ncol = n_regime-1, byrow=TRUE), rep(1,n_regime))^2
  P <- P / apply(P, 1, sum)  
  
  #stationary probabilities
  if (n_regime == 2) {
      pi_0 <- c(P[2,1] / (P[1,2] + P[2,1]),
                  P[1,2] / (P[1,2] + P[2,1]))
  } else {    
    pi_0 <- try( solve(t(diag(n_regime) - P + 1), rep(1,n_regime)), silent=TRUE )
    #if error
    if(length(attr(pi_0, "class"))!=0) pi_0 <- rep(1/n_regime, n_regime)
  }
  
  #pi_init represents regime probabilities at t=1
  if(init$pi$method=='stationary') {
     #pi_0 gives regime probs at t=1 (and t=0 implicitly)
     pi_init <- pi_0
   } else if(init$pi$method=='ident'){
     #init$pi$value gives regime at t=0
     pi_init <- P[init$pi$value,]
   } else if(init$pi$method=='fixed'){
     #init$pi$value gives regime probs at t=1
     pi_init <- init$pi$value
   }

  #initial value in MEM recursion
  if(init$v$method=='estim') {
    v1 <- all_vars^2 #>0
  } else {
    v1 <- init$v$value
  }

  return(list(a0=a0, a1=a1, b=b, g1=g1, o=o, P=P, pi_init=pi_init, v1=v1))
}
        
RSMEM <- function(all_vars,
                n_regime=2,
                #there are 5 parameters that can change regimes
                #(a0, a1, b, g1, o)
                #equal_by_reg specifies whether each one should be held
                #constant in all regimes
                equal_by_reg=rep(FALSE, 5),
                Vt=list(model = c('MEM', 'AMEM'),
                          pow = 1,
                           RV = c('RVol', 'RVar'),
                         mult = 100^2 #for realized variance
                         ),
                dist = c('gamma', 'igamma', 'logN', 'loglogistic'),                         
                #init specifies how initial values should be treated
                init=list( 'pi'=list(method=c('stationary', 'ident', 'fixed'),
                                     value),
                           #'estim': regime probs at t=1 are estimated, this
                           #is known to give a vector with 0s and one 1
                            #value specifies the initial guess of the starting
                            #regime (not vector of probs)
                           #'stationary': regime probs at t=1 (and implicitly
                            #at t=0) are constrained to correspond to the
                            #stationary distribution of P
                            #value is unusued
                           #'ident': starting regime is assumed to be known
                           #at t=0, this can act as an identification constraint
                            #value gives regime position at t=0 (not vector of probs)
                           #'fixed': regime probs at t=1 are given and fixed
                            #value gives the vector of regime probs
                           'v'=list(method=c('fixed', 'estim'), value) ),
                           #'fixed': value gives variance at t=0, this is so
                           #to preserve consistency with previous code
                           #'estim': variance at t=1 is estimated and value
                           #gives the initial choice for that variance
                 Data){
  #------------------------          
  #load returns and RV data
  #------------------------          
  mult   <- Vt$mult
  R  <- Data$RV$r         #returns
  RV <- Data$RV$rv * mult #realized variance
  n  <- length(RV)
  if(Vt$RV=='RVol') RV <- sqrt(RV)

  #------------------------          
  #load parameters
  #------------------------          
  params <- working_to_real_RSMEM(all_vars, n_regime, equal_by_reg, Vt, dist, init)
  a0 <- params$a0; a1 <- params$a1; b <- params$b; g1 <- params$g1; 
  P <- params$P; pi_init <- params$pi_init; v1 <- params$v1; o <- params$o;
  pow <- Vt$pow
  #we match length of parameters with the number of regimes
  a0 <- rep(a0, n_regime-length(a0)+1)
  o  <- rep(o , n_regime-length(o)+1)

  #matrix of one-step ahead conditional means in each regime
  vt_ <- matrix(ncol=n_regime, nrow=n) 

  #------------------------          
  #t = 1
  #------------------------          
  w.pred <- pi_init #pi_init represents regime probabilities at t=1

  if(init$v$method=='estim'){
    vt_pow <- (a0-a0) + v1^pow
    #v1 represents the conditional mean at t=1 if it IS estimated
  } else {
    vt_pow <- a0 + a1*v1^pow + b*v1^pow
    #v1 represents the conditional mean at t=0 if it IS NOT estimated
  }

  vt <- vt_pow^(1/pow)
  vt_[1,] <- vt

  if(dist=="gamma"){
    #shape, scale > 0 -> o > 0
    a_j <- w.pred * dgamma(RV[1]/vt, shape=o, scale=1/o) / vt
  } else if(dist=="igamma"){
    #shape, scale > 0, but mean exists only for shape > 1 -> o > 1
    #https://cran.r-project.org/web/packages/actuar/vignettes/lossdist.pdf
    #f(x) = (theta/x)^shape * e^(-theta/x) / ( x * gamma(shape )
    #mean = theta / (shape - 1) for shape > 1
    #moments: minvgamma(order=1, shape=o, scale=o-1)
    a_j <- w.pred * dinvgamma(RV[1]/vt, shape=o, scale=o-1) / vt
  } else if(dist=='logN'){
    #sdlog > 0 -> o > 0
    a_j <- w.pred * dlnorm(RV[1]/vt, meanlog=-o^2/2, sdlog=o) / vt
  } else if(dist=='loglogistic'){
    #The Loglogistic distribution with parameters shape = a and scale = s
    #has density: f(x) = a (x/s)^a / (x [1 + (x/s)^a]^2), for x > 0, a > 0 and b > 0.
    #E(X) = a * pi/s / sin(pi/s), if a > 1 -> o > 1
    #mllogis(1, shape=2, scale = sin(pi/2)/(pi/2))
    a_j <- w.pred * dllogis(RV[1]/vt, shape=o, scale=sin(pi/o)/(pi/o)) / vt
  } 

  a <- sum(a_j)
  w <- a_j / a
  loglik <- log(a)
  if(is.na(loglik)) return(Inf)

  #------------------------          
  #t = 2, 3, ...
  #------------------------          

  for (i in 2:n){

    if(Vt$model=="MEM"){
      vt_pow <- a0 + a1*RV[i-1]^pow + b*vt_pow
    } else if(Vt$model=="AMEM"){
      vt_pow <- a0 + a1*RV[i-1]^pow + b*vt_pow + g1*(R[i-1]<0)*RV[i-1]^pow
    }
    vt <- vt_pow^(1/pow)
    vt_[i,] <- vt

    w.pred <- as.vector(w %*% P)

    if(dist=="gamma"){
      a_j <- w.pred * dgamma(RV[i]/vt, shape=o, scale=1/o) / vt
    } else if(dist=="igamma"){
      a_j <- w.pred * dinvgamma(RV[i]/vt, shape=o, scale=o-1) / vt
    } else if(dist=='logN'){
      a_j <- w.pred * dlnorm(RV[i]/vt, meanlog=-o^2/2, sdlog=o) / vt
    } else if(dist=='loglogistic'){
      a_j <- w.pred * dllogis(RV[i]/vt, shape=o, scale=sin(pi/o)/(pi/o)) / vt
    } 
      
    a <-  sum(a_j)
    w <- a_j / a
    loglik <- loglik + log(a)
    if(is.na(loglik)) return(Inf)
  }
  
  loglik <- -loglik
  attr(loglik, "Vt")  <- vt_
  attr(loglik, "VT")  <- vt
  attr(loglik, "RVT") <- RV[n]
  attr(loglik, "RT")  <- R[n]
  attr(loglik, "w")   <- w

  
  return(loglik)

}


#-------------------------------------------------------------------------------
# SIMULATION OF MS-AMEM BASED HAAS ET AL.
# unconstrained / two regimes / gamma distribution
#-------------------------------------------------------------------------------

RSMEM_sim <- function(params, dist, l_sim, n_sim, pi_0, vt_0, RV_0, R_0){

  a0 <- params$a0; a1 <- params$a1; b <- params$b; g1 <- params$g1;
  o <- params$o; P <- params$P; 
  a1 <- rep(a1, n_regime-length(a1)+1)
  b  <- rep(b, n_regime-length(b)+1)
  g1 <- rep(g1, n_regime-length(g1)+1)
  o  <- rep(o, n_regime-length(o)+1)

  #----------
  #simulation
  #----------
  RV <- matrix(NA, nrow=l_sim, ncol=n_sim) 
  pi_1 <- as.vector(pi_0 %*% P)
  MC_sim <- matrix(sim.mc(pi_1, P, rep(l_sim, n_sim)), ncol=n_sim)
  
  for(j in 1:n_sim){

    now <- MC_sim[1,j]
    vt  <- a0 + (a1 + g1*(R_0<0))*RV_0 + b*vt_0
    RV[1,j] <- rgamma(1, shape=o, scale=1/o) * vt[now]
      
    for(i in 2:l_sim){
      now <- MC_sim[i,j]
      vt  <- a0 + (a1 + g1*0.5)*RV[i-1,j] + b*vt
      RV[i,j] <- rgamma(1, shape=o, scale=1/o) * vt[now]
    }
  }

  return(t(RV))

}

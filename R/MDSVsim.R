#' @title MDSV Simulation
#' @description Method for simulating from a variety of MDSV model (uniquely or jointly).
#' @param N An integer designing the number of components for the MDSV process
#' @param K An integer designing the number of states of each MDSV process component
#' @param para A vector of parameters use for the MDSV simulation.
#' @param ModelType An integer designing the type of model to be fit. \eqn{0} for univariate log-returns, \eqn{1} for univariate realized variances and \eqn{2} for joint log-return and realized variances.
#' @param LEVIER if \code{TRUE}, estime the MDSV model with leverage.
#' @param n.sim The simulation horizon.
#' @param n.start The burn-in sample.
#' @param m.sim The number of simulations.
#' @param rseed An integer use to initialize the random number generator for the resampling with replacement method (if not supplied take randomly).
#' 
#' @return A list consisting of:
#' \itemize{
#'     \item ModelType : type of model to be fitted.
#'     \item LEVIER : wheter the fit take the leverage effect into account or not.
#'     \item N : number of components for the MDSV process.
#'     \item K : number of states of each MDSV process component.
#'     \item parameters : parameters used for the simulation.
#'     \item n.sim : simulation horizon.
#'     \item n.start : burn-in sample.
#'     \item m.sim : number of simulations.
#'     \item r_t : simulated log-returns (only if ModelType == 0 or ModelType == 2).
#'     \item RV_t : simulated realized variances (only if ModelType == 1 or ModelType == 2).
#' }
#' 
#' @details 
#' n.start is the simulation horizon to be delete (from the start). The leverage effect is taken into 
#' account according to the FHMV model (see Augustyniak et al., 2019). While simulating an
#' univariate realized variances data, this package does not allow to add leverage effect.
#' 
#' @references  
#' Augustyniak, M., Bauwens, L., & Dufays, A. (2019). A new approach to volatility modeling: the factorial hidden Markov volatility model. 
#' \emph{Journal of Business & Economic Statistics}, 37(4), 696-709. \url{https://doi.org/10.1080/07350015.2017.1415910}
#' 
#' @seealso For fitting \code{\link{MDSVfit}}, filtering \code{\link{MDSVfilter}}, bootstrap forecasting \code{\link{MDSVboot}} and rolling estimation and forecast \code{\link{MDSVroll}}.
#' 
#' @examples 
#' \dontrun{
#' # MDSV(N=2,K=3) without leverage on univariate log-returns S&P500
#' N         <- 2      # Number of components
#' K         <- 3      # Number of states
#' ModelType <- 0      # Univariate log-returns
#' para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72)
#' LEVIER    <- FALSE  # No leverage effect
#' n.sim     <- 2000   # length of the simulation
#' n.start   <- 70      # No burn-in
#' m.sim     <- 1      # Number of simulation
#' 
#' # simulation
#' out       <- MDSVsim(K = K, N = N, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
#'                       n.start = n.start, m.sim = m.sim, rseed = NA)
#' 
#' 
#' # MDSV(N=3,K=3) with leverage on joint log-returns and realized variances NASDAQ
#' N         <- 3      # Number of components
#' K         <- 3      # Number of states
#' ModelType <- 2      # Univariate log-returns
#' para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72, 
#'               xi = -0.5, varphi = 0.93, delta1 = 0.93, delta2 = 0.04, shape = 2.10,
#'               l = 0.78, theta = 0.876)
#' LEVIER    <- TRUE  # No leverage effect
#' n.sim     <- 2000   # length of the simulation
#' n.start   <- 500      # No burn-in
#' m.sim     <- 3      # Number of simulation
#' 
#' # simulation
#' out       <- MDSVsim(K = K, N = N, para = para, ModelType = ModelType, LEVIER = LEVIER, n.sim = n.sim, 
#'                       n.start = n.start, m.sim = m.sim, rseed = NA)
#' 
#' }

#' @export
#' @importFrom mhsmm sim.mc 
MDSVsim<-function(N, K, para, ModelType = 0, LEVIER = FALSE, n.sim = 1000, n.start = 0, m.sim = 1, rseed = NA){
  
  if ( (!is.numeric(N)) || (!is.numeric(K)) ) {
    stop("MDSVsim(): input N and K must all be numeric!")
  }else if(!(N%%1==0) || !(K%%1==0)){
    stop("MDSVsim(): input N and K must all be integer!")
  }else if(K<2){
    stop("MDSVfit(): input K must be greater than one!")
  }else if(N<1){
    stop("MDSVfit(): input N must be positive!")
  }
  
  if ( (!is.numeric(n.sim)) || (!is.numeric(n.start)) || (!is.numeric(m.sim))) {
    stop("MDSVsim(): input n.sim, n.start and m.sim must all be numeric!")
  }else if(!(n.sim%%1==0) || !(n.start%%1==0) || !(m.sim%%1==0)){
    stop("MDSVsim(): input n.sim, n.start and m.sim must all be integer!")
  }else if((n.sim<1) || (m.sim<1)){
    stop("MDSVfit(): input n.sim and m.sim must be positive!")
  }
  
  if(n.start<0){
    stop("MDSVfit(): input n.start must be positive or 0!")
  }
  
  if(!is.numeric(ModelType)) {
    stop("MDSVsim(): input ModelType must be numeric!")
  }else if(!(ModelType %in% c(0,1,2))){
    stop("MDSVsim(): input ModelType must be 0, 1 or 2!")
  }
  
  if ((!is.numeric(rseed)) & !is.na(rseed)) {
    print("MDSVsim() WARNING: input rseed must be numeric! rseed set to random")
    rseed <- sample.int(.Machine$integer.max,1)
  }else if(is.numeric(rseed)){ 
    if(!(rseed%%1==0)){
      rseed <- floor(rseed)
      print(paste0("MDSVsim() WARNING: input rseed must be an integer! rseed set to ",rseed))
    }
    set.seed(rseed)
  }
  
  if(!is.logical(LEVIER)) {
    stop("MDSVsim(): input LEVIER must be logical!")
  }
  
  if ((ModelType == 1) & LEVIER) {
    print("MDSVsim() WARNING: can not simulate univariate realized variance with leverage effect! Set to FALSE")
    LEVIER <- FALSE
  }
  
  if(LEVIER) n.start <- n.start+70
  
  if(!is.numeric(para)) {
    stop("MDSVsim(): input para must be numeric!")
  }else if(!is.vector(para)) {
    stop("MDSVsim(): input para must be vector!")
  }else if(((!LEVIER) & (ModelType==0) & !(length(para)==5)) ||
           ((!LEVIER) & (ModelType==1) & !(length(para)==6)) ||
           ((!LEVIER) & (ModelType==2) & !(length(para)==10)) ||
           ((LEVIER) & (ModelType==0) & !(length(para)==7)) ||
           ((LEVIER) & (ModelType==1) & !(length(para)==8)) ||
           ((LEVIER) & (ModelType==2) & !(length(para)==12))){
    stop("MDSVsim(): incorrect input para!")
  }
  
  if((para[1]>1) || (para[1]<0)) {
    stop("MDSVsim(): input para[omega] must be between 0 and 1!")
  }else if((para[2]>1) || (para[2]<0)) {
    stop("MDSVsim(): input para[a] must be between 0 and 1!")
  }else if((para[3]<=1)) {
    stop("MDSVsim(): input para[b] must be greater than 1!")
  }else if((para[4]<=0)) {
    stop("MDSVsim(): input para[sigma] must be greater than 0!")
  }else if((para[5]>1) || (para[5]<0)) {
    stop("MDSVsim(): input para[v0] must be between 0 and 1!")
  }else if((ModelType==1) & (para[6]<=0)) {
    stop("MDSVsim(): input para[shape] must be greater than 0!")
  }else if((ModelType==2) & (para[10]<=0)){
    stop("MDSVsim(): input para[shape] must be greater than 0!")
  }else if(LEVIER){
    if(ModelType==0){
      if(para[6]<=0){
        stop("MDSVsim(): input para[l] must be greater than 0!")
      }else if((para[7]>1) || (para[7]<0)){
        stop("MDSVsim(): input para[theta_l] must be between 0 and 1!")
      }
    }else if(ModelType==1){
      if(para[7]<=0){
        stop("MDSVsim(): input para[l] must be greater than 0!")
      }else if((para[8]>1) || (para[8]<0)){
        stop("MDSVsim(): input para[theta_l] must be between 0 and 1!")
      }
    }else if(ModelType==2){
      if(para[11]<=0){
        stop("MDSVsim(): input para[l] must be greater than 0!")
      }else if((para[12]>1) || (para[12]<0)){
        stop("MDSVsim(): input para[theta_l] must be between 0 and 1!")
      }
    }
  }
  
  vars<-c("omega","a","b","sigma","v0")
  if(ModelType==1) vars <- c(vars,"shape")
  if(ModelType==2) vars <- c(vars,"xi","varphi","delta1","delta2","shape")
  if(LEVIER)       vars <- c(vars,"l","theta")
  names(para)           <- vars
  
  out <- NULL
  
  v               <- volatilityVector(para = para, K = K, N = N)
  MatrixP         <- P(para = para, K = K, N = N)
  stationnaryDist <- probapi(omega = para["omega"], K = K, N = N)
  V_t             <- v[sim.mc(stationnaryDist,MatrixP,n.sim+n.start)]
  
  if(!LEVIER){
    for(sim in 1:m.sim){
      tmp <- NULL
      if(ModelType == 0) {
        r_t         <- rnorm(n.sim+n.start,0,sqrt(V_t))
        r_t         <- r_t[(n.start+1):(n.sim+n.start)]
        tmp         <- list(r_t = r_t)
      }else if(ModelType == 2) {
        r_t         <- rnorm(n.sim+n.start,0,sqrt(V_t))
        e_t         <- r_t/sqrt(V_t)
        RV_t        <- exp(para["xi"] + para["varphi"]*log(V_t) + para["delta1"]*e_t + para["delta2"]*(e_t^2-1) + para["shape"]*rnorm(n.sim+n.start))
        RV_t        <- RV_t[(n.start+1):(n.sim+n.start)]
        r_t         <- r_t[(n.start+1):(n.sim+n.start)]
        tmp         <- list(r_t = r_t, RV_t = RV_t)
      }else if(ModelType == 1) {
        RV_t        <- V_t*rgamma(n.sim+n.start,shape = para["shape"], rate = 1/para["shape"])
        RV_t        <- RV_t[(n.start+1):(n.sim+n.start)]
        tmp         <- list(RV_t = RV_t)
      }
      out             <- c(out,list(tmp))
      names(out)[sim] <- paste("sim",sim,sep = ".")
    }
  }else{
    for(sim in 1:m.sim){
      tmp <- NULL
      if(ModelType == 0) {
        r_t         <- rnorm(70,0,sqrt(V_t[1:70]))
        for(t in 1:(n.sim+n.start)){
          lt        <- levierVolatility(ech = r_t, para = para, Model_type = ModelType)$levier
          r_t       <- c(r_t, rnorm(1,0,sqrt(lt*V_t[70+t])))
        }
        r_t         <- r_t[(n.start+1):(n.sim+n.start)]
        tmp         <- list(r_t = r_t)
      }else if(ModelType == 2) {
        r_t         <- rnorm(70,0,sqrt(V_t[1:70]))
        e_t         <- r_t/sqrt(V_t[1:70])
        for(t in 1:(n.sim+n.start-70)){
          lt        <- levierVolatility(ech = r_t, para = para, Model_type = ModelType)$levier
          r_t       <- c(r_t, rnorm(1,0,sqrt(lt*V_t[70+t])))
          e_t       <- c(e_t,r_t[t]/sqrt(lt*V_t[70+t]))
        }
        RV_t        <- exp(para["xi"] + para["varphi"]*log(V_t) + para["delta1"]*e_t + para["delta2"]*(e_t^2-1) + para["shape"]*rnorm(n.sim+n.start))
        RV_t        <- RV_t[(n.start+1):(n.sim+n.start)]
        r_t         <- r_t[(n.start+1):(n.sim+n.start)]
        tmp         <- list(r_t = r_t, RV_t = RV_t)
      }
      out             <- c(out,list(tmp))
      names(out)[sim] <- paste("sim",sim,sep = ".")
    }
  }

  if(ModelType==0) Model_type <- "Univariate log-return"
  if(ModelType==1) Model_type <- "Univariate realized variances"
  if(ModelType==2) Model_type <- "Joint log-return and realized variances"
  
  out<-c(list(ModelType     = Model_type, 
              LEVIER        = LEVIER, 
              N             = N, 
              K             = K, 
              parameters    = para,
              n.sim         = n.sim, 
              n.start       = n.start, 
              m.sim         = m.sim),
         out)
  
  class(out) <- "MDSVsim"
  
  return(out)
}
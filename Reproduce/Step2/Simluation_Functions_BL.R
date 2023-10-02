------------------------------------------------------------------------------------------------------
  # The "comp.y()" function is used to compute the observed logit sensitivities and logit specificities
  #------------------------------------------------------------------------------------------------------
comp.y <- function(Data){
  # Data is the kby4 DTA data provided as a data frame or matrix specifying the column names "TP", "FP", "TN" and "FN" clearly.
  
  library(mada); library(faraway)
  k <- nrow(Data)
  # A continuity correction of 0.5 will be added to a cell which appears to contain 0, since logit(p) won't be defined otherwise
  
  y <- matrix(c(logit(madad(Data, correction.control = "single")$sens$sens), 
                logit(madad(Data, correction.control = "single")$spec$spec)), nrow=k, ncol=2)
  return(y)
}


#------------------------------------------------------------------------------------------------
# The "comp.Psi()" function is used to compute the within-study variances of logit Se and logit Sp
# using the delta-method
#------------------------------------------------------------------------------------------------
comp.Psi <- function(Data){
  # "Data" is the kby4 data frame containing the number of "TP", "TN", "FP" and "FN".
  library(mada); library(faraway)
  k <- nrow(Data); p <- 2
  correction <- 0.5
  Y <- comp.y(Data=Data)
  
  Psi <- list()
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + correction
      Data$FP[i] = Data$FP[i] + correction
      Data$FN[i] = Data$FN[i] + correction
      Data$TN[i] = Data$TN[i] + correction
    }
    n1 <- Data$TP+Data$FN; n2 <- Data$FP+Data$TN
    Psi[[i]] <- matrix(c((1/(n1[i]*(ilogit(Y[i,1])*(1-ilogit(Y[i,1]))))), 0, 0, (1/(n2[i]*(ilogit(Y[i,2])*(1-ilogit(Y[i,2])))))), 2,2)
  }
  
  return(Psi)
}



#--------------------------------------------------------------------------------------------
# The "fit_bnn()" function is used to fit the BNN model of Reitsma et al. (2005)
#--------------------------------------------------------------------------------------------
fit_bnn <- function(Data, method=c("reml", "ml")){
  # "Data" is a list of data frames consisting of the 2X2 DTA table from several studies.
  # The column names of "Data" must include the frequencies "TP", "FN", "FP" and "TN".
  
  # Load the "metafor" package to compute the study-specific logit-transformed Se and Sp
  library(metafor)
  
  # Load the "mvmeta" package to fit the BNN model
  require(mvmeta)
  
  # Specify the number of studies
  k <- nrow(Data)
  
  # Compute the observed logit-transformed Se and Sp and the within-study covariance matrices
  Y <- comp.y(Data = Data)
  Psi <- comp.Psi(Data = Data)
  
  # Since the logit transformation requires a continuity correction, add the typical continuity
  # correction of 0.5 only to the 2x2 tables that contains a cell with 0 count.
  correction <- 0.5
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + correction
      Data$FP[i] = Data$FP[i] + correction
      Data$FN[i] = Data$FN[i] + correction
      Data$TN[i] = Data$TN[i] + correction
    }
  }
  
  # Store the frequency of correct and incorrect test results
  TP <- Data$TP; FP <- Data$FP; FN <- Data$FN; TN <- Data$TN
  
  # Store the number of diseased and non-diseased subjects and total number of patients
  n1 <- TP + FN; n2 <- FP + TN; npat <- n1 + n2
  
  # Store the within-study covariance vector (a zero vector) for each study
  cov.sesp.logit <- rep(0, k)
  
  # To store the logit transformed Se/(1-Sp);
  trsesp.logit <- matrix(0, nrow = k, ncol = 2)
  trsesp.logit <- cbind(tsens = escalc(measure="PLO", xi=TP, ni=n1)[,1],
                        tfpr = escalc(measure="PLO", xi=TN, ni=n2)[,1])
  
  # To store the within-study variances of logit Se and logit (1-Sp);
  var.se.logit = var.sp.logit <- rep(0, k)
  var.se.logit <- escalc(measure="PLO", xi=TP, ni=n1)[,2]
  var.sp.logit <- escalc(measure="PLO", xi=TN, ni=n2)[,2]
  
  # To store the data frame that is used as an input in the "mvmeta()" function;
  Data.mvmeta.logit <- matrix(0, nrow = k, ncol=6)
  Data.mvmeta.logit <- data.frame(pat.num=npat, trse=trsesp.logit[,1], trsp=trsesp.logit[,2],
                                  varse=var.se.logit, covsefpr=cov.sesp.logit, varsp=var.sp.logit)
  
  # Finally, fit the BNN LMM using the logit transformation
  fit.logit <- list()
  mvmeta.logit <- mvmeta(cbind(trse, trsp)~1, S=Data.mvmeta.logit[, 4:6], method = method, data = Data.mvmeta.logit, control=list(maxiter=1000))
  
  # Backtransformed results
  require(faraway)
  summary.logit <- summary(mvmeta.logit)$coefficients
  SeSp <- c(ilogit(summary.logit[1,1]), ilogit(summary.logit[2,1]))
  SeCI <- as.numeric(ilogit(summary.logit[1,5:6]))
  SpCI <- as.numeric(ilogit(summary.logit[2,5:6]))
  
  loglik.bnn <- mvmeta.logit$logLik
  
  fit.logit$mvmeta.logit <- summary(mvmeta.logit)
  fit.logit$SeSp <- SeSp
  fit.logit$SeSp.logit <- logit(SeSp)
  fit.logit$Psi <- mvmeta.logit$Psi
  fit.logit$vcov <- mvmeta.logit$vcov
  fit.logit$CIs <- c("Se.lb"=SeCI[1],"Se.ub"=SeCI[2], "Sp.lb"=SpCI[1],"Sp.ub"=SpCI[2])
  fit.logit$Stats <- c(logLik=loglik.bnn, AIC=-2*loglik.bnn + 2*5, BIC=-2*loglik.bnn + 5*log(k*2))
  
  return(fit.logit)
  
}


#--------------------------------------------------------------------------------------------
# The "fit_bbn()" function is used to fit the BBN model of Chu and Cole (2006)
#--------------------------------------------------------------------------------------------
fit_bbn <- function(Data, Study){
  # "Data" is a data frame consisting of the 2X2 DTA table from several studies.
  # "index" represents the location of the DTA 'Data' in a list of such datasets.
  # The column names of "Data" must include the frequencies "TP", "FN", "FP" and "TN".
  # "study" represents "study names or identities. Vector of characters, need to be speficied
  #                     either directly or by referring to a variable in data frame."
  if("metafor" %in% (.packages())) detach(package:metafor)
  require(Metatron); require(faraway)
  
  k <- nrow(Data)
  
  # Compute the observed logit-transformed Se and Sp and the within-study covariance matrices
  Y <- comp.y(Data = Data)
  Psi <- comp.Psi(Data = Data)
  
  # Now fit the BBN RE model of Chu et al. (2006)
  metatron.bbn <- fit.bivar(TP=TP, FN=FN, TN=TN, FP=FP, study=Study, data=Data)
  
  # Extract the parameter estimates and their covariance matrix or CIs when possible
  loglik.bbn <- metatron.bbn$logLik[1]
  
  fit.bbn <- list()
  
  fit.bbn$metatron.bbn <- metatron.bbn
  fit.bbn$SeSp <- ilogit(as.numeric(metatron.bbn$fixef))
  fit.bbn$SeSp.logit <- as.numeric(metatron.bbn$fixef)
  fit.bbn$Psi <- metatron.bbn$rancoef[1:2,1:2]
  fit.bbn$vcov <- metatron.bbn$vcov
  fit.bbn$CIs <- c("Se.lb"=ilogit(fit.bbn$SeSp.logit[1]-qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[1,1])), "Se.ub"=ilogit(fit.bbn$SeSp.logit[1]+qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[1,1])),
                   "Sp.lb"=ilogit(fit.bbn$SeSp.logit[2]-qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[2,2])), "Sp.ub"=ilogit(fit.bbn$SeSp.logit[2]+qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[2,2])))
  fit.bbn$Stats <- c(logLik=loglik.bbn, AIC= -2*loglik.bbn + 2*5, BIC = -2*loglik.bbn + 5*log(k*2))
  
  return(fit.bbn)
  
}




#----------------------------------------------------------------------------------
# The MCQ.mu() function returns the first term of the full MCQ() function. We
# maximize this function to get the MLEs of mu.
#----------------------------------------------------------------------------------
MCQ.mu <- function(Y, MCb, mu, Psi){
  # Y is a kx2 matrix of logit Se and logit Sp of the k independent studies
  # MCb is a k-dimensional list of Rx2 matrices of random-effects, where R is the
  # MC sample size which varies from iteration to iteration.
  # mu is a 2x1 mean vector. 
  # Psi is a k-dimensional list of 2x2 diagonal matrices of within-study covariance
  # matrix assumed to be known but estimated from the data using the delta method.
  
  k <- nrow(Y)
  R <- nrow(MCb[[k]])
  out1 <- rep(0, k)
  out2 <- list()
  
  for(i in 1:k){
    out2[[i]] <- vector(mode = "numeric", length = R)
    for(r in 1:R){
      out2[[i]][r] <- -0.5*log(2*pi)-0.5*log(det(Psi[[i]]))-0.5*(t(Y[i,]-mu-MCb[[i]][r,])%*%solve(Psi[[i]])%*%(Y[i,]-mu-MCb[[i]][r,]))
    }
    out1[i] <- mean(out2[[i]])
  }
  return(sum(out1))
}


#--------------------------------------------------------------------------------#
# The vec2full() and vec() functions are taken as is from the OpenMx's R package #
# functions vech2full() and vech() functions, respectively.                      #
#--------------------------------------------------------------------------------#
vec2full <- function (x) {
  if (is.matrix(x)) {
    if (nrow(x) > 1 && ncol(x) > 1) {
      stop("Input to the full vech2full must be a (1 x n) or (n x 1) matrix.")
    }
    dimension <- max(dim(x))
  }
  else if (is.vector(x)) {
    dimension <- length(x)
  }
  else {
    stop("Input to the function vech2full must be either a matrix or a vector.")
  }
  k <- sqrt(2 * dimension + 0.25) - 0.5
  ret <- matrix(0, nrow = k, ncol = k)
  if (nrow(ret) != k) {
    stop("Incorrect number of elements in vector to construct a matrix from a half-vectorization.")
  }
  ret[lower.tri(ret, diag = TRUE)] <- as.vector(x)
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
  return(ret)
}


vect <- function (x) {
  return(x[lower.tri(x, diag = TRUE)])
}


diag2vech <- function (x) {
  return(as.matrix(diag(as.matrix(x))))
}


#----------------------------------------------------------------------------------
# The MCQ.sigma() function returns the second term of the full MCQ() function. We
# maximize this function to get the MLEs of Sigma.
#----------------------------------------------------------------------------------
MCQ.sigma <- function(MCb, sigma){
  # MCb is a k-dimensional list of Rx2 matrices of random-effects, where R is the
  # MC sample size which varies from iteration to iteration.
  # sigma is a 3x1 vector containing the unique elements of the between-study 
  # covariance matrix--sigma11, sigma12 and sigma22, respectively. 
  
  k <- length(MCb)
  R <- nrow(MCb[[k]])
  Sigma <- vec2full(sigma) #equivalent to matrix(c(sigma[1:2], sigma[2:3]),2,2)
  out1 <- rep(0, k)
  out2 <- list()
  
  library(matrixcalc)
  
  if(is.positive.definite(Sigma)==FALSE){# To ensure that the estimated covariace matrix is positive definite
    # Positive definite Sigma is required to simulate samples from the bivariate Laplace distribution (see rmvl() function)
    return(-Inf)
  }else{
    
    for(i in 1:k){
      out2[[i]] <- vector(mode = "numeric", length = R)
      for(r in 1:R){
        out2[[i]][r] <- -0.5*log(2^(3/2)*pi)-0.5*log(det(Sigma))-0.25*log(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])-sqrt(2)*sqrt(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])
      }
      out1[i] <- mean(out2[[i]])
    }
    return(sum(out1))
  }
}


#-----------------------------------------------------------------------------------
# The MH.mcmc() function returns a Metropolis-Hastings Markov Chain Monte Carlo 
# random samples from the conditional distribution of the random-effects(b) given
# the response (y) while using the marginal distribution of the random-effects
# (Bivariate Laplace) as a candidate density.
#-----------------------------------------------------------------------------------
MH.mcmc <- function(Y, mu, Sigma, Psi, R){
  # Y is a kx2 matrix of response vector of the k independent studies included in the MA
  # mu is a 2x1 mean vector
  # Sigma is a 2x2 between-study covariance matrix
  # Psi is a k-dimensional list of 2x2 diagonal matrices of within-study covariance
  # matrix assumed to be known but estimated from the data using the delta method.
  # R is the MC sample size
  
  library(LaplacesDemon)
  library(mvtnorm)
  k <- length(Psi)
  set.seed(010)
  b0 <- rmvl(1, mu=c(0,0), Sigma = Sigma)
  b.out <- list()
  b.cand <- list()
  probs <- list()
  U <- list()
  
  for(i in 1:k){
    b.out[[i]] <- matrix(0, ncol = 2, nrow = R)
    b.cand[[i]] <- matrix(0, ncol = 2, nrow = R)
    probs[[i]] <- vector(mode = "numeric", length = R)
    U[[i]] <- vector(mode = "numeric", length = R)
    for(r in 1:R){
      # Generate the candidate densities
      set.seed(r*i^2)
      b.cand[[i]] <- matrix(rmvl(R, mu=c(0,0), Sigma=Sigma), ncol = 2)
      U[[i]] <- runif(R)
      
      if(r==1){ # the current value of b
        b.curr <- b0[1,]
      }
      else{
        b.curr <- b.out[[i]][r-1,]  
      }
      V.curr <- b.cand[[i]][r,] 
      
      probs[[i]][r] <- min(1, ifelse(dmvnorm(Y[i,], mean = mu+V.curr, sigma = Psi[[i]])==0, 1e-323, dmvnorm(Y[i,], mean = mu+V.curr, sigma = Psi[[i]]))/
                             ifelse(dmvnorm(Y[i,], mean = mu+b.curr, sigma = Psi[[i]])==0, 1e-323, dmvnorm(Y[i,], mean = mu+b.curr, sigma = Psi[[i]])))
      
      if(U[[i]][r]<=probs[[i]][r]){
        b.out[[i]][r,] <- V.curr
      }else{
        b.out[[i]][r,] <- b.curr
      }
    }
  }
  return(b.out)
}


#---------------------------------------------------------------------------------------
# The FIM.Oakes() function approximates the observed Fisher information matrix (FIM) 
# using the method proposed by Oakes(1999).
#---------------------------------------------------------------------------------------
FIM.Oakes <- function(Y, MCb, Psi, mu, Sigma, R){
  # Y is a kx2 matrix of response vector of the k independent studies included in the MA
  # MCb is an k-dimensional list of Rx2 matrices of random-effects, where R is the
  # MC sample size which varies from iteration to iteration.
  # Psi is an k-dimesional list of 2x2 diagonal within-study covariance matrix assumed to
  # be known.
  # mu is the 2x1 vector of means
  # Sigma is the 2x2 between-study covariance matrix
  # R is the Monte Carlo sample size
  library(matrixcalc)
  k <- nrow(Y)
  p <- length(mu) + length(vect(Sigma))
  H.out <- matrix(0, nrow = p, ncol = p) 
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigma11 <- Sigma[1,1]
  sigma12 <- Sigma[1,2]
  sigma22 <- Sigma[2,2]
  dmu12 <- 0
  dmu22 <- 0
  dmu1mu2 <- 0
  dsigma112 <- 0
  dsigma222 <- 0
  dsigma122 <- 0
  dsigma11sigma22 <- 0
  dsigma11sigma12 <- 0
  dsigma22sigma12 <- 0
  vmu1 <- list() 
  vmu2 <- list()
  v1sigma11 <- list()
  v2sigma11 <- list()
  v1sigma22 <- list()
  v2sigma22 <- list()
  v1sigma12 <- list()
  v2sigma12 <- list()
  vsigma11sigma22 <- list()
  vsigma11sigma12 <- list()
  vsigma22sigma12 <- list()
  
  for(i in 1:k){
    vmu1[[i]] <- vector(mode= "numeric", length = R) 
    vmu2[[i]] <- vector(mode= "numeric", length = R)
    v1sigma11[[i]] <- vector(mode= "numeric", length = R)
    v2sigma11[[i]] <- vector(mode= "numeric", length = R)
    v1sigma22[[i]] <- vector(mode= "numeric", length = R)
    v2sigma22[[i]] <- vector(mode= "numeric", length = R)
    v1sigma12[[i]] <- vector(mode= "numeric", length = R)
    v2sigma12[[i]] <- vector(mode= "numeric", length = R)
    vsigma11sigma22[[i]] <- vector(mode= "numeric", length = R)
    vsigma11sigma12[[i]] <- vector(mode= "numeric", length = R)
    vsigma22sigma12[[i]] <- vector(mode= "numeric", length = R)
    
    for(r in 1:R){
      vmu1[[i]][r] <- -1*(-Y[i,1]/Psi[[i]][1,1]+mu1/Psi[[i]][1,1]+MCb[[i]][r,][1]/Psi[[i]][1,1])
      vmu2[[i]][r] <- -1*(-Y[i,2]/Psi[[i]][2,2]+mu2/Psi[[i]][2,2]+MCb[[i]][r,][2]/Psi[[i]][2,2])
      
      
      v1sigma11[[i]][r] <- 0.5*sigma22^2/det(Sigma)^2 + 0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(2*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][2]^2-sigma22^2*MCb[[i]][r,][1]^2)^2/det(Sigma)^4 +
        0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma22^2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma22*sigma12^2*MCb[[i]][r,][2]^2-0.5*sigma22^3*MCb[[i]][r,][1]^2)^2/det(Sigma)^3 +
        (1/sqrt(2))*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(2*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][2]^2-sigma22^2*MCb[[i]][r,][1]^2)^2/det(Sigma)^4 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(2*sigma22^2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma22*sigma12^2*MCb[[i]][r,][2]^2-sigma22^3*MCb[[i]][r,][1]^2)/det(Sigma)^3 
      
      v2sigma11[[i]][r] <- -0.5*sigma22/det(Sigma)-0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(2*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][2]^2-sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^2 -
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma12*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][2]^2-0.5*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^2
      
      v1sigma22[[i]][r] <- 0.5*sigma11^2/det(Sigma)^2 + 0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(2*sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]^2-sigma11^2*MCb[[i]][r,][2]^2)^2/det(Sigma)^4 +
        0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma11^2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma11*sigma12^2*MCb[[i]][r,][1]^2-0.5*sigma11^3*MCb[[i]][r,][2]^2)/det(Sigma)^3 +
        (1/sqrt(2))*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][1]^2-0.5*sigma11^2*MCb[[i]][r,][2]^2)^2/det(Sigma)^4 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(2*sigma11^2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma11*sigma12^2*MCb[[i]][r,][2]^2-sigma11^3*MCb[[i]][r,][2]^2)/det(Sigma)^3
      
      v2sigma22[[i]][r] <- -0.5*sigma11/det(Sigma)-0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(2*sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]^2-sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^2 -
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][1]^2-0.5*sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^2
      
      v1sigma12[[i]][r] <- (sigma11*sigma22+sigma12^2)/det(Sigma)^2 + 0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])^2/det(Sigma)^4 - 
        0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma22*MCb[[i]][r,][1]^2*+sigma11*MCb[[i]][r,][2]^2-2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^2 -
        2*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma22*sigma12^2*MCb[[i]][r,][1]^2+sigma11*sigma12^2*MCb[[i]][r,][2]^2-sigma11*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^3*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^3 +
        (1/sqrt(2))*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])^2/det(Sigma)^4 - 
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma22*MCb[[i]][r,][1]^2+sigma11*MCb[[i]][r,][2]^2-2*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^2 - 
        sqrt(2)*4*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma22*sigma12^2*MCb[[i]][r,][1]^2+sigma11*sigma12^2*MCb[[i]][r,][2]^2-sigma11*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^3*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^3
      
      v2sigma12[[i]][r] <- sigma12/det(Sigma)-0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^2 -
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])/det(Sigma)^2
      
      vsigma11sigma12[[i]][r] <- -sigma22*sigma12/det(Sigma)^2 - 0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12*MCb[[i]][r,][2]^2)/det(Sigma)^2 -
        2*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma12^2*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^3*MCb[[i]][r,][2]^2-0.5*sigma12*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^3 + 
        0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(2*sigma22*sigma12*MCb[[i]][r,][1]^2+2*sigma11*sigma12*MCb[[i]][r,][2]^2-2*sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-2*sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])*(2*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][2]^2-sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^4 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])*(sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][2]^2-0.5*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^4 -  
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12*MCb[[i]][r,][2]^2)/det(Sigma)^2 - 
        4*sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma22*sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^3*MCb[[i]][r,][2]^2-0.5*sigma22^2*sigma12*MCb[[i]][r,][1]^2)/det(Sigma)^3
      
      vsigma11sigma22[[i]][r] <- -0.5*(det(Sigma)-sigma22*sigma11)/det(Sigma)^2-0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma22*MCb[[i]][r,][1]^2)/det(Sigma)^2 +
        (t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma11*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma11*sigma12^2*MCb[[i]][r,][2]^2-0.5*sigma11*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^3 +
        0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(2*sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]^2-sigma11^2*MCb[[i]][r,][2]^2)*(2*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][2]^2-sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^4 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][1]^2-0.5*sigma11^2*MCb[[i]][r,][2]^2)*(sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][2]^2-0.5*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^4 -
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma22*MCb[[i]][r,][1]^2)/det(Sigma)^2 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(2*sigma11*sigma22*sigma12*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma11*sigma12^2*MCb[[i]][r,][2]^2-sigma11*sigma22^2*MCb[[i]][r,][1]^2)/det(Sigma)^3
      
      vsigma22sigma12[[i]][r] <- -sigma12*sigma11/det(Sigma)^2-0.5*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12*MCb[[i]][r,][1]^2)/det(Sigma)^2 -
        2*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-1)*(sigma12^2*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^3*MCb[[i]][r,][1]^2-0.5*sigma12*sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^3 +
        0.25*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-2)*(2*sigma22*sigma12*MCb[[i]][r,][1]^2+2*sigma11*sigma12*MCb[[i]][r,][2]^2-2*sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-2*sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])*(2*sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]^2-sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^4 +
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-3/2)*(sigma22*sigma12*MCb[[i]][r,][1]^2+sigma11*sigma12*MCb[[i]][r,][2]^2-sigma11*sigma22*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12^2*MCb[[i]][r,][1]*MCb[[i]][r,][2])*(sigma12*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^2*MCb[[i]][r,][1]^2-0.5*sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^4 - 
        sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-sigma12*MCb[[i]][r,][1]^2)/det(Sigma)^2 - 
        4*sqrt(2)*(t(MCb[[i]][r,])%*%solve(Sigma, tol=.Machine$double.eps*1e-50)%*%MCb[[i]][r,])^(-0.5)*(sigma12^2*sigma11*MCb[[i]][r,][1]*MCb[[i]][r,][2]-0.5*sigma12^3*MCb[[i]][r,][1]^2-0.5*sigma12*sigma11^2*MCb[[i]][r,][2]^2)/det(Sigma)^3
      
    }
    dmu12 <- dmu12 - 1/Psi[[i]][1,1] + mean((vmu1[[i]] - mean(vmu1[[i]]))^2)
    dmu22 <- dmu22 - 1/Psi[[i]][2,2] + mean((vmu2[[i]] - mean(vmu2[[i]]))^2)
    dmu1mu2 <- dmu1mu2 + mean((vmu1[[i]] - mean(vmu1[[i]]))*(vmu2[[i]] - mean(vmu2[[i]])))
    dsigma112 <- dsigma112 + mean(v1sigma11[[i]]) + mean((v2sigma11[[i]] - mean(v2sigma11[[i]]))^2)
    dsigma222 <- dsigma222 + mean(v1sigma22[[i]]) + mean((v2sigma22[[i]] - mean(v2sigma22[[i]]))^2)
    dsigma122 <- dsigma122 + mean(v1sigma12[[i]]) + mean((v2sigma12[[i]] - mean(v2sigma12[[i]]))^2)
    dsigma11sigma12 <- dsigma11sigma12 + mean(vsigma11sigma12[[i]])  + mean((v2sigma11[[i]] -mean(v2sigma11[[i]]))*(v2sigma12[[i]] - mean(v2sigma12[[i]])))
    dsigma11sigma22 <- dsigma11sigma22 + mean(vsigma11sigma22[[i]])  + mean((v2sigma11[[i]] - mean(v2sigma11[[i]]))*(v2sigma22[[i]] - mean(v2sigma22[[i]])))
    dsigma22sigma12 <- dsigma22sigma12 + mean(vsigma22sigma12[[i]])  + mean((v2sigma22[[i]] - mean(v2sigma22[[i]]))*(v2sigma12[[i]] - mean(v2sigma12[[i]])))
  }
  H.out[1,1:2] <- c(dmu12, dmu1mu2)
  H.out[2,1:2] <- c(dmu1mu2, dmu22)
  H.out[3,3:p] <- c(dsigma112, dsigma11sigma12, dsigma11sigma22)
  H.out[4,3:p] <- c(dsigma11sigma12, dsigma122, dsigma22sigma12)
  H.out[p,3:p] <- c(dsigma11sigma22, dsigma22sigma12, dsigma222)
  
  return(-H.out)
  
}



#-------------------------------------------------------------------------------------------
# The function WA_Cov() function computes the weighted average estimate and its corresponding
# asymptotic covariance matrix.
#-------------------------------------------------------------------------------------------
WA_Cov <- function(Data, sigma){
  # Data is a kby4 data frame containing the frequencies corresponding to TP, TN, FP, FN
  # sigma is a vector of length three comprising the unique between-study covariance matrix,
  # sigma11, sigma12, and sigma22
  
  k <- nrow(Data)
  Y <- comp.y(Data)
  Psi <- comp.Psi(Data)
  Sigma <- vec2full(sigma)
  Out <- list()
  Wi <- list()
  WiYi <- list()
  for(i in 1:k){
    Wi[[i]] <- solve(Sigma + Psi[[i]])
    WiYi[[i]] <- Wi[[i]]%*%Y[i,]
  }
  Muhat.W <- as.numeric(solve(Reduce("+", Wi))%*%Reduce("+", WiYi))
  Cov_Muhat.W <- solve(Reduce('+', Wi))
  
  Out$Muhat.W <- Muhat.W
  Out$Cov_Muhat.W <- Cov_Muhat.W
  
  return(Out)
}


#-------------------------------------------------------------------------------------
# The MCloglik.y() function approximates the marginal loglikelihood function using 
# samples generated from the marginal distribution of the random-effects.
#-------------------------------------------------------------------------------------
MCloglik.y <- function(Y, mu, sigma, Psi, R){
  # Y is a kx2 matrix of response vector of the k independent studies included in the MA
  # mu is a 2x1 mean vector
  # sigma is a 3x1 vector containing the unique elements of the between-study 
  # covariance matrix--sigma11, sigma12 and sigma22, respectively. 
  # Psi is a k-dimensional list of 2x2 diagonal matrices of within-study covariance
  # matrix assumed to be known but estimated from the data using the delta method.
  # R is the Monte Carlo sample size
  
  library(LaplacesDemon)
  library(mvtnorm)
  k <- nrow(Y)
  Sigma <- vec2full(sigma)
  out1 <- rep(0, k)
  out2 <- list()
  MCb <- list()
  
  for(i in 1:k){
    set.seed(2*i^2)
    out2[[i]] <- vector(mode = "numeric", length = R)
    MCb[[i]] <- matrix(0, nrow=R, ncol=2)
    MCb[[i]] <- rmvl(R, mu=c(0,0), Sigma = Sigma)
    for(r in 1:R){
      out2[[i]][r] <- dmvnorm(x=Y[i,], mean = mu+MCb[[i]][r,], sigma = Psi[[i]])
    }
    out1[i] <- log(mean(out2[[i]]))
  }
  return(sum(out1))
}


#-------------------------------------------------------------------------------------------
# The Sim_DTAData() function is used to generate DTA datasets assuming the bivariate laplace 
# distribution for the observed study-specific logit(Se) and logit(Sp).
#-------------------------------------------------------------------------------------------
Sim_DTAData <- function(mu, Sigma, n1, n2, k, seed){
  # 'mu' is the vector of true means
  # 'Sigma' is the true between-study covariance matrix
  # 'n1' is the study-specific number of diseased subjects
  # 'n2' is the study-specific number of non-diseased subjects
  # 'k' is the number of studies in the MA
  # 'seed' is a random seed to be used when generating 
  # random variates from a given distribuion function
  
  require(LaplacesDemon)
  require(faraway)
  
  transSeSp = btransSeSp <- matrix(0, nrow = k, ncol = 2)
  TP = FP = TN = FN <- numeric(0)
  Data <- matrix(0, nrow = k, ncol = 5)
  
  set.seed(seed)
  transSeSp <- rmvl(n=k, mu = mu, Sigma = Sigma)
  btransSeSp <- ilogit(transSeSp)
  
  for(i in 1:k){
    set.seed(i+seed)
    TP[i] <- sum(rbinom(n=1, size = n1[i], prob = btransSeSp[i,1]))
    FN[i] <- n1[i] - TP[i]
    TN[i] <- sum(rbinom(n=1, size = n2[i], prob = btransSeSp[i,2]))
    FP[i] <- n2[i] - TN[i]
  }
  Data <- data.frame(Study=as.factor(paste("Study", sep=" ", 1:k)), TP=TP, FN=FN, FP=FP, TN=TN)
  
  return(Data)
}


#----------------------------------------------------------------------------------
# The mcemDTA() function returns parameter estimates of the Bivariate Normal-Laplace
# RE model along with their standard errors using the Monte Carlo EM algorithm. Note
# that the starting values are obtained from the Reitsma et al. (2005) model. The 
# starting MC sample size is defaulted to 10. However, users can specify their own
# R value too.
#----------------------------------------------------------------------------------
mcemDTA <- function(Scen, mu0=c(0,0), sigma0=c(0,0,0), R=10, delta = 0.001, eps = 0.005, maxit=150){
  # Scen is a number denoting the simulation scenario
  # mu0 is a 2x1 vector of starting values denoting the mean vectors.
  # Sigma0 is a 3x1 vector of starting values denoting the three unique values of the
  # between-study covariance matrix.
  # By default, the starting MC sample size is set to be 10, although it would increase
  # as the algorithm progresses.
  # delta and epsilon are parameters used to control convergence of the algorithm
  
  tryCatch({
    
    library(mada)
    library(matrixcalc)
    library(faraway)
    
    # load the data that corresponds to the simulation scenario and also
    # 'outl', which is a vector representing the study numbers that are 
    # identified as outlying.
    Data <- Data.nr[[Scen]]
    out <- outl[[Scen]]
    k <- nrow(Data)
    
    # Set starting values
    if(all(mu0==0)){
      if(all(out==0)){
        Dat <- Data.nr[[Scen]]
      }else{
        Dat <- Data.nr[[Scen]][-out,]
      }
      fit0 <- fit_bnn(Data = Dat, method = "ml")
      mu0 <- fit0$SeSp.logit 
    }
    if(all(sigma0==0)){# We test for the positive definiteness of sigma0 since the rmvl()
      # function shown below requires that.
      if(all(out==0)){
        Dat <- Data.nr[[Scen]]
      }else{
        Dat <- Data.nr[[Scen]][-out,]
      }
      fit0 <- fit_bnn(Data = Dat, method = "ml")
      if(is.positive.definite(fit0$Psi)==TRUE){
        sigma0 <- c(vect(fit0$Psi))
      }else{
        sigma0 <- c(vect(diag(1,2)))
      }
    }
    
    # Standard errors of the mean vectors from the standard model
    mu0.se <- c(sqrt(fit0$vcov[1,1]), sqrt(fit0$vcov[2,2]))
    
    Psi <- comp.Psi(Data = Data)
    Y <- comp.y(Data = Data)
    b <- list()
    EM.out <- matrix(0, nrow=maxit+1, ncol = 13)
    EM.out[1,] <- c(mu0, sigma0, mu0.se, NA, NA, NA, NA, NA, fit0$Stats[1])
    Out <- list()
    j <- 1
    cont <- TRUE
    
    while(cont){
      mu.old <- EM.out[j,1:2]
      sigma.old <- EM.out[j, 3:5]
      muold.se <- EM.out[j, 6:7]
      j <- j+1
      sum1 <- matrix(0, nrow = 2, ncol = 2) 
      sum2 <- c(0, 0)
      
      #set.seed(280918)
      b <- MH.mcmc(Y=Y, mu=mu.old, Sigma = vec2full(sigma.old), Psi = Psi, R=R)
      
      for(i in 1:k){
        sum1 <- sum1 + solve(Psi[[i]])
        for(r in 1:R){
          sum2 <- sum2 + solve(Psi[[i]])%*%(Y[i,]-b[[i]][r,])
        }
      }
      
      mu.new <- as.vector(solve(sum1)%*%sum2)/R
      sigma.new <- optim(par = sigma.old, fn = MCQ.sigma, MCb=b, method = "Nelder-Mead",
                         control=list(fnscale=-1, maxit=20000))$par
      
      EM.out[j,1:7] <- c(mu.new, sigma.new, NA, NA)
      EM.out[j,8] <- abs(mu.new[1]-mu.old[1])/(abs(mu.new[1])+delta) 
      EM.out[j,9] <- abs(mu.new[2]-mu.old[2])/(abs(mu.new[2])+delta)
      EM.out[j,10] <- abs(sigma.new[1]-sigma.old[1])/(abs(sigma.old[1])+delta)
      EM.out[j,11] <- abs(sigma.new[2]-sigma.old[2])/(abs(sigma.old[2])+delta)
      EM.out[j,12] <- abs(sigma.new[3]-sigma.old[3])/(abs(sigma.old[3])+delta)
      EM.out[j,13] <- MCloglik.y(Y=Y, mu=mu.new, sigma = sigma.new, Psi=Psi, R=R)
      
      FIM <- FIM.Oakes(Y=Y, MCb=b, Psi = Psi, mu=mu.new, Sigma = vec2full(sigma.new), R=R)
      
      if(j==2 | j==3){
        cont <- TRUE
      }else{
        cont <- (max(round(EM.out[(j-2),8:12],3))>eps | max(round(EM.out[(j-1),8:12],3))>eps | max(round(EM.out[j,8:12],3))>eps | is.positive.definite(FIM[1:2,1:2])==FALSE) & j<=maxit
      }
      
      if(cont==TRUE & j<=40){
        R <- R + floor(R/5)
      }
      
      if(cont==TRUE & j>=41){
        R <- 10000
      }
      
      if(cont==FALSE){
        mu.new.se <- sqrt(diag2vech(solve(FIM[1:2,1:2], tol=.Machine$double.eps*1e-50)))
        EM.out[j,6:7] <- c(mu.new.se)
        WACov <- WA_Cov(Data = Data, sigma = sigma.new)
        mu.new.WA <- WACov$Muhat.W
        mu.new.se.WA <- as.numeric(sqrt(diag2vech(WACov$Cov_Muhat.W)))
      }
      
      #print(list(EM.out[((j-2):j),], "R"=R, "j"=j-1))
    }
    
    Converged <- j<maxit
    Out$Converged <- Converged
    Out$R <- R
    Out$EM <- EM.out[1:j,]
    Out$Estimates <- c("mu1"=EM.out[j,1], "mu2"=EM.out[j,2], "sig11"=EM.out[j,3], "sig12"=EM.out[j,4], "sig22"=EM.out[j,5]) #EM.out[j,1:5]
    Out$SeSp <- ilogit(EM.out[j,1:2])
    Out$FIM <- FIM
    Out$StdErrors <- EM.out[j,6:7]
    Out$CIs <- ilogit(c("Se.lb"=EM.out[j,1]-1.96*EM.out[j,6], "Se.ub"=EM.out[j,1]+1.96*EM.out[j,6],
                        "Sp.lb"=EM.out[j,2]-1.96*EM.out[j,7], "Sp.ub"=EM.out[j,2]+1.96*EM.out[j,7]))
    Out$WACoV <- WACov$Cov_Muhat.W
    Out$mu.new.WA <- mu.new.WA
    Out$mu.new.se.WA <- mu.new.se.WA
    Out$logLik <- EM.out[,13]
    Out$Stats <- c("logLik" = Out$logLik[j],
                   "AIC" = -2*(Out$logLik[j]) + 2*5,
                   "BIC" = -2*(Out$logLik[j]) + log(2*k)*5)
    
    save(Out, file = paste("Scen", Scen, "SimResNov1.RData", sep = "_"))
    
    return(Out)
    
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  
}



#---------------------------------------------------------------------------------------------
# II. The "Create.Figure()" function is used to produce Figures presented in
# Section 4 of the paper.
#---------------------------------------------------------------------------------------------

Create.Figure <- function(mat.scen, Outl, Sigsize, fit.bnn, fit.bbn, fit.bnlm){
  # 'mat.scen' is the 16 by 14 matrix representing the matrix of scenario 
  # 'Outl' a character representing the parameter in which outlying studies are introduced (i.e., "Se and Sp", "Se", "Sp", "None") 
  # 'Sigsize' which is either 'large' or 'small' denotes the type of Sigma used to simulate outliers 
  # 'fit.bnn' is a list containing the BNN model's fit to the data
  # 'fit.bbn' is a list containing the BBN model's fit to the data
  # 'fit.bnlm' is a list containing the BNLM model's fit to the data
  
  
  #filenames <- c("OutlSeSp.pdf", "OutlSe.pdf", "OutlSp.pdf", "NoOutl.pdf")
  filenames <- c("OutlSeSp_large.tif", "OutlSe_large.tif", "OutlSp_large.tif", "OutlSeSp_small.tif", "OutlSe_small.tif", "OutlSp_small.tif", "NoOutl.tif")
  ylims <- list(c(0,1,0,1), c(0,1,0,1), c(0,1,0,1), c(0,1,0,1), c(0,1,0,1), c(0,1,0,1), c(0,1,0,1))
  ylabs <- list(c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"),
                c("Sensitivity", "Specificity"))
  
  legdims <- list(c(10,0.4,10,0.4),c(10,.4,10,.4),c(10,.4,10,.4),c(10,.4,10,.4),c(10,.4,10,.4),c(10,.4,10,.4),c(10,.4,10,.4))
  
  if(Outl=="Se and Sp" & Sigsize=="large"){
    i=1
  }else if(Outl=="Se and Sp" & Sigsize=="small"){
    i=13
  }else if(Outl=="Se" & Sigsize=="large"){
    i=5
  }else if(Outl=="Se" & Sigsize=="small"){
    i=17
  }else if(Outl=="Sp" & Sigsize=="large"){
    i=9
  }else if(Outl=="Sp" & Sigsize=="small"){
    i=21
  }else if(Outl=="None" & Sigsize=="None"){
    i=25
  }
  
  windowsFonts(Font=windowsFont(family = "Times New Roman"))
  
  # Plot results
  tiff(filename = filenames[ceiling(i/4)], width = 12, height = 6, compression = "lzw", res=300, units="in")
  par(mfrow=c(1,2))
  
  # For Se
  plot(mat.scen[(i:(i+3)),13], c(fit.bnlm[[i]]$SeSp[1],fit.bnlm[[(i+1)]]$SeSp[1],fit.bnlm[[(i+2)]]$SeSp[1],fit.bnlm[[(i+3)]]$SeSp[1]),
       type="b", pch=0, lty=1, ylim = ylims[[ceiling(i/4)]][1:2], ylab = "", xlab = "", xaxt="n")
  title(ylab = ylabs[[ceiling(i/4)]][1], xlab = "Number of studies", mgp=c(2,2,0), cex.lab=1.2, main = "Estimates of Sensitivity", mgp=c(2,2,0), cex.lab=1.2, family="Font")
  axis(1, at=mat.scen[(i:(i+3)),13])
  
  lines(mat.scen[(i:(i+3)),13], c(fit.bnn[[i]]$SeSp[1],fit.bnn[[(i+1)]]$SeSp[1],fit.bnn[[(i+2)]]$SeSp[1],fit.bnn[[(i+3)]]$SeSp[1]), type = "b", pch=1, lty=1)
  lines(mat.scen[(i:(i+3)),13], c(fit.bbn[[i]]$SeSp[1],fit.bbn[[(i+1)]]$SeSp[1],fit.bbn[[(i+2)]]$SeSp[1],fit.bbn[[(i+3)]]$SeSp[1]), type = "b", pch=2, lty=1)
  
  legend(x=legdims[[ceiling(i/4)]][1], y=legdims[[ceiling(i/4)]][2], legend = c("BBN", "BNN","BNL"), pch = c(2:0), bty="n")
  
  
  
  # For Sp
  plot(mat.scen[(i:(i+3)),13], c(fit.bnlm[[i]]$SeSp[2],fit.bnlm[[(i+1)]]$SeSp[2],fit.bnlm[[(i+2)]]$SeSp[2],fit.bnlm[[(i+3)]]$SeSp[2]),
       type = "b", pch=0, lty=1, ylim = ylims[[ceiling(i/4)]][3:4], ylab = "", xlab = "", xaxt="n")
  title(ylab = ylabs[[ceiling(i/4)]][2], xlab = "Number of studies", mgp=c(2,2,0), cex.lab=1.2, main = "Estimates of Specificity", family="Font")
  axis(1, at=mat.scen[(i:(i+3)),13])
  
  lines(mat.scen[(i:(i+3)),13], c(fit.bnn[[i]]$SeSp[2],fit.bnn[[(i+1)]]$SeSp[2],fit.bnn[[(i+2)]]$SeSp[2],fit.bnn[[(i+3)]]$SeSp[2]), type = "b", pch=1, lty=1)
  lines(mat.scen[(i:(i+3)),13], c(fit.bbn[[i]]$SeSp[2],fit.bbn[[(i+1)]]$SeSp[2],fit.bbn[[(i+2)]]$SeSp[2],fit.bbn[[(i+3)]]$SeSp[2]), type = "b", pch=2, lty=1)
  
  legend(x=legdims[[ceiling(i/4)]][3],y=legdims[[ceiling(i/4)]][4], legend = c("BBN", "BNN","BNL"), pch = c(2:0), bty="n")
  
  dev.off()
  
}


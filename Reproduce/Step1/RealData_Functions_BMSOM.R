#----------------------------------------------------------------------------------------------
# This R script contains code to demonstrate how to identify outlying studies in meta-analysis
# of diagnostic test accuracy studies using a new bivariate mean shift outlier model (BMSOM). 
#----------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------
# The "comp.y()" function is used to compute the observed logit sensitivities and logit specificities
#------------------------------------------------------------------------------------------------------
comp.y <- function(Data){
  # Data is the kby4 DTA data provided as a data frame or matrix specifying the column names "TP", "FP", "TN" and "FN" clearly.
  
  require(mada); require(faraway)
  k <- nrow(Data)
  # A continuity correction of 0.5 will be added to a cell which appears to contain 0, since logit(p) won't be defined otherwise
  
  y <- matrix(c(logit(madad(Data, correction.control = "single")$sens$sens), 
                logit(madad(Data, correction.control = "single")$spec$spec)), nrow=k, ncol=2)
  return(y)
}


#------------------------------------------------------------------------------------------------
# The "comp.Psi()" function is used to compute the within-study variances of logit Se and logit Sp
# using the theory of delta-method.
#------------------------------------------------------------------------------------------------
comp.Psi <- function(Data){
  # "Data" is the kby4 data frame containing the number of "TP", "TN", "FP" and "FN".
  library(mada); library(faraway)
  k <- nrow(Data); p <- 2
  correction <- 0.5
  Y <- comp.y(Data=Data)
  
  Psi <- array(0, dim = c(p,p,k))
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + correction
      Data$FP[i] = Data$FP[i] + correction
      Data$FN[i] = Data$FN[i] + correction
      Data$TN[i] = Data$TN[i] + correction
    }
  }
  
  n1 <- Data$TP+Data$FN; n2 <- Data$FP+Data$TN
  Psi[1,1,] <- 1/(n1*(ilogit(Y[,1])*(1-ilogit(Y[,1]))))
  Psi[2,2,] <- 1/(n2*(ilogit(Y[,2])*(1-ilogit(Y[,2]))))
  
  return(Psi)
}


#--------------------------------------------------------------------------------------------
# I. The "fit_bnn()" function is used to fit the BNN model of Reitsma et al. (2005)
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
  
  # Obtain backtransformed proportions
  Se <- ilogit(Y[,1]); Sp <- ilogit(Y[,2])
  
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
  
  # Get the Jacobian of the logit transformation that needs to be added to
  # the likelihood of the BNN-logit model to compare the three models using
  # information criteria (AIC/BIC).
  Jacobian.logit <- matrix(0, ncol = 2, nrow = k)
  Jacobian.logit <- cbind((1/(Se*(1-Se))), (1/(Sp*(1-Sp))))
  
  # Finally, fit the BNN LMM using the logit transformation
  fit.logit <- list()
  mvmeta.logit <- mvmeta(cbind(trse, trsp)~1, S=Data.mvmeta.logit[, 4:6], method = method, data = Data.mvmeta.logit, control=list(maxiter=1000))
  
  # Backtransformed results
  require(faraway)
  summary.logit <- summary(mvmeta.logit)$coefficients
  SeSp <- c(ilogit(summary.logit[1,1]), ilogit(summary.logit[2,1]))
  SeCI <- as.numeric(ilogit(summary.logit[1,5:6]))
  SpCI <- as.numeric(ilogit(summary.logit[2,5:6]))
  
  loglik.bnn <- mvmeta.logit$logLik + sum(log(Jacobian.logit[,1])+log(Jacobian.logit[,2]))
  
  fit.logit$mvmeta.logit <- summary(mvmeta.logit)
  fit.logit$SeSp <- SeSp
  fit.logit$SeSp.logit <- logit(SeSp)
  fit.logit$Psi <- mvmeta.logit$Psi
  fit.logit$vcov <- mvmeta.logit$vcov
  fit.logit$CIs <- c("Se.lb"=SeCI[1],"Se.ub"=SeCI[2], "Sp.lb"=SpCI[1],"Sp.ub"=SpCI[2])
  fit.logit$Stats <- c(logLik=loglik.bnn, AIC=-2*loglik.bnn + 2*5)
  
  return(fit.logit)
  
}


#---------------------------------------------------------------------------------------
# The ml_logLik() function computes log-likelihood function for the proposed BMSOM
#---------------------------------------------------------------------------------------
ml_logLik <- function(pars, Data, j){
  # 'pars' is a vector of unknown parameters (mu1,mu2,eta1,eta2,sigma11,sigma12,sigma22)
  # of the proposed BCSOM.
  # 'Data' is the kby4 data frame containing the columns TP,FP,FN, and TN
  # 'j' is an index that represents for which study to shift the mean vector  
  
  if(!length(pars)==7) stop("The length of the starting value parameters should be seven")
  
  k <- nrow(Data)
  mu <- pars[1:2]
  eta <- pars[3:4]
  Sigma <- matrix(c(pars[5:6],pars[6:7]),2,2)
  
  
  require(matrixcalc)
  if(is.positive.definite(Sigma)==FALSE){# To ensure that the estimated covariace matrix is positive definite
    return(-Inf)
  }else{
    
    Y <- comp.y(Data)
    Phi <- comp.Psi(Data)
    etai <- matrix(0, nrow = k, ncol = 2)
    ll <- numeric(0)
    
    library(mvtnorm)
    for (i in 1:k) {
      if(i==j){
        etai[i,] <- eta
      }
      
      Sigmai <- Sigma + Phi[,,i]
      
      ll[i] <- dmvnorm(x=Y[i,], mean = mu+etai[i,], sigma = Sigmai, log = TRUE)
    }
    
    Out <- sum(ll)
    Out
    
  }
  
}


#-------------------------------------------------------------------------------------------
# The Sim_DTAData() function is used to generate DTA datasets assuming the bivariate normal 
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
  
  require(mvtnorm)
  require(faraway)
  
  transSeSp = btransSeSp <- matrix(0, nrow = k, ncol = 2)
  TP = FP = TN = FN <- numeric(0)
  Data <- matrix(0, nrow = k, ncol = 5)
  
  set.seed(seed)
  transSeSp <- rmvnorm(n=k, mean = mu, sigma = Sigma)
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
# The fit_bmsom() function is used to fit the proposed BMSOM using the ML method
#----------------------------------------------------------------------------------
fit_bmsom <- function(Data){
  # Data is the kbyn data frame where the columns TP, FN, FP and TN are specificied
  require(faraway); require(numDeriv)
  
  k <- nrow(Data)
  Y <- comp.y(Data = Data)
  Phi <- comp.Psi(Data = Data)
  
  # Obtain backtransformed proportions
  Se <- ilogit(Y[,1]); Sp <- ilogit(Y[,2])
  
  # Get the Jacobian of the logit transformation that needs to be added to
  # the likelihood of the BMSOM model to compare models using
  # information criteria (AIC/BIC).
  Jacobian.logit <- matrix(0, ncol = 2, nrow = k)
  Jacobian.logit <- cbind((1/(Se*(1-Se))), (1/(Sp*(1-Sp))))
  
  
  fit.bnn <- fit_bnn(Data=Data, method = "ml")
  muhat <- fit.bnn$SeSp.logit
  
  require(matrixcalc)
  if(is.positive.definite(fit.bnn$Psi)==TRUE){
    Sigmahat <- fit.bnn$Psi
  } else{
    Sigmahat <- diag(2)
  }
  
  pars <- c(muhat[1:2], 0, 0, Sigmahat[1,1],Sigmahat[1,2],Sigmahat[2,2])
  
  fit.bmsom <- list(); muhat.bmsom=etai <- matrix(0, nrow = k, ncol = 2) 
  Sigmahat.bmsom <- matrix(0, nrow = k, ncol = 3) 
  logLik=LRT=converged=pval <- rep(0,k)
  out= matrix(0, nrow = k, ncol = 10)
  
  for(i in 1:k){
    
    fit.bmsom[[i]] <- try(optim(par = pars, fn = ml_logLik, Data=Data, j=i,
                                method = "BFGS", control = list(fnscale=-1, maxit=1000)), silent = TRUE)
    
    if(class(fit.bmsom[[i]])=="try-error"){
      fit.bmsom[[i]] <- optim(par = pars, fn = ml_logLik, Data=Data, j=i,
                              method = "Nelder-Mead", control = list(fnscale=-1, maxit=10000))
    }
    
    converged[i] <- fit.bmsom[[i]]$convergence
    muhat.bmsom[i,] <- fit.bmsom[[i]]$par[1:2]
    etai[i,] <- fit.bmsom[[i]]$par[3:4]
    Sigmahat.bmsom[i,] <- c(fit.bmsom[[i]]$par[5:7])
    logLik[i] <- fit.bmsom[[i]]$value + sum(log(Jacobian.logit[,1])+log(Jacobian.logit[,2]))
    LRT[i] <- -2*(fit.bnn$Stats[1] - logLik[i])
    pval[i] <- pchisq(q=LRT[i], df=2, lower.tail = FALSE)
  }
  
  out <- round(cbind(muhat.bmsom, etai, Sigmahat.bmsom, logLik, LRT, converged, pval), 4)
  
  colnames(out) <- c("muhat1", "muhat2", "etahat1", "etahat2", "sigmahat11", "sigmahat12", "sigmahat22", "logLik", "LRT", "convergence", "p_value")
  return(list(fit.bmsom, out))
  
}



#-----------------------------------------------------------------------------------------------
# The fit_bmsom_boot() function is used to fit the BMSOM to B bootstrap samples, and it returns 
# the empirical distribution of the LRT statistic, LRT_j, j=1,2,...,k.
#-----------------------------------------------------------------------------------------------
fit_bmsom_boot <- function(Data, B){
  # Data is a kbyn dataframe containing the columns TP, FP, FN and TN
  # B stands for the bootstrap size
  
  k <- nrow(Data)
  p.vals <- rep(0,k)
  n1 <- Data$TP+Data$FN; n2 <- Data$TN+Data$FP
  LRT.boot=LRT.boot.sort <- matrix(0, nrow = B, ncol = k)
  out <- list()
  
  fit.bnn <- fit_bnn(Data = Data, method = "ml")
  fit.bmsom <- fit_bmsom(Data = Data)
  
  for(b in 1:B){
    DTA.boot <- Sim_DTAData(mu=fit.bnn$SeSp.logit, Sigma = fit.bnn$Psi, n1=n1, n2=n2, k=k, seed=b)
    fit.bnn.boot <- fit_bnn(Data = DTA.boot, method = "ml")
    fit.bmsom.boot <- fit_bmsom(Data = DTA.boot)
    for(i in 1:k){
      LRT.boot[b,i] <- -2*(fit.bnn.boot$Stats[1] - fit.bmsom.boot[[2]][i,8])
    }
    LRT.boot.sort[b,] <- sort(LRT.boot[b,], decreasing = TRUE)
  }
  
  out$LRT.obs <- fit.bmsom[[2]][,9]
  out$LRT.boot <- LRT.boot
  out$LRT.boot.sort <- LRT.boot.sort
  out
}


#-----------------------------------------------------------------------------------
# The plot_bmsom() function is used to produce different plots used to assess the 
# presence or absence of outlying studies in meta-analysis of diagnostic test 
# accuracy studies using the proposed BMSOM method.
#-----------------------------------------------------------------------------------
plot_bmsom <- function(model.boot, model.bmsom, prob=0.95, qs = 3){
  # 'model.boot' is an object returned by the fit_bmsom.boot() that comprises of 
  # LRT statistics obtained after fitting the BMSOM to the bootstrap samples
  # 'model.bmsom' is an object returned by the fit_bmsom() function.
  # 'prob' is the probability for which the to be computed quantile corresponds to 
  
  k <- nrow(model.bmsom[[2]])
  
  # Plot the BMSOM model parameters
  q <- numeric()
  for(i in 1:qs){
    q[i] <- quantile(x=model.boot$LRT.boot.sort[,i], probs = prob, type = 2)
  }
  
  par(mfrow=c(1,2))
  
  # 1. Plot the estimated shifted mean for logit(Se)
  plot(0, type="n", xlim = c(1,k), ylim = c(min(model.bmsom[[2]][,3])-1,max(model.bmsom[[2]][,3])+1), ylab = "", xlab = "", main = "The estimated mean shift for logit(Se)", xaxt="n")
  title(ylab = expression(paste(hat(eta)[1])), xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.bmsom[[2]][,3], type = "p")
  abline(h=0, lty = 3)
  axis(side = 1, at = 1:k, las=2)
  
  # 2. Plot the estimated shifted mean for logit(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(min(model.bmsom[[2]][,4])-1,max(model.bmsom[[2]][,4])+1), ylab = "", xlab = "", main = "The estimated mean shift for logit(Sp)", xaxt="n")
  title(ylab = expression(paste(hat(eta)[2])), xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.bmsom[[2]][,4], type = "p")
  abline(h=0, lty = 3)
  axis(side = 1, at = 1:k, las=2)
  
  par(mfrow=c(1,1))
  
  # 3. Plot the observed LRT statistics and the corresponding thresholds obtained from bootstrap
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(model.boot$LRT.boot.sort[,1])), ylab = "", xlab = "", main = "Distribution of the likelihood ratio test statistic", xaxt="n")
  title(ylab = expression(paste("LRT"[i])), xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.bmsom[[2]][,9], type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q[i], q[i]), lty=qs+1-i)
  }
  
  
  tiff(filename = "LRT_Plot.tif", width=14, height=6, compression = "lzw", res=300, units="in")
  par(mfrow=c(1,1))
  # 3. Plot the observed LRT statistics and the corresponding thresholds obtained from bootstrap
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(model.boot$LRT.boot.sort[,1])), ylab = "", xlab = "", main = "Distribution of the likelihood ratio test statistic", xaxt="n")
  title(ylab = expression(paste("LRT"[i])), xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.bmsom[[2]][,9], type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q[i], q[i]), lty=qs+1-i)
  }
  dev.off()
  
  
}

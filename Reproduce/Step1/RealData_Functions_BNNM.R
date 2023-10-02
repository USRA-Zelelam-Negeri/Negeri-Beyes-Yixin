#------------------------------------------------------------------------------------------
# This R script contains code that demonstrates how to identify outlying and influential
# studies in meta-analysis of diagnostic test accuracy studies using new proposed methods
# that are based on the BNN model of Reitsma et al. (2005).
#------------------------------------------------------------------------------------------


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


#--------------------------------------------------------------------------------------------------
# The "comp.Psi()" function is used to compute the within-study variances of logit Se and logit Sp
# which is based on the theory of the delta-method
#--------------------------------------------------------------------------------------------------
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
# The "fit_bnn()" function is used to fit the BNN model of Reitsma et al. (2005)
#--------------------------------------------------------------------------------------------
fit_bnn <- function(Data, method=c("reml", "ml")){
  # "Data" is a list of data frames consisting of the 2X2 DTA table from several studies.
  # The column names of "Data" must include the frequencies "TP", "FN", "FP" and "TN".
  
  Data1 <- Data
  
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
  fit.logit$Data <- Data1
  fit.logit$Y <- Y
  fit.logit$n1 <- n1
  fit.logit$n2 <- n2
  fit.logit$SeSp <- SeSp
  fit.logit$SeSp.logit <- logit(SeSp)
  fit.logit$Psi <- mvmeta.logit$Psi
  fit.logit$vcov <- mvmeta.logit$vcov
  fit.logit$CIs <- c("Se.lb"=SeCI[1],"Se.ub"=SeCI[2], "Sp.lb"=SpCI[1],"Sp.ub"=SpCI[2])
  fit.logit$Stats <- c(logLik=loglik.bnn, AIC=-2*loglik.bnn + 2*5)
  
  return(fit.logit)
  
}


#-----------------------------------------------------------------------------------------
# The function diagnostics.dta() computes the standardized deleted residuals the proposed
# influence metrics meta-analysis of DTA studies using the BNN model.
#---------------------------------------------------------------------------------------
diagnostics.dta <- function(model, method=c("reml", "ml")){
  # "model" is an object returned by the "fit.bnn()" function.
  # "method" specifies the method of parameter estimation to be used
  
  require(faraway)
  
  Data <- model$Data
  k <- nrow(Data)
  Y <- model$Y
  n1 <- model$n1; n2 <- model$n2
  SeSp.logit <- model$SeSp.logit
  
  # Compute the within-study covariance matrices
  Phi <- array(data = 0, dim = c(2,2,k))
  Phi[1,1,] <- 1/(n1*(ilogit(Y[,1])*(1-ilogit(Y[,1]))))
  Phi[2,2,] <- 1/(n2*(ilogit(Y[,2])*(1-ilogit(Y[,2]))))
  
  # "mat.sqrt()" computes the square root of a positive-definite matrix
  mat.sqrt <- function(M){
    eig <- eigen(M)
    if(all(eig$values>c(0,0)))
      M.sqrt <- eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
    else
      stop("The matrix is not positive-definite")
    return(M.sqrt)
  }
  
  #----------------------------
  # I. Compute residuals
  #----------------------------
  
  # 1. Compute the raw residuals
  rawres <- sweep(Y, 2, model$SeSp.logit, "-")
  colnames(rawres) <- c("rawres_trSe", "rawres_trSp")
  
  # 2. Compute the internally studentized (or standardized) residuals
  # Compute covariance of the raw residuals
  covrawres <- list(); studres <- matrix(0, nrow=k, ncol=2)
  for(i in 1:k){
    covrawres[[i]] <- model$Psi + Phi[,,i] - model$vcov
    studres[i,] <- c(solve(mat.sqrt(covrawres[[i]]))%*%rawres[i,])
  }
  colnames(studres) <- c("Studres_trSe", "Studres_trSp")
  
  # 3. Compute the externally studentized (or deleted) residuals
  delres=del <- matrix(0, nrow = k, ncol = 2)
  fit=covdelres <- list()
  for(i in 1:k){
    fit[[i]] <- fit_bnn(Data = Data[-i,], method = method)
    del[i,] <- Y[i,] - fit[[i]]$SeSp.logit
    covdelres[[i]] <- fit[[i]]$Psi + Phi[,,i] + fit[[i]]$vcov
    delres[i,] <- c(solve(mat.sqrt(covdelres[[i]]))%*%del[i,])
  }
  colnames(delres) <- c("delres_trSe", "delres_trSp")
  pval.se <- 2*(1-pnorm(abs(delres[,1])))
  pval.sp <- 2*(1-pnorm(abs(delres[,2])))
  
  
  #----------------------------------
  # II. Compute influence measures
  #----------------------------------
  
  # 1. Compute the overall variance-covariance ratio (VARCOVRATIO)
  VARCOVRATIOi <- rep(0, k)
  for(i in 1:k){
    VARCOVRATIOi[i] <- det(fit[[i]]$vcov)/det(model$vcov)
  }
  names(VARCOVRATIOi) <- c("VARCOVRATIOi")
  
  # 2. Compute the between-study covariance ratio (SIGMARATIO)
  
  SIGMARATIOi <- rep(0, k)
  for(i in 1:k){
    SIGMARATIOi[i] <- det(fit[[i]]$Psi)/det(model$Psi)
  }
  names(SIGMARATIOi) <- c("SIGMARATIOi")
  
  return(list(Se=data.frame("Study"=1:k, "rawres"=rawres[,1], "studres"=studres[,1], "delres"=delres[,1], "pValSe" = pval.se,  "VARCOVRATIO"=VARCOVRATIOi, "SIGMARATIO"=SIGMARATIOi),
              Sp=data.frame("Study"=1:k, "rawres"=rawres[,2], "studres"=studres[,2], "delres"=delres[,2], "pValSp" = pval.sp,  "VARCOVRATIO"=VARCOVRATIOi, "SIGMARATIO"=SIGMARATIOi)))
}



#----------------------------------------------------------------------------------------------
# The "plotdiag.dta()" function plots the proposed residuals and influence measures 
# that are useful for identify potential outlying and influential studies of MA of DTA studies
#----------------------------------------------------------------------------------------------
plotdiag.dta <- function(Summary){
  # Summary is an object returned by the "diagnostics.dta()" function
  
  # We produce the following figures
  
  par(mfrow=c(1,2))
  
  # 1. Plot standardized deleted residuals for tr(Se)
  plot(Summary$Se$Study, Summary$Se$delres, ylim = c(min(Summary$Se$delres)-1,max(Summary$Se$delres)+1), type = 'b', col='black', pch=16, xlab = "", ylab = "", xaxt = "n")
  title(main = "Standardized residuals for logit(Se)", xlab = "Study number", ylab = "Studentized deleted residuals")
  axis(side = 1, at = 1:length(Summary$Se$delres), las=2)
  abline(h=c(-1.96, 1.96), lty=3)
  abline(h=0, lty="dashed")
  
  # 2. Plot standardized deleted residuals for tr(Sp)
  plot(Summary$Sp$Study, Summary$Sp$delres, ylim = c(min(Summary$Sp$delres)-1,max(Summary$Sp$delres)+1), type = 'b', col='black', pch=16, xlab = "", ylab = "", xaxt = "n")
  title(xlab = "Study number", ylab = "Studentized deleted residuals", main = "Standardized residuals for logit(Sp)")
  axis(side = 1, at = 1:length(Summary$Sp$delres), las=2)
  abline(h=c(-1.96, 1.96), lty=3)
  abline(h=0, lty="dashed")
  
  # 3. Plot VARCOVRATIO
  plot(Summary$Se$Study, Summary$Se$VARCOVRATIO, ylim = c(min(Summary$Se$VARCOVRATIO)-1,max(Summary$Se$VARCOVRATIO)+1), type = 'b', col='black', pch=16, xlab = "", ylab = "", xaxt = "n")
  title(xlab = "Study number", ylab = "VARCOVRATIO", main = "VARCOVRATIO for logit(Se) and logit(Sp)")
  axis(side = 1, at = 1:length(Summary$Se$VARCOVRATIO), las=2)
  abline(h=1, lty="dashed")
  
  # 4. Plot SIGMARATIO 
  plot(Summary$Se$Study, Summary$Se$SIGMARATIO, ylim = c(min(Summary$Se$SIGMARATIO)-1,max(Summary$Se$SIGMARATIO)+1), type = 'b', col='black', pch=16, xlab="", ylab="", xaxt="n")
  title(xlab = "Study number", ylab = "SIGMARATIO", main = "SIGMARATIO for logit(Se) and logit(Sp)")
  axis(side = 1, at = 1:length(Summary$Se$SIGMARATIO), las=2)
  abline(h=1, lty="dashed")
  
}


#------------------------------------------------------------------------------------------
# The "Sim_DTAData()" function is used to generate DTA datasets assuming bivariate normal
# distribution for the study-specific observed pair of logit(Se) and logit(Sp). 
#------------------------------------------------------------------------------------------
Sim_DTAData <- function(mu, Sigma, n1, n2, k, seed){
  # mu is the vector of true means
  # Sigma is the true between-study covariance matrix
  # n1 is the study-specific number of diseased subjects
  # n2 is the study-specific number of non-diseased subjects
  # k is the number of studies in the MA
  # seed is a random seed to be used when generating 
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


#--------------------------------------------------------------------------------------
# The "influence.bootstrp()" function implements a parametric bootstrap approach 
# to testing whether study i, i=1,2,...,k, is potentially outlying and/or influential.
#--------------------------------------------------------------------------------------
influence.bootstrp <- function(stats, Data, method, B){
  # "stats" is an R object returned by the "diagnostics.dta()" function
  # "Data" is a kbyp dataframe comprising of k 2by2 DTA data
  # "method" stands for the method of estimation (ML or REML) to be used
  # "B" is the bootstrap size
  
  k <- nrow(Data)
  
  # Extract the observed "delres", "VARCOVRATIO" and "SIGMARATIO" for each study
  delres.se.obs <- stats$Se$delres; delres.sp.obs <- stats$Sp$delres; 
  VARCOVRATIOi.obs <- stats$Se$VARCOVRATIO; SIGMARATIOi.obs <- stats$Se$SIGMARATIO
  
  # Extract the estimated mean vector and between-study covariance matrix
  fit.bnn.obs <- fit_bnn(Data = Data, method = method)
  muhat.obs <- fit.bnn.obs$SeSp.logit
  Sigmahat.obs <- fit.bnn.obs$Psi
  n1 <- Data$TP + Data$FN; n2 <- Data$FP + Data$TN
  
  # Next generate the bootsrap samples, fit the BNN model to the generated data
  # and compute the bootstrap "VARCOVRATIO" and "SIGMARATIO" test statistics to 
  # obtain their respective sampling distributions and p-values.
  
  VARCOVRATIOi.boot = SIGMARATIOi.boot = delres.se.boot = delres.sp.boot <- matrix(0, nrow = B, ncol = k)
  for (b in 1:B){
    set.seed(b^2)
    DTA.boot <- Sim_DTAData(mu=muhat.obs, Sigma = Sigmahat.obs, n1=n1, n2=n2, k=k, seed=b)
    fit.bnn.boot <- fit_bnn(Data = DTA.boot, method = method)
    stats.boot <- diagnostics.dta(model = fit.bnn.boot, method = method)
    VARCOVRATIOi.boot[b,] <- stats.boot$Se$VARCOVRATIO
    SIGMARATIOi.boot[b,] <- stats.boot$Se$SIGMARATIO
    delres.se.boot[b,] <- stats.boot$Se$delres
    delres.sp.boot[b,] <- stats.boot$Sp$delres
  }
  
  # When the estimated between-study covariance matrix is not positive-definite, SIGMARATIO
  # becomes undefined since det(\Sigma) equals 0 or a negative number. For these scenario, we
  # remove the corresponding SIGMARATIO.
  if(!identical(which(SIGMARATIOi.boot[,1]==Inf | SIGMARATIOi.boot[,1]==-Inf | SIGMARATIOi.boot[,1]=="NaN"), integer(0))){
    SIGMARATIOi.boot <- SIGMARATIOi.boot[-c(which(SIGMARATIOi.boot[,1]==Inf  | SIGMARATIOi.boot[,1]==-Inf | SIGMARATIOi.boot[,1]=="NaN")),]
  }
  
  
  Out <- list()
  Out$VARCOVRATIOi.obs <- VARCOVRATIOi.obs
  Out$SIGMARATIOi.obs <- SIGMARATIOi.obs
  Out$delres.se.obs <- delres.se.obs
  Out$delres.sp.obs <- delres.sp.obs
  Out$VARCOVRATIOi.boot <- VARCOVRATIOi.boot
  Out$SIGMARATIOi.boot <- SIGMARATIOi.boot
  Out$delres.se.boot <- delres.se.boot
  Out$delres.sp.boot <- delres.sp.boot
  Out
}


#-----------------------------------------------------------------------------------
# The plot_influence() function produces different plots used to assess the presence
# or absence of outlying or influential studies in meta-analysis of diagnostic test 
# accuracy studies.
#-----------------------------------------------------------------------------------
plot_influence <- function(model.boot, prob, qs = 3){
  # 'model.boot' is an object returned by the influence.bootstp() that comprises of LRT statistics obtained
  # after fitting the BNN to the bootstrap samples
  # 'prob' is the probability for which the to be computed quantile corresponds to 
  # 'qs' denotes the number of thresholds to be plotted to identify up to 'qs' number of outliers
  
  k <- length(model.boot$VARCOVRATIOi.obs)
  
  delres.se.boot.sort=delres.sp.boot.sort=varcovratio.boot.sort <- matrix(0, nrow = nrow(model.boot$VARCOVRATIOi.boot), ncol = k)
  sigmaratio.boot.sort <- matrix(0, nrow = nrow(model.boot$SIGMARATIOi.boot), ncol = k)
  for(b in 1:nrow(model.boot$VARCOVRATIOi.boot)){
    delres.se.boot.sort[b,] <- sort(model.boot$delres.se.boot[b,], decreasing = TRUE)
    delres.sp.boot.sort[b,] <- sort(model.boot$delres.sp.boot[b,], decreasing = TRUE)
    varcovratio.boot.sort[b,] <- sort(model.boot$VARCOVRATIOi.boot[b,], decreasing = FALSE)
  }
  
  for (b in 1:nrow(model.boot$SIGMARATIOi.boot)) {
    sigmaratio.boot.sort[b,] <- sort(model.boot$SIGMARATIOi.boot[b,], decreasing = FALSE)
  }
  
  q.sedelres=q.spdelres=q.varcov=q.sigma <- numeric()
  for(i in 1:qs){
    q.sedelres[i] <- quantile(x=abs(delres.se.boot.sort[,i]), probs = prob, type = 2)
    q.spdelres[i] <- quantile(x=abs(delres.sp.boot.sort[,i]), probs = prob, type = 2)
    q.varcov[i] <- quantile(x=varcovratio.boot.sort[,i], probs = prob, type = 2)
    q.sigma[i] <- quantile(x=sigmaratio.boot.sort[,i], probs = prob, type = 2)
  }
  
  
  tiff(filename = "Std_dlt_residuals.tif", width=14, height=6, compression = "lzw", res=300, units="in")
  par(mfrow=c(1,2))
  # 1. Plot standardized deleted residuals for tr(Se)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(abs(delres.se.boot.sort[,1]))), ylab = "", xlab = "", main = "Standardized deleted residuals for logit(Se)", xaxt="n")
  title(ylab = "Studentized deleted residuals", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, abs(model.boot$delres.se.obs), type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.sedelres[i], q.sedelres[i]), lty=qs+1-i)
  }
  
  # 2. Plot standardized deleted residuals for tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(abs(delres.sp.boot.sort[,1]))), ylab = "", xlab = "", main = "Standardized deleted residuals for logit(Sp)", xaxt="n")
  title(ylab = "Studentized deleted residuals", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, abs(model.boot$delres.sp.obs), type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.spdelres[i], q.spdelres[i]), lty=qs+1-i)
  }
  dev.off()
  
  
  tiff(filename = "VARCOV_SIGMA_RATIO.tif", width=14, height=6, compression = "lzw", res=300, units="in")
  par(mfrow=c(1,2))
  # 3. Plot VARCOVRATIOi for tr(Se) and tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,2), ylab = "", xlab = "", main = "VARCOVRATIO for logit(Se) and logit(Sp)", xaxt="n")
  title(ylab = "VARCOVRATIO", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.boot$VARCOVRATIOi.obs, type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.varcov[i], q.varcov[i]), lty=qs+1-i)
  }
  
  # 4. Plot SIGMARATIOi for tr(Se) and tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,2), ylab = "", xlab = "", main = "SIGMARATIO for logit(Se) and logit(Sp)", xaxt="n")
  title(ylab = "SIGMARATIO", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.boot$SIGMARATIOi.obs, type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.sigma[i], q.sigma[i]), lty=qs+1-i)
  }
  dev.off()
  
  
  par(mfrow=c(1,2))
  
  # 1. Plot standardized deleted residuals for tr(Se)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(abs(delres.se.boot.sort[,1]))), ylab = "", xlab = "", main = "Standardized deleted residuals for logit(Se)", xaxt="n")
  title(ylab = "Studentized deleted residuals", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, abs(model.boot$delres.se.obs), type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.sedelres[i], q.sedelres[i]), lty=qs+1-i)
  }
  
  # 2. Plot standardized deleted residuals for tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,max(abs(delres.sp.boot.sort[,1]))), ylab = "", xlab = "", main = "Standardized deleted residuals for logit(Sp)", xaxt="n")
  title(ylab = "Studentized deleted residuals", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, abs(model.boot$delres.sp.obs), type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.spdelres[i], q.spdelres[i]), lty=qs+1-i)
  }
  
  
  # 3. Plot VARCOVRATIOi for tr(Se) and tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,2), ylab = "", xlab = "", main = "VARCOVRATIO for logit(Se) and logit(Sp)", xaxt="n")
  title(ylab = "VARCOVRATIO", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.boot$VARCOVRATIOi.obs, type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.varcov[i], q.varcov[i]), lty=qs+1-i)
  }
  
  # 4. Plot SIGMARATIOi for tr(Se) and tr(Sp)
  plot(0, type="n", xlim = c(1,k), ylim = c(0,2), ylab = "", xlab = "", main = "SIGMARATIO for logit(Se) and logit(Sp)", xaxt="n")
  title(ylab = "SIGMARATIO", xlab = "Study number", mgp=c(2,2,0), cex.lab=1.2)
  lines(1:k, model.boot$SIGMARATIOi.obs, type = "p")
  axis(side = 1, at = 1:k, las=2)
  for (i in 1:qs){
    lines(c(1-1,k+1), c(q.sigma[i], q.sigma[i]), lty=qs+1-i)
  }
  
}



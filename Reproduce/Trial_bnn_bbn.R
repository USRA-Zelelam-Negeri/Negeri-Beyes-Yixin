NEE <-read.csv(file = "/Users/a/Desktop/RA_Aggregate data meta analysis/SupplementaryMaterial/NEE.csv", header = TRUE)
NEE <- data.frame("Study"=paste("Study", 1:nrow(NEE), sep = " "), NEE)

MMSE <- read.table(file = "/Users/a/Desktop/RA_Aggregate data meta analysis/SupplementaryMaterial/MMSE.txt", header = TRUE)
MMSE <- data.frame("Study"=paste("Study", 1:nrow(MMSE), sep = " "), MMSE)

library("mada")

bnn <- function(Data, method=c("reml", "ml")){
  # Each row of data is a study containing TP,FN,FP and TN
  # detach("package:mada", unload=TRUE)
  library(metafor)
  require(mvmeta)
  
  k <- nrow(Data) # number of studies
  
  add_correct <- 0.1
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + add_correct
      Data$FP[i] = Data$FP[i] + add_correct
      Data$FN[i] = Data$FN[i] + add_correct
      Data$TN[i] = Data$TN[i] + add_correct
    }
  }
  
  TP <- Data$TP
  FP <- Data$FP
  FN <- Data$FN
  TN <- Data$TN
  n1 <- TP + FN
  n2 <- FP + TN
  npat <- n1 + n2
  
  # Store the within-study covariance vector (a zero vector) for each study
  within_var <- rep(0,k)

  # store the logit transformed Se/(1-Sp)
  logit_sesp <- matrix(0, nrow = k, ncol = 2)
  logit_sesp <- cbind(tsens = escalc(measure="PLO", xi=TP, ni=n1)[,1],
                        tfpr = escalc(measure="PLO", xi=TN, ni=n2)[,1])
  
  # store the within-study variances of logit(Se) and logit (1-Sp)
  var.se.logit = var.sp.logit <- rep(0, k)
  var.se.logit <- escalc(measure="PLO", xi=TP, ni=n1)[,2]
  var.sp.logit <- escalc(measure="PLO", xi=TN, ni=n2)[,2]
  
  # df for mvmeta()
  df_mvmeta <- matrix(0, nrow = k, ncol=6)
  df_mvmeta <- data.frame(pat.num=npat, trse=logit_sesp[,1], trsp=logit_sesp[,2],
                          varse=var.se.logit, covsefpr=within_var, varsp=var.sp.logit)
  
  # fit the BNN using the logit transformation
  fit.logit <- list()
  mvmeta_logit <- mvmeta(cbind(trse, trsp)~1, S=df_mvmeta[, 4:6], method = method, 
                         data = df_mvmeta, control=list(maxiter=1000))
  
  require(faraway)
  summary.logit <- summary(mvmeta_logit)$coefficients
  SeSp <- c(ilogit(summary.logit[1,1]), ilogit(summary.logit[2,1]))
  SeCI <- as.numeric(ilogit(summary.logit[1,5:6]))
  SpCI <- as.numeric(ilogit(summary.logit[2,5:6]))
  
  loglik.bnn <- mvmeta_logit$logLik
  
  fit.logit$mvmeta_logit <- summary(mvmeta_logit)
  fit.logit$SeSp <- SeSp
  fit.logit$SeSp.logit <- logit(SeSp)
  fit.logit$Psi <- mvmeta_logit$Psi
  fit.logit$vcov <- mvmeta_logit$vcov
  fit.logit$CIs <- c("Se.lb"=SeCI[1],"Se.ub"=SeCI[2], "Sp.lb"=SpCI[1],"Sp.ub"=SpCI[2])
  fit.logit$Stats <- c(logLik=loglik.bnn, AIC=-2*loglik.bnn + 2*5, BIC=-2*loglik.bnn + 5*log(k*2))
  
  return(fit.logit)
}


bnn.nee <- bnn(Data = NEE, method = "ml")
bnn.mmse <- bnn(Data = MMSE, method = "ml")

data_prepare <- function(Data, Study){
  n <- length(Study)
  study_num <- rep(1:n, each=2)
  
  TP <- Data$TP
  FP <- Data$FP
  FN <- Data$FN
  TN <- Data$TN
  
  well_classified = c()
  mis_classified = c()
  groups = c()
  
  for (i in 1:n){
    well_classified <- append(well_classified, TP[i])
    well_classified <- append(well_classified, TN[i])
    mis_classified <- append(mis_classified, FN[i])
    mis_classified <- append(mis_classified, FP[i])
    groups <- append(groups, 'diseased')
    groups <- append(groups, 'healthy')
  }
  
  result <- data.frame(study_num, well_classified, mis_classified, groups)
  return(result)
}

bbn <- function(Data, Study){
  # study: study number from data

  require(metafor)
  require(mvmeta)
  require(Metatron)
  require(faraway)
  require(lme4)
  
  k <- nrow(Data)
  
  TP <- Data$TP
  FP <- Data$FP
  FN <- Data$FN
  TN <- Data$TN
  n1 <- TP + FN
  n2 <- FP + TN
  npat <- n1 + n2

  # fit the BBN model
  metatron.bbn <- fit.bivar(TP=TP, FN=FN, TN=TN, FP=FP, study=Study, data=Data)
  
  # Extract the parameter estimates and their covariance matrix / CIs when possible
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


bbn_new <- function(Data, Study){
  # study: study number from data
  
  require(metafor)
  require(mvmeta)
  require(faraway)
  require(lme4)
  
  k <- nrow(Data)
  
  TP <- Data$TP
  FP <- Data$FP
  FN <- Data$FN
  TN <- Data$TN
  n1 <- TP + FN
  n2 <- FP + TN
  npat <- n1 + n2
  
  # fit the BBN model
  r1 <- data_prepare(Data, Study)
  r2 <- glmer(cbind(well_classified, mis_classified) ~ groups-1+(groups-1|study_num), r1, binomial)
  
  loglik.bbn <- logLik(r2)
  fixed_coeff <- fixef(r2)
  logit_sensitivity <- fixed_coeff['groupsdiseased']
  logit_speciticity <- fixed_coeff['groupshealthy']
  fixed <- data.frame(logit_sensitivity, logit_speciticity)
  
  random_coeff <- ranef(r2)$study_num
  random_mat <- cov(random_coeff)
  
  fit.bbn <- list()
  
  fit.bbn$model.bbn <- r2
  fit.bbn$SeSp <- ilogit(as.numeric(fixed))
  fit.bbn$SeSp.logit <- as.numeric(fixed)
  fit.bbn$Psi <- random_mat[1:2,1:2]
  fit.bbn$vcov <- vcov(r2)
  fit.bbn$CIs <- c("Se.lb"=ilogit(fit.bbn$SeSp.logit[1]-qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[1,1])), "Se.ub"=ilogit(fit.bbn$SeSp.logit[1]+qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[1,1])),
                   "Sp.lb"=ilogit(fit.bbn$SeSp.logit[2]-qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[2,2])), "Sp.ub"=ilogit(fit.bbn$SeSp.logit[2]+qt(p=0.975, df=k-2)*sqrt(fit.bbn$vcov[2,2])))
  fit.bbn$Stats <- c(logLik=loglik.bbn, AIC= -2*loglik.bbn + 2*5, BIC = -2*loglik.bbn + 5*log(k*2))
  
  return(fit.bbn)
}



bbbn.nee <- bbn(Data = NEE, Study = NEE$Study)
bbn.mmse <- bbn(Data = MMSE, Study = MMSE$Study)

# Results -- Table 1
bnn.nee$SeSp
bbn.nee$SeSp

bnn.nee$CIs
bbn.nee$CIs
 
bnn.nee$Stats
bbn.nee$Stats

# Results -- Table2
bnn.mmse$SeSp
bbn.mmse$SeSp

bnn.mmse$CIs
bbn.mmse$CIs

bnn.mmse$Stats
bbn.mmse$Stats


# Results -- From own implementation
bbn.nee_new <- bbn_new(Data = NEE, Study = NEE$Study)
bbn.mmse_new <- bbn_new(Data = MMSE, Study = MMSE$Study)

bbn.nee_new$SeSp
bbn.nee_new$CIs
bbn.nee_new$Stats

bbn.mmse_new$SeSp
bbn.mmse_new$CIs
bbn.mmse_new$Stats

 
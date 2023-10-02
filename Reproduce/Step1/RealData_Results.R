#-----------------------------------------------------------------------------#
# This R script illustrates the results in Section 2 and 5 of the manuscript  #
#-----------------------------------------------------------------------------#

##----------##
## Section 2##
##----------##

### 2.1 - NEE data
NEE <-read.csv(file = "/Users/a/Desktop/RA_Aggregate data meta analysis/SupplementaryMaterial/NEE.csv", header = TRUE)
colnames(NEE)[1] <- "Author"
NEE <- data.frame("Study"=paste("Study", 1:nrow(NEE), sep = " "), NEE)

### 2.1.1 Detecting outlying studies - residual-based approach
source("/Users/a/Desktop/RA_Aggregate data meta analysis/Reproduce/Step1/RealData_Functions_BNNM.R")

fit.nee <- fit_bnn(Data = NEE, method = "reml")
diag.nee <- diagnostics.dta(model = fit.nee, method = "reml")
system.time(infl.nee.boot <- influence.bootstrp(stats = diag.nee, Data = NEE, method = "reml", B=1000))
plot_influence(model.boot = infl.nee.boot, prob = 0.95, qs=3)
# Study 2 and 6 are identified to be outlying and influential in Sp; study
# 3 is identified as outlying and influential in Se.

### 2.1.2 Detecting outlying studies - LRT-based approach
source("/Users/a/Desktop/RA_Aggregate data meta analysis/Reproduce/Step1/RealData_Functions_BMSOM.R")

system.time(fit.bmsom.nee <- fit_bmsom(Data = NEE))
system.time(fit.bmsom.boot.nee <- fit_bmsom_boot(Data = NEE, B=100))
plot_bmsom(model.boot = fit.bmsom.boot.nee, model.bmsom = fit.bmsom.nee, prob = 0.95, qs=3)
# Study 2 and 6 are identified as potential outliers.
outl.nee <- c(2,6)

### 2.1.3 Forest Plots
library(mada)
tiff(filename = "NEE_Forestplot.tif", width=12, height=8, compression = "lzw", res=500, units="in")
#pdf("NEE_Forestplot.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
forest(madad(x=NEE, correction.control = "single", method = "wald"), type = "sens", main = "Forest plot of sensitivities", snames=NEE$Study)
text(-0.4, 10.15, "Study", pos = 2)
text(1.07, 10.15, "Sensitivity [95% CI]", pos=4)

forest(madad(x=NEE, correction.control = "single", method = "wald"), type = "spec", main = "Forest plot of specificities", snames=NEE$Study)
text(-0.8, 10.15, "Study", pos=2)
text(1.07, 10.15, "Specificity [95% CI]", pos=4)
dev.off()


### 2.2 Mini-mental state examination (MMSE) data
MMSE <- read.table(file = "MMSE.txt", header = TRUE)
MMSE <- data.frame("Study"=paste("Study", 1:nrow(MMSE), sep = " "), MMSE)

### 2.2.1 Detecting outlying studies - residual-based approach
source("RealData_Functions_BNNM.R")

fit.mmse <- fit_bnn(Data = MMSE, method = "reml")
diag.mmse <- diagnostics.dta(model = fit.mmse, method = "reml")
system.time(infl.mmse.boot <- influence.bootstrp(stats = diag.mmse, Data = MMSE, method = "reml", B=1000))
plot_influence(model.boot = infl.mmse.boot, prob = 0.95, qs=3)


### 2.2.2 Detecting outlying studies - LRT-based approach
source("RealData_Functions_BMSOM.R")

system.time(fit.bmsom.mmse <- fit_bmsom(Data = MMSE))
system.time(fit.bmsom.boot.mmse <- fit_bmsom_boot(Data = MMSE, B=100))
plot_bmsom(model.boot = fit.bmsom.boot.mmse, model.bmsom = fit.bmsom.mmse, prob = 0.95, qs=3)
outl.mmse <- c(1,2)

### 2.2.3 Forest Plots
tiff(filename = "MMSE_Forestplot.tif", width=12, height=8, compression = "lzw", res=500, units="in")
#pdf("MMSE_Forestplot.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
forest(madad(x=MMSE, correction.control = "single", method = "wald"), type = "sens", main = "Forest plot of sensitivities", snames=MMSE$Study)
text(-0.5, 9.15, "Study", pos = 2)
text(1.0, 9.15, "Sensitivity [95% CI]", pos=4)

forest(madad(x=MMSE, correction.control = "single", method = "wald"), type = "spec", main = "Forest plot of specificities", snames=MMSE$Study)
text(-0.3, 9.15, "Study", pos=2)
text(1.0, 9.15, "Specificity [95% CI]", pos=4)
dev.off()


##----------##
## Section 5##
##----------##

### 5.1 - NEE data
source("/Users/a/Desktop/RA_Aggregate data meta analysis/Reproduce/Step1/RealData_Functions.R")
fit.bnlm.nee <- mcemDTA(Data = NEE, outl = outl.nee, maxit = 100)
fit.bnn.nee <- fit_bnn(Data = NEE, method = "ml")
fit.bbn.nee <- fit_bbn(Data = NEE, Study = NEE$Study)

fit.bnlm.nee$Converged

fit.bnlm.nee$SeSp; fit.bnn.nee$SeSp; fit.bbn.nee$SeSp

(fit.bnlm.nee$SeSp-fit.bnn.nee$SeSp)*100; (fit.bnlm.nee$SeSp-fit.bbn.nee$SeSp)*100
(fit.bnlm.nee$SeSp[1]-fit.bnn.nee$SeSp[1])/fit.bnlm.nee$SeSp[1]*100; (fit.bnlm.nee$SeSp[2]-fit.bnn.nee$SeSp[2])/fit.bnlm.nee$SeSp[2]*100
(fit.bnlm.nee$SeSp[1]-fit.bbn.nee$SeSp[1])/fit.bnlm.nee$SeSp[1]*100; (fit.bnlm.nee$SeSp[2]-fit.bbn.nee$SeSp[2])/fit.bnlm.nee$SeSp[2]*100

fit.bnlm.nee$CIs; fit.bnn.nee$CIs; fit.bbn.nee$CIs
fit.bnlm.nee$CIs[2]-fit.bnlm.nee$CIs[1]; fit.bnn.nee$CIs[2]-fit.bnn.nee$CIs[1]; fit.bbn.nee$CIs[2]-fit.bbn.nee$CIs[1]
fit.bnlm.nee$CIs[4]-fit.bnlm.nee$CIs[3]; fit.bnn.nee$CIs[4]-fit.bnn.nee$CIs[3]; fit.bbn.nee$CIs[4]-fit.bbn.nee$CIs[3]

fit.bnlm.nee$Stats; fit.bnn.nee$Stats; fit.bbn.nee$Stats


### 5.2 - MMSE data
fit.bnlm.mmse <- mcemDTA(Data = MMSE, outl = outl.mmse, maxit = 100)
fit.bnn.mmse <- fit_bnn(Data = MMSE, method = "ml")
fit.bbn.mmse <- fit_bbn(Data = MMSE, Study = MMSE$Author)

fit.bnlm.mmse$Converged

fit.bnlm.mmse$SeSp; fit.bnn.mmse$SeSp; fit.bbn.mmse$SeSp

fit.bnlm.mmse$SeSp-fit.bnn.mmse$SeSp; fit.bnlm.mmse$SeSp-fit.bbn.mmse$SeSp
(fit.bnlm.mmse$SeSp[1]-fit.bnn.mmse$SeSp[1])/fit.bnlm.mmse$SeSp[1]*100; (fit.bnlm.mmse$SeSp[2]-fit.bnn.mmse$SeSp[2])/fit.bnlm.mmse$SeSp[2]*100
(fit.bnlm.mmse$SeSp[1]-fit.bbn.mmse$SeSp[1])/fit.bnlm.mmse$SeSp[1]*100; (fit.bnlm.mmse$SeSp[2]-fit.bbn.mmse$SeSp[2])/fit.bnlm.mmse$SeSp[2]*100

fit.bnlm.mmse$CIs; fit.bnn.mmse$CIs; fit.bbn.mmse$CIs
fit.bnlm.mmse$CIs[2]-fit.bnlm.mmse$CIs[1]; fit.bnn.mmse$CIs[2]-fit.bnn.mmse$CIs[1]; fit.bbn.mmse$CIs[2]-fit.bbn.mmse$CIs[1]
fit.bnlm.mmse$CIs[4]-fit.bnlm.mmse$CIs[3]; fit.bnn.mmse$CIs[4]-fit.bnn.mmse$CIs[3]; fit.bbn.mmse$CIs[4]-fit.bbn.mmse$CIs[3]

fit.bnlm.mmse$Stats; fit.bnn.mmse$Stats; fit.bbn.mmse$Stats

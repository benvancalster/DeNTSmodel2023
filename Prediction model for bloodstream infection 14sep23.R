##############################
##############################
## PREDICTION MODEL FOR NTS ##
##############################
##############################

# Ben Van Calster, KU Leuven
# June-September 2023

# I'm not the world's best coder, but maybe that makes it easier to follow ;-)


###############################
# LOAD PACKAGES AND FUNCTIONS #
###############################

library(DescTools)
library(gmodels)
library(plyr)
library(rms) # logistic regression, restricted cubic splines
library(logistf) # LR using Firth's correction
library(pmsampsize)
source("C:\\Ben\\Methodological research\\R2 from c for binary outcomes\\approximate_R2_function.R")
library(auRoc) # AUC/c statistic using 95% CI based on the logit transform
library(CalibrationCurves) # Calibration
library(boot) # bootstrapping, logit
library(locfit) # expit
library(ggplot2)
library(rmda) # Net Benefit and Decision Curve Analysis



###############
# DATA IMPORT #
###############

dents <- read.delim("C://Users//u0022832//OneDrive - KU Leuven//Datasets//Datasets Bieke Tack//230523_BT_HDT_dataset pred model.txt")



####################
# DATA PREPARATION #
####################


# Recode categorical variables to 0/1 or 0/1/2

dents$sex = 1*(dents$sex == "Feminin")
dents$season = 1*(dents$season == "wet")
dents$mh_diarrh = 1*(dents$mh_diarrh == "Oui")
dents$malnutrition = as.numeric(mapvalues(dents$malnutrition, c("No acute malnutrition","MAM","MAS"), c(0,1,2)))
dents$tachypnea = 1*(dents$tachypnea == "tachypnea")
dents$distress = 1*(dents$distress == "resp distress")
dents$hepaspl = 1*(dents$pe_splen=="Oui" | dents$pe_hepa =="Oui")
dents$malar = 0*(dents$negmicro=="Negative microscopy" & dents$allrecentpf=="No (very) recent Pf malaria") + 
              1*(dents$negmicro=="Negative microscopy" & dents$allrecentpf=="(Very) recent Pf malaria") +
              2*(dents$negmicro=="Positive microscopy")
dents$malar[dents$negmicro=="" | dents$allrecentpf==""] = NA
dents$anemia_corr = as.numeric(mapvalues(dents$anemia_corr, c("No anemia","Moderate anemia","Severe anemia"), c(0,1,2)))
dents$anemia_corr[dents$anemia_corr == ""] = NA
dents$n_antibiotics = 1*(dents$n_antibiotics == "At least one antibiotic")
dents$nts = 1*(dents$nts=="NTS")


# Make categorical versions of age and fever duration

dents$agemcat = 0 + 1*(dents$age_m_calc>=6) + 1*(dents$age_m_calc>=24)
dents$durcat = 0 + 1*(dents$dur_fev>3) 


# Log transform fever duration, positively skewed (experience tells me it works better like that for spline modeling)

dents$ldurfev = log(dents$dur_fev+1)
dents$l2durfev = log2(dents$dur_fev+1)


# Make dummy variables from categorical variables

dents$agef = as.factor(dents$agemcat)
dents$agef1 = 1*(dents$agemcat==1)
dents$agef2 = 1*(dents$agemcat==2)
dents$malnutf = as.factor(dents$malnutrition)
dents$malnutf1 = 1*(dents$malnutrition==1)
dents$malnutf2 = 1*(dents$malnutrition==2)
dents$malarf = as.factor(dents$malar)
dents$malarf1 = 1*(dents$malar==1)
dents$malarf2 = 1*(dents$malar==2)
dents$an1 = 1*(dents$anemia_corr==1)
dents$an2 = 1*(dents$anemia_corr==2)
dents$anf = as.factor(dents$anemia_corr)


# Make dataset with complete cases (very few missings: 4/2327, 0.2%)

dents$cc = 1 - 1*(is.na(dents$malar) | is.na(dents$anemia_corr))
dents_cc = dents[dents$cc==1,]


# Make extra variable to model continuous variables using RCS with 3 knots (only needed on CC dataset)

dents_cc$agercs = rcspline.eval(dents_cc$age_m_calc,nk=3)
dents_cc$durrcs = rcspline.eval(dents_cc$dur_fev,nk=3)
dents_cc$l2durrcs = rcspline.eval(dents_cc$l2durfev,nk=3)
dents_cc$resprcs = rcspline.eval(dents_cc$resp_man,nk=3)
dents_cc$hbrcs = rcspline.eval(dents_cc$hb_corr,nk=3)



#########################################
# DESCRIPTIVE STATISTICS FOR PREDICTORS #
#########################################

Desc(dents$age_m_calc)
Desc(dents$dur_fev)
Desc(dents$resp_man)
Desc(dents$hb_corr)

CrossTable(dents$agemcat)
CrossTable(dents$sex)
CrossTable(dents$durcat)
CrossTable(dents$malnutrition)
CrossTable(dents$tachypnea)
CrossTable(dents$distress)
CrossTable(dents$n_antibiotics)
CrossTable(dents$malar)
CrossTable(dents$season)
CrossTable(dents$mh_diarrh)
CrossTable(dents$hepaspl)
CrossTable(dents$anemia_corr)
dents$hb_corr[is.na(dents$anemia_corr)]



#######################################
# HISTOGRAMS OF CONTINUOUS PREDICTORS #
#######################################

hist(dents$age_m_calc,xlab="Age (months)",main=NULL,breaks=60,cex.lab=1.5,cex.axis=0.8)

hist(dents$dur_fev,xlab="Fever duration (days)",main=NULL,breaks=32,cex.lab=1.5,cex.axis=0.8)
hist(log(dents$dur_fev+1),xlab="Log(Fever duration (days+1))",main=NULL,breaks=32,cex.lab=1.5,cex.axis=0.8)

hist(dents$resp_man,xlab="Respiratory rate (breaths/min)",main=NULL,breaks=22,cex.lab=1.5,cex.axis=0.8)

hist(dents$hb_corr,xlab="Hemoglobin (g/dl)",main=NULL,breaks=38,cex.lab=1.5,cex.axis=0.8)



##########################################
# DESCRIPTIVE STATISTICS FOR THE OUTCOME #
##########################################

CrossTable(dents$def6) # negative vs {possible bsi, bsi, nts}
CrossTable(dents$bsi) # {negative, possible bsi} vs {bsi, nts}
CrossTable(dents$nts) # {negative, possible bsi, bsi} vs nts



###################################
# CORRELATIONS BETWEEN PREDICTORS #
###################################


# Spearman: used for two continuous/ordinal variables

round(cor(dents[,c("age_m_calc","season","dur_fev","mh_diarrh","malnutrition","resp_man","tachypnea","distress",
                   "hepaspl","malarf1","malarf2","hb_corr","anemia_corr","n_antibiotics","sex")],
          use="pairwise.complete.obs",method='spearman'),digits=2)


# Pearson: used to calculate point biserial correlation between continuous/ordinal and binary variables,
#          and to calculate phi correlation between two binary variables

round(cor(dents[,c("age_m_calc","season","dur_fev","mh_diarrh","malnutrition","resp_man","tachypnea","distress",
                   "hepaspl","malarf1","malarf2","hb_corr","anemia_corr","n_antibiotics","sex")],
          use="pairwise.complete.obs",method='pearson'),digits=2)



############################################
# SAMPLE SIZE CALCULATIONS FOR NTS VS REST #
############################################

pmsampsize(type = "b", rsquared = 0.08, parameters = 17, shrinkage = 0.9, prevalence = 0.13) 
# max_rsq 0.54, 15% of that is 0.08 (C approx 0.733); RESULT 1827



########################################
# UNIVARIABLE AUC'S FOR THE PREDICTORS #
########################################

round(auc.nonpara.mw(dents_cc$age_m_calc[dents_cc$nts==1],dents_cc$age_m_calc[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$season[dents_cc$nts==1],dents_cc$season[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$dur_fev[dents_cc$nts==1],dents_cc$dur_fev[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$mh_diarrh[dents_cc$nts==1],dents_cc$mh_diarrh[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$malnutrition[dents_cc$nts==1],dents_cc$malnutrition[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$resp_man[dents_cc$nts==1],dents_cc$resp_man[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$distress[dents_cc$nts==1],dents_cc$distress[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$hepaspl[dents_cc$nts==1],dents_cc$hepaspl[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$malarf1[dents_cc$nts==1],dents_cc$malarf1[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$malarf2[dents_cc$nts==1],dents_cc$malarf2[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$hb_corr[dents_cc$nts==1],dents_cc$hb_corr[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$anemia_corr[dents_cc$nts==1],dents_cc$anemia_corr[dents_cc$nts==0],method="pepe"), digits=2)
round(auc.nonpara.mw(dents_cc$n_antibiotics[dents_cc$nts==1],dents_cc$n_antibiotics[dents_cc$nts==0],method="pepe"), digits=2)



##################################################################################
# UNIVARIABLE DESCRIPTIVE ASSESSMENT OF FUNCTIONAL FORM FOR CONTINUOUS VARIABLES #
##################################################################################


# AGE

agel = lrm(nts ~ age_m_calc, data=dents_cc)
agercs = lrm(nts ~ rcs(age_m_calc,3), data=dents_cc)
agelp0 = predict(agel, newdata=c(1:59), se.fit=T)
agercsp0 = predict(agercs, newdata=c(1:59), se.fit=T)
agep = as.data.frame(cbind(c(1:59),expit(cbind(agercsp0[[1]],agercsp0[[1]]-1.96*agercsp0[[2]],agercsp0[[1]]+1.96*agercsp0[[2]])),
                            expit(agelp0[[1]])))

ggplot(agep, aes(x=V1)) + 
  geom_line(aes(y = V2), color = "black", lwd=2) + 
  geom_line(aes(y = V3), color="steelblue", linetype="twodash", lwd=1.5) + 
  geom_line(aes(y = V4), color="steelblue", linetype="twodash", lwd=1.5) +
  geom_line(aes(y = V5), color = "black", linetype="dotted", lwd=1) + 
  labs(x = "Age (months)", y="Probability NTS")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


# FEVER DURATION

durl = lrm(nts ~ dur_fev, data=dents_cc)
durl2 = lrm(nts ~ l2durfev, data=dents_cc)
durrcs = lrm(nts ~ rcs(l2durfev,3), data=dents_cc)
durlp0 = predict(durl, newdata=c(0:40), se.fit=T)
durl2p0 = predict(durl2, newdata=log2(c(1:41)), se.fit=T)
durrcsp0 = predict(durrcs, newdata=log2(c(1:41)), se.fit=T)
durp = as.data.frame(cbind(log(c(1:41)),expit(cbind(durrcsp0[[1]],durrcsp0[[1]]-1.96*durrcsp0[[2]],
                      durrcsp0[[1]]+1.96*durrcsp0[[2]])),c(0:40),expit(durlp0[[1]]),expit(durl2p0[[1]])))

ggplot(durp, aes(x=V5)) + 
  geom_line(aes(y = V2), color = "black", lwd=2) + 
  geom_line(aes(y = V3), color="steelblue", linetype="twodash", lwd=1.5) + 
  geom_line(aes(y = V4), color="steelblue", linetype="twodash", lwd=1.5) +
  geom_line(aes(y = V6), color = "black", linetype="dotted", lwd=1) + 
  geom_line(aes(y = V7), color = "black", linetype="dashed", lwd=1) + 
  labs(x = "Fever duration in days", y="Probability NTS")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


# RESPIRATORY RATE

resl = lrm(nts ~ resp_man, data=dents_cc)
resrcs = lrm(nts ~ rcs(resp_man,3), data=dents_cc)
reslp0 = predict(resl, newdata=c(6:108), se.fit=T)
resrcsp0 = predict(resrcs, newdata=c(6:108), se.fit=T)
resp = as.data.frame(cbind(c(6:108),expit(cbind(resrcsp0[[1]],resrcsp0[[1]]-1.96*resrcsp0[[2]],
                      resrcsp0[[1]]+1.96*resrcsp0[[2]])),expit(reslp0[[1]])))

ggplot(resp, aes(x=V1)) + 
  geom_line(aes(y = V2), color = "black", lwd=2) + 
  geom_line(aes(y = V3), color="steelblue", linetype="twodash", lwd=1.5) + 
  geom_line(aes(y = V4), color="steelblue", linetype="twodash", lwd=1.5) +
  geom_line(aes(y = V5), color = "black", linetype="dotted", lwd=1) + 
  labs(x = "Respiratory rate (breaths/min)", y="Probability NTS")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


# HEMOGLOBIN

hbl = lrm(nts ~ hb_corr, data=dents_cc)
hbrcs = lrm(nts ~ rcs(hb_corr,3), data=dents_cc)
hblp0 = predict(hbl, newdata=c(14:198)/10, se.fit=T)
hbrcsp0 = predict(hbrcs, newdata=c(14:198)/10, se.fit=T)
hbp = as.data.frame(cbind(c(14:198)/10,expit(cbind(hbrcsp0[[1]],hbrcsp0[[1]]-1.96*hbrcsp0[[2]],
                                                     hbrcsp0[[1]]+1.96*hbrcsp0[[2]])),expit(hblp0[[1]])))

ggplot(hbp, aes(x=V1)) + 
  geom_line(aes(y = V2), color = "black", lwd=2) + 
  geom_line(aes(y = V3), color="steelblue", linetype="twodash", lwd=1.5) + 
  geom_line(aes(y = V4), color="steelblue", linetype="twodash", lwd=1.5) +
  geom_line(aes(y = V5), color = "black", linetype="dotted", lwd=1) + 
  labs(x = "Hemoglobin (g/dl)", y="Probability NTS")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))



##########################
# MULTIVARIABLE MODELING #
##########################


# FIRST: BACKWARD ELIMINATION (ALPHA 1%) FOR RCS TERMS
# THEN: BACKWARD ELIMINATION ON WHAT IS LEFT (ALPHA 1%) - THIS MAY ELIMINATE LINEAR TERM BUT NOT RCS TERM FOR A CONTINUOUS VARIABLE
# FINALLY: FIT A FIRTH MODEL ON THE REMAINING TERMS
mod1 = lrm(nts ~ age_m_calc + season + ldurfev + mh_diarrh + malnutf + resp_man + distress + hepaspl + malarf + 
             hb_corr + n_antibiotics + agercs + l2durrcs + resprcs + hbrcs, data=dents_cc)
mod2 = fastbw(mod1, rule="p", sls = 0.01, type="individual", force=c(1:13))
mod3 <- lrm(nts~., data=dents_cc[,c("nts",mod2$names.kept)])
mod4 = fastbw(mod3, rule="p", sls = 0.01, type="individual")
mod5 = logistf(nts ~ ., flic=TRUE, data=dents_cc[,c("nts",mod4$names.kept)])

# Odds ratios (also get coeffs and ORs for age by 12 months and respiratory rate per 10 breaths/min)
exp(cbind(mod5$coefficients,mod5$coefficients-1.96*sqrt(diag(mod5$var)),mod5$coefficients+1.96*sqrt(diag(mod5$var))))
dents_cc$age12 = dents_cc$age_m_calc/12
dents_cc$resp10 = dents_cc$resp_man/10
mod5b = logistf(nts ~ age12 + season + l2durfev + mh_diarrh + malnutf + resp10 + hepaspl + malarf + hb_corr, 
                 flic=TRUE, data=dents_cc)
or=exp(cbind(mod5b$coefficients,mod5b$coefficients-1.96*sqrt(diag(mod5b$var)),mod5b$coefficients+1.96*sqrt(diag(mod5b$var))))


# GET ESTIMATED RISKS FOR FINAL MODEL

estrisks = predict(mod5, type="response")
estrisks = as.data.frame(cbind(estrisks,dents_cc$nts))


# APPARENT C STATISTIC, CALIBRATION PLOT, NET BENEFIT, Sensitivity and Specificity

#appc = Cstat(estrisks$estrisks,dents_cc$nts)
appc = auc.nonpara.mw(estrisks$estrisks[dents_cc$nts==1],estrisks$estrisks[dents_cc$nts==0],method="pepe")
val.prob.ci.2(estrisks$estrisks,y=dents_cc$nts,smooth="loess",dostats = T) # use loess fit, another option for flexible (nonlinear) modeling
appdca <- decision_curve(V2~estrisks,
                               data = estrisks,
                               fitted.risk = TRUE,
                               thresholds = seq(0, 1, by = .05),
                               confidence.intervals = NA) 
plot_decision_curve( list(appdca), 
                     curve.names = c("DENTS model"),
                     col = c("blue"), 
                     ylim = c(-.005, 0.125), #set ylim
                     xlim = c(0,1),
                     lty = c(1), confidence.intervals = F,
                     standardize = FALSE, #plot Net benefit instead of standardized net benefit
                     legend.position = "topright",xlab="Risk threshold",cost.benefit.xlab="Harm to benefit ratio") 

appsesp = round(matrix(data=c(binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.2),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.25),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.3),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.35),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.4),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.45),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==1]>=0.5),sum(dents_cc$nts==1),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.2),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.25),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.3),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.35),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.4),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.45),sum(dents_cc$nts==0),method="wilson"),
                    binconf(sum(estrisks$estrisks[dents_cc$nts==0]<0.5),sum(dents_cc$nts==0),method="wilson")),
                    nrow=14,ncol=3,byrow=T),digits=2) # Apparent sensitivity and specificity

appnb = appdca[[1]][5:11,c(1,6,7)]



# VIOLIN PLOTS OF ESTIMATED RISKS BY OUTCOME

estrisks$ntsc[estrisks$V2==0] = "No NTS"
estrisks$ntsc[estrisks$V2==1] = "NTS"
p <- ggplot(estrisks, aes(x=ntsc, y=estrisks)) + 
  geom_violin() + theme(legend.position="none") + ylim(0,1) + 
  labs(x = "NTS status", y="Estimated risk of NTS") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
p



# BOOTSTRAPPING FOR OPTIMISM CORRECTION, SELECTION FREQUENCY, AND RISK INSTABILITY

nboot=200
simres = matrix(data = NA, nrow = nboot, ncol = 24)
dents_boot = dents_cc[,c("nts", "age_m_calc", "season", "l2durfev", "mh_diarrh", "malnutf", "malnutf1", "malnutf2", "resp_man", 
                         "distress", "hepaspl", "malarf", "malarf1", "malarf2", "hb_corr", "n_antibiotics", "agercs", "l2durrcs", 
                         "resprcs", "hbrcs")]
betaboot = as.data.frame(matrix(data = NA, nrow = 0, ncol = 17))
colnames(betaboot) = c("age_m_calc", "season", "l2durfev", "mh_diarrh", "malnutf1", "malnutf2", "resp_man", 
                       "distress", "hepaspl", "malarf1", "malarf2", "hb_corr", "n_antibiotics", "agercs", "l2durrcs", 
                       "resprcs", "hbrcs")
riskboot = as.data.frame(matrix(data = NA, nrow = 0, ncol = 2323))

set.seed(1234)
for (bootnr in 1:nboot){
  
  train_data <- dents_boot[sample(row.names(dents_boot),replace=TRUE),]
  
  # do the modeling
  mod1 = lrm(nts ~ age_m_calc + season + l2durfev + mh_diarrh + malnutf + resp_man + distress + hepaspl + malarf + hb_corr + 
                   n_antibiotics + agercs + l2durrcs + resprcs + hbrcs, data=train_data)
  mod2 = fastbw(mod1, rule="p", sls = 0.01, type="individual", force=c(1:13))
  mod3 <- lrm(nts~., data=train_data[,c("nts",mod2$names.kept)])
  mod4 = fastbw(mod3, rule="p", sls = 0.01, type="individual")
  mod5 = logistf(nts ~ ., flic=TRUE, data=train_data[,c("nts",mod4$names.kept)])
  betaboot = rbind.fill(betaboot,as.data.frame(t(mod5$coefficients[-1])))
  
  # predict the values on the bootstrap data
  predictboot <- predict(mod5)
  predictboot_risk <- predict(mod5, type="response")

  # predict the values on the original, unresampled data
  predictorig <- predict(mod5, newdata = dents_boot[,c(mod5$terms[-1])])
  predictorig_risk <- predict(mod5, newdata = dents_boot[,c(mod5$terms[-1])],type="response")
  riskboot = rbind.fill(riskboot,as.data.frame(t(predictorig_risk)))

  # return a vector of summary results
  simres[bootnr,1:24] <- c(
    Cstat(predictboot,train_data$nts), # AUC of bootstrap model on bootstrap dataset
    Cstat(predictorig,dents_boot$nts), # AUC of bootstrap model on original dataset
    glm(nts ~ predictorig,  family=binomial(link='logit'), data=dents_boot)$coefficients[2], # calibration slope
    # Optimism in NB at relevant thresholds:
    ((sum(predictboot_risk[train_data$nts==1]>=0.2) - (0.2/(1-0.2))*sum(predictboot_risk[train_data$nts==0]>=0.2)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.2) - (0.2/(1-0.2))*sum(predictorig_risk[dents_boot$nts==0]>=0.2)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.25) - (0.25/(1-0.25))*sum(predictboot_risk[train_data$nts==0]>=0.25)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.25) - (0.25/(1-0.25))*sum(predictorig_risk[dents_boot$nts==0]>=0.25)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.3) - (0.3/(1-0.3))*sum(predictboot_risk[train_data$nts==0]>=0.3)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.3) - (0.3/(1-0.3))*sum(predictorig_risk[dents_boot$nts==0]>=0.3)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.35) - (0.35/(1-0.35))*sum(predictboot_risk[train_data$nts==0]>=0.35)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.35) - (0.35/(1-0.35))*sum(predictorig_risk[dents_boot$nts==0]>=0.35)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.4) - (0.4/(1-0.4))*sum(predictboot_risk[train_data$nts==0]>=0.4)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.4) - (0.4/(1-0.4))*sum(predictorig_risk[dents_boot$nts==0]>=0.4)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.45) - (0.45/(1-0.45))*sum(predictboot_risk[train_data$nts==0]>=0.45)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.45) - (0.45/(1-0.45))*sum(predictorig_risk[dents_boot$nts==0]>=0.45)))/dim(dents_boot)[1],
    ((sum(predictboot_risk[train_data$nts==1]>=0.5) - (0.5/(1-0.5))*sum(predictboot_risk[train_data$nts==0]>=0.5)) - 
       (sum(predictorig_risk[dents_boot$nts==1]>=0.5) - (0.5/(1-0.5))*sum(predictorig_risk[dents_boot$nts==0]>=0.5)))/dim(dents_boot)[1],
    # Optimism in sensitivity at relevant thresholds:
    sum(predictboot_risk[train_data$nts==1]>=0.2)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.2)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.25)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.25)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.3)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.3)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.35)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.35)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.4)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.4)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.45)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.45)/sum(dents_boot$nts==1),
    sum(predictboot_risk[train_data$nts==1]>=0.5)/sum(train_data$nts==1)-sum(predictorig_risk[dents_boot$nts==1]>=0.5)/sum(dents_boot$nts==1),
    # Optimism in specificity at relevant thresholds:
    sum(predictboot_risk[train_data$nts==0]<0.2)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.2)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.25)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.25)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.3)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.3)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.35)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.35)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.4)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.4)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.45)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.45)/sum(dents_boot$nts==0),
    sum(predictboot_risk[train_data$nts==0]<0.5)/sum(train_data$nts==0)-sum(predictorig_risk[dents_boot$nts==0]<0.5)/sum(dents_boot$nts==0)
  ) # Net Benefit for thresholds (0.2:0.5) by 0.05

}

bootresults = list(simres,betaboot,t(riskboot))
save(bootresults, file = "C://Ben//Development and Regeneration//Bieke Tack//Data March2023//bootresults.RData") 
load("C://Ben//Development and Regeneration//Bieke Tack//Data March2023//bootresults.RData")


# OPTIMISM CORRECTED PERFORMANCE

# C statistic
optc = mean(bootresults[[1]][,1]) - mean(bootresults[[1]][,2])
ccorr = appc - optc

# Calibration slope
mean(bootresults[[1]][,3])

# Net Benefit, Sensitivity, Specificity
NBopt = colMeans(bootresults[[1]][,4:10])
Sensopt = colMeans(bootresults[[1]][,11:17])
Specopt = colMeans(bootresults[[1]][,18:24])
NBcorr = appnb[,2] - NBopt
Senscorr = round(cbind(appsesp[1:7,1] - Sensopt, appsesp[1:7,2] - Sensopt, appsesp[1:7,3] - Sensopt), digits = 2)
Speccorr = round(cbind(appsesp[8:14,1] - Specopt, appsesp[8:14,2] - Specopt, appsesp[8:14,3] - Specopt), digits = 2)

# Decision curve
plot(100*appdca[[1]][22:42,1], appdca[[1]][22:42,6], type="l",lty=2,col="gray",xlim=c(0,100),ylim=c(-0.01,0.15),
     ylab="Net Benefit", xlab="Decision threshold (% risk of NTS)")
polygon(c(20,50,50,20), c(-0.0165,-0.0165,0.156,0.156),
        col = "lightgray", border=NA)
lines(c(20,25,30,35,40,45,50), NBcorr, lwd=2,col="black")
abline(h=0,col="gray")
lines(100*appdca[[1]][1:21,1], appdca[[1]][1:21,6])
legend("topright", inset=.02, legend=c("DeNTS model (uncorrected)", "DeNTS model (optimism-corrected)", 
                                       "Treat all with modified antibiotics", "Treat none (SoC antibiotics for all)"), 
       col=c("black","black","gray","gray"), lty=c(1,1,2,1), lwd=c(1,2,1,1), cex=0.8)


# SELECTION FREQUENCY

selfreq = round(colMeans(1*!is.na(bootresults[[2]])), digits = 2)


# RISK INSTABILITY PLOT

dim(bootresults[[3]])
meanrisk = rowMeans(bootresults[[3]])
meanriskr = rank(meanrisk)
plot(meanriskr,bootresults[[3]][,1],type="p", pch=".", cex=0.1, xlab="Patient", ylab="Estimated risk of NTS", ylim=c(0,1))
for (i in 2:200){
  points(meanriskr,bootresults[[3]][,i], pch=".", cex=0.1)
}
points(meanriskr, meanrisk, pch=".", col="red", cex=1)



######################################################
# MULTIVARIABLE MODELING WITH CATEGORIZED PREDICTORS #
######################################################

# FIT MODEL WITH BINARY VERSION FOR RESPIRATORY RATE (TACHYPNEA) AND HB_CORR (ANEMIA_CORR)

modc = lrm(nts ~ age_m_calc + season + l2durfev + mh_diarrh + malnutf + tachypnea + hepaspl + malarf + anf, data=dents_cc)
modc$coefficients

# MAKE SCORE SYSTEM USING METHOD FROM SULLIVAN ET AL (STAT MED 2004)

points = round(c(0, modc$coefficients[2]*12, modc$coefficients[2]*24, modc$coefficients[2]*36, modc$coefficients[2]*48,
           0, modc$coefficients[4]*1, modc$coefficients[4]*1.58, modc$coefficients[4]*2, modc$coefficients[4]*3,
           modc$coefficients[4]*4, 0, modc$coefficients[3], 0, modc$coefficients[5], 0, modc$coefficients[6], 0,
           modc$coefficients[7], 0, modc$coefficients[8], 0, modc$coefficients[9], 0, modc$coefficients[10], 0, 
           modc$coefficients[11], 0, modc$coefficients[12], 0, modc$coefficients[13])/abs(modc$coefficients[2]*12),digits = 0)

dents_cc$tpts = points[2]*(dents_cc$age_m_calc>=12 & dents_cc$age_m_calc<=23) + 
                points[3]*(dents_cc$age_m_calc>=24 & dents_cc$age_m_calc<=35) +
                points[4]*(dents_cc$age_m_calc>=36 & dents_cc$age_m_calc<=47) +
                points[5]*(dents_cc$age_m_calc>=48) + 
                points[7]*(dents_cc$dur_fev==1) + 
                points[8]*(dents_cc$dur_fev==2) + 
                points[9]*(dents_cc$dur_fev==3) +
                points[10]*(dents_cc$dur_fev>=4 & dents_cc$dur_fev<=10) +
                points[11]*(dents_cc$dur_fev>10) +
                points[13]*(dents_cc$season==1) + 
                points[15]*(dents_cc$mh_diarrh==1) + 
                points[17]*(dents_cc$malnutrition==1) + 
                points[19]*(dents_cc$malnutrition==2) +
                points[21]*(dents_cc$tachypnea==1) +
                points[23]*(dents_cc$hepaspl==1) + 
                points[25]*(dents_cc$malarf1==1) +
                points[27]*(dents_cc$malarf2==1) +
                points[29]*(dents_cc$an1==1) +
                points[31]*(dents_cc$an2==1)
dents_cc$tpts2 = modc$coefficients[2]*12*dents_cc$tpts + modc$coefficients[1] + modc$coefficients[2]*6


# GET ESTIMATED RISKS FOR SCORES

# risks for patients in dataset
ptsmod = lrm(nts ~ tpts2, data=dents_cc)
ptsrisk = expit(predict(ptsmod))
scoredata = as.data.frame(cbind(ptsrisk,dents_cc$nts))
colnames(scoredata) = c("ptsrisk","nts")

# risks for each score between -8 and 22 (to get score table)
scoreref = as.data.frame(c(-8:22))
colnames(scoreref)=c("tpts")
scoreref$tpts2 = modc$coefficients[2]*12*scoreref$tpts + modc$coefficients[1] + modc$coefficients[2]*6
scoreref$ptsrisk = expit(predict(ptsmod,newdata=scoreref))


# APPARENT PERFORMANCE

# C statistic
appc2 = auc.nonpara.mw(ptsrisk[dents_cc$nts==1],ptsrisk[dents_cc$nts==0],method="pepe")

# Calibration
val.prob.ci.2(ptsrisk,y=dents_cc$nts,smooth="loess",dostats = T) # use loess fit, another option for flexible (nonlinear) modeling

# Net Benefit
appdca2 <- decision_curve(nts~ptsrisk,
                          data=scoredata,
                         fitted.risk = TRUE,
                         thresholds = seq(0, 1, by = .05),
                         confidence.intervals = NA) 
plot_decision_curve( list(appdca2), 
                     curve.names = c("DENTS score"),
                     col = c("blue"), 
                     ylim = c(-.005, 0.125), #set ylim
                     xlim = c(0,1),
                     lty = c(1), confidence.intervals = F,
                     standardize = FALSE, #plot Net benefit instead of standardized net benefit
                     legend.position = "topright",xlab="Risk threshold",cost.benefit.xlab="Harm to benefit ratio") 

appnb2 = appdca2[[1]][5:11,c(1,6,7)]

# Sensitivity and specificity
appsesp2 = round(matrix(data=c(binconf(sum(ptsrisk[dents_cc$nts==1]>=0.2),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.25),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.3),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.35),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.4),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.45),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==1]>=0.5),sum(dents_cc$nts==1),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.2),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.25),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.3),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.35),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.4),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.45),sum(dents_cc$nts==0),method="wilson"),
                              binconf(sum(ptsrisk[dents_cc$nts==0]<0.5),sum(dents_cc$nts==0),method="wilson")),
                        nrow=14,ncol=3,byrow=T),digits=2) # Apparent sensitivity and specificity


# SCATTER PLOT OF RISKS FROM MODEL VS SCORE SYSTEM

plot(estrisks$estrisks, ptsrisk, type="p", xlab="Estimated risk from risk model", ylab="Estimated risk from score system", 
     xlim=c(0,1), ylim=c(0,1), cex.lab=1.5)
abline(a=0, b=1, col="red")


# BOOTSTRAP CORRECTION FOR SCORE SYSTEM (NOTE: THIS CONDITIONS ON VARIABLE SELECTION SO NO FULL OPTIMISM CORRECTION!)

nboot=200
simres2 = matrix(data = NA, nrow = nboot, ncol = 24)
dents_boot2 = dents_cc[,c("nts", "age_m_calc", "season", "l2durfev", "dur_fev", "mh_diarrh", "malnutf", "malnutf1", "malnutf2", "tachypnea", 
                         "hepaspl", "malarf", "malarf1", "malarf2", "anf", "an1", "an2")]

set.seed(1234)
for (bootnr in 1:nboot){
  
  train_data <- dents_boot2[sample(row.names(dents_boot2),replace=TRUE),]
  boots = as.data.frame(matrix(data = NA, nrow = 0, ncol = 2323))
  origs = as.data.frame(matrix(data = NA, nrow = 0, ncol = 2323))
  
  # do the modeling
  modc = lrm(nts ~ age_m_calc + season + l2durfev + mh_diarrh + malnutf + tachypnea + hepaspl + malarf + anf, data=train_data)
  points = round(c(0, modc$coefficients[2]*12, modc$coefficients[2]*24, modc$coefficients[2]*36, modc$coefficients[2]*48,
                   0, modc$coefficients[4]*1, modc$coefficients[4]*1.58, modc$coefficients[4]*2, modc$coefficients[4]*3,
                   modc$coefficients[4]*4, 0, modc$coefficients[3], 0, modc$coefficients[5], 0, modc$coefficients[6], 0,
                   modc$coefficients[7], 0, modc$coefficients[8], 0, modc$coefficients[9], 0, modc$coefficients[10], 0, 
                   modc$coefficients[11], 0, modc$coefficients[12], 0, modc$coefficients[13])/abs(modc$coefficients[2]*12),digits = 0)

  # predict the values on the bootstrap data
  boots = as.data.frame(t(rbind.fill(boots,as.data.frame(t(points[2]*(train_data$age_m_calc>=12 & train_data$age_m_calc<=23) + 
                                             points[3]*(train_data$age_m_calc>=24 & train_data$age_m_calc<=35) + 
                                             points[4]*(train_data$age_m_calc>=36 & train_data$age_m_calc<=47) +
                                             points[5]*(train_data$age_m_calc>=48) + points[7]*(train_data$dur_fev==1) + 
                                             points[8]*(train_data$dur_fev==2) + points[9]*(train_data$dur_fev==3) + 
                                             points[10]*(train_data$dur_fev>=4 & train_data$dur_fev<=10) +
                                             points[11]*(train_data$dur_fev>10) + points[13]*(train_data$season==1) + 
                                             points[15]*(train_data$mh_diarrh==1) + 
                                             points[17]*(train_data$malnutf1==1) + points[19]*(train_data$malnutf2==1) + 
                                             points[21]*(train_data$tachypnea==1) +
                                             points[23]*(train_data$hepaspl==1) + points[25]*(train_data$malarf1==1) + 
                                             points[27]*(train_data$malarf2==1) +
                                             points[29]*(train_data$an1==1) + points[31]*(train_data$an2==1))))))
  colnames(boots) = c("tpts_b")
  boots$tpts2 = modc$coefficients[2]*12*boots$tpts_b + modc$coefficients[1] + modc$coefficients[2]*6
  boots$nts = train_data$nts
  ptsmod = lrm(nts ~ tpts2, data=boots)
  ptsrisk_b = expit(predict(ptsmod))
  
  # predict the values on the original, unresampled data
  origs = as.data.frame(t(rbind.fill(origs,as.data.frame(t(points[2]*(dents_boot2$age_m_calc>=12 & dents_boot2$age_m_calc<=23) + 
                                                             points[3]*(dents_boot2$age_m_calc>=24 & dents_boot2$age_m_calc<=35) + 
                                                             points[4]*(dents_boot2$age_m_calc>=36 & dents_boot2$age_m_calc<=47) +
                                                             points[5]*(dents_boot2$age_m_calc>=48) + points[7]*(dents_boot2$dur_fev==1) + 
                                                             points[8]*(dents_boot2$dur_fev==2) + points[9]*(dents_boot2$dur_fev==3) + 
                                                             points[10]*(dents_boot2$dur_fev>=4 & dents_boot2$dur_fev<=10) +
                                                             points[11]*(dents_boot2$dur_fev>10) + points[13]*(dents_boot2$season==1) + 
                                                             points[15]*(dents_boot2$mh_diarrh==1) + points[17]*(dents_boot2$malnutf1==1) + 
                                                             points[19]*(dents_boot2$malnutf2==1) + points[21]*(dents_boot2$tachypnea==1) +
                                                             points[23]*(dents_boot2$hepaspl==1) + points[25]*(dents_boot2$malarf1==1) + 
                                                             points[27]*(dents_boot2$malarf2==1) + points[29]*(dents_boot2$an1==1) + 
                                                             points[31]*(dents_boot2$an2==1))))))
  colnames(origs) = c("tpts_o")

  origs$tpts2 = modc$coefficients[2]*12*origs$tpts_o + modc$coefficients[1] + modc$coefficients[2]*6
  origs$nts = dents_boot2$nts
  ptsrisk_o = expit(predict(ptsmod, newdata = origs))
  
  # return a vector of summary results
  simres2[bootnr,1:24] <- c(
    Cstat(ptsrisk_b,train_data$nts), # AUC of bootstrap model on bootstrap dataset
    Cstat(ptsrisk_o,dents_boot2$nts), # AUC of bootstrap model on original dataset
    glm(nts ~ logit(ptsrisk_o),  family=binomial(link='logit'), data=dents_boot2)$coefficients[2], # calibration slope
    ((sum(ptsrisk_b[train_data$nts==1]>=0.2) - (0.2/(1-0.2))*sum(ptsrisk_b[train_data$nts==0]>=0.2)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.2) - (0.2/(1-0.2))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.2)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.25) - (0.25/(1-0.25))*sum(ptsrisk_b[train_data$nts==0]>=0.25)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.25) - (0.25/(1-0.25))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.25)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.3) - (0.3/(1-0.3))*sum(ptsrisk_b[train_data$nts==0]>=0.3)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.3) - (0.3/(1-0.3))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.3)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.35) - (0.35/(1-0.35))*sum(ptsrisk_b[train_data$nts==0]>=0.35)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.35) - (0.35/(1-0.35))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.35)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.4) - (0.4/(1-0.4))*sum(ptsrisk_b[train_data$nts==0]>=0.4)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.4) - (0.4/(1-0.4))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.4)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.45) - (0.45/(1-0.45))*sum(ptsrisk_b[train_data$nts==0]>=0.45)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.45) - (0.45/(1-0.45))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.45)))/dim(dents_boot2)[1],
    ((sum(ptsrisk_b[train_data$nts==1]>=0.5) - (0.5/(1-0.5))*sum(ptsrisk_b[train_data$nts==0]>=0.5)) - 
       (sum(ptsrisk_o[dents_boot2$nts==1]>=0.5) - (0.5/(1-0.5))*sum(ptsrisk_o[dents_boot2$nts==0]>=0.5)))/dim(dents_boot2)[1],
    sum(ptsrisk_b[train_data$nts==1]>=0.2)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.2)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.25)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.25)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.3)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.3)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.35)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.35)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.4)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.4)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.45)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.45)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==1]>=0.5)/sum(train_data$nts==1)-sum(ptsrisk_o[dents_boot2$nts==1]>=0.5)/sum(dents_boot2$nts==1),
    sum(ptsrisk_b[train_data$nts==0]<0.2)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.2)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.25)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.25)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.3)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.3)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.35)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.35)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.4)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.4)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.45)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.45)/sum(dents_boot2$nts==0),
    sum(ptsrisk_b[train_data$nts==0]<0.5)/sum(train_data$nts==0)-sum(ptsrisk_o[dents_boot2$nts==0]<0.5)/sum(dents_boot2$nts==0)
  ) # Net Benefit for thresholds (0.2:0.5) by 0.05
  
}

bootresults2 = list(simres2)
save(bootresults2, file = "C://Ben//Development and Regeneration//Bieke Tack//Data March2023//bootresults2.RData") 
load("C://Ben//Development and Regeneration//Bieke Tack//Data March2023//bootresults2.RData")


# OPTIMISM CORRECTED PERFORMANCE FOR SCORE SYSTEM
# I add the optimism calculation of the optimal model and that of the score system to approximate the full optimism
# of variable selection followed by score system derivation
# This is an ad hoc approach

# C statistic
summary(bootresults2[[1]])
optc2 = mean(bootresults2[[1]][,1]) - mean(bootresults2[[1]][,2])
optc12 = optc + optc2
ccorr2 = appc2 - optc12

# Calibration slope
1 - ((1 - mean(bootresults[[1]][,3])) + (1 - mean(bootresults2[[1]][,3])))
#round(mean(bootresults[[1]][,3])*mean(bootresults2[[1]][,3]), digits = 2)

# Net Benefit, Sensitivity, Specificity
NBopt2 = colMeans(bootresults2[[1]][,4:10])
Sensopt2 = colMeans(bootresults2[[1]][,11:17])
Specopt2 = colMeans(bootresults2[[1]][,18:24])
NBopt12 = NBopt + NBopt2
Sensopt12 = Sensopt + Sensopt2
Specopt12 = Specopt + Specopt2
NBcorr2 = appnb2[,2] - NBopt12
Senscorr2 = round(cbind(appsesp2[1:7,1] - Sensopt12, appsesp2[1:7,2] - Sensopt12, appsesp2[1:7,3] - Sensopt12), digits = 2)
Speccorr2 = round(cbind(appsesp2[8:14,1] - Specopt12, appsesp2[8:14,2] - Specopt12, appsesp2[8:14,3] - Specopt12), digits = 2)

# Decision curve
plot(100*appdca[[1]][22:42,1], appdca[[1]][22:42,6], type="l",lty=2,col="gray",xlim=c(0,100),ylim=c(-0.01,0.15),
     ylab="Net Benefit", xlab="Decision threshold (% risk of NTS)")
polygon(c(20,50,50,20), c(-0.0165,-0.0165,0.156,0.156),
        col = "lightgray", border=NA)
lines(c(20,25,30,35,40,45,50), NBcorr, lwd=2,col="black")
lines(c(20,25,30,35,40,45,50), NBcorr2, lwd=2, lty=2, col="black")
abline(h=0,col="gray")
lines(100*appdca[[1]][1:21,1], appdca[[1]][1:21,6])
legend("topright", inset=.02, legend=c("DeNTS model (uncorrected)", "DeNTS model (optimism-corrected)","DeNTS score (optimism-corrected)", 
                                       "Treat all with modified antibiotics", "Treat none (SoC antibiotics for all)"), 
       col=c("black","black","black","gray","gray"), lty=c(1,1,2,2,1), lwd=c(1,2,2,1,1), cex=0.8)
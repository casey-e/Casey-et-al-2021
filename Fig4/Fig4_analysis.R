rm(list=ls())

df <- read.csv("https://raw.githubusercontent.com//casey-e//Casey-et-al-2022//main//Fig4//Fig4_dataframe.csv")
df$Cohort<-factor(df$Cohort,levels=c("1","2"))
df$Group<-factor(df$Group,levels=c("Ctrl","CeADrd2KO"))
View(df)



#### Perform individual univariate tests and verify asumptions ####

## OF: Time in center ##
library(nlme)
OF_Center<-gls(OF_time_in_center_seconds~Group, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(OF_Center, type='pearson')
f<-fitted(OF_Center)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(OF_Center, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(OF_time_in_center_seconds~Group, data=df)

#Evaluate significance
summary(OF_Center)



## DLBT: time in light ##
library(nlme)
DLBT_TimeInLight<-gls(DLBT_time_in_light_seconds~Group, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(DLBT_TimeInLight, type='pearson')
f<-fitted(DLBT_TimeInLight)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(DLBT_TimeInLight, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(DLBT_time_in_light_seconds~Group, data=df)

#Evaluate significance
summary(DLBT_TimeInLight)




## DLBT: latency ##
library(nlme)
DLBT_latency<-gls(DLBT_latency_seconds~Group, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(DLBT_latency, type='pearson')
f<-fitted(DLBT_latency)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(DLBT_latency, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(DLBT_latency_seconds~Group, data=df)

#Evaluate significance
summary(DLBT_latency)



## DLBT: number_of_entries ##
library(lme4)
DLBT_entries<-glm(DLBT_number_of_entries ~ Group, df, family=poisson)

#Check assumptions
#Evaluate residuals 
e<-residuals(DLBT_entries, type='pearson')
f<-fitted(DLBT_entries)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

disp<-sum(e^2)/df.residual(DLBT_entries)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}

#Evaluate significance with Wald's test
summary(DLBT_entries)


## EPM: percentage of time on open arms ##
library(nlme)
EPM_TimeOpen<-gls(EPM_percentage_time_on_open~Group, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_TimeOpen, type='pearson')
f<-fitted(EPM_TimeOpen)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(EPM_TimeOpen, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(EPM_percentage_time_on_open~Group, data=df)

#Residuals doensn't show normal distribution

#Evaluate 2 ways model including Cohort as factor
EPM_TimeOpen_TwoWayANOVA<-gls(EPM_percentage_time_on_open~Group*Cohort, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_TimeOpen_TwoWayANOVA, type='pearson')
f<-fitted(EPM_TimeOpen_TwoWayANOVA)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(EPM_TimeOpen_TwoWayANOVA, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(EPM_percentage_time_on_open~Group*Cohort, data=df)

#Residuals doesn't have normal distribution. Besides Levene's test is not significant, residuals plot indicate heterocedasticity between cohorts
#Evaluate 2 ways model including Cohort as factor, with Variance modeling using VarIdent for cohort
EPM_TimeOpen_TwoWayANOVAVarIdent<-gls(EPM_percentage_time_on_open~Group*Cohort, weights=varIdent(form=~1|Cohort),data=df)
#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_TimeOpen_TwoWayANOVAVarIdent, type='pearson')
f<-fitted(EPM_TimeOpen_TwoWayANOVAVarIdent)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(EPM_TimeOpen_TwoWayANOVAVarIdent, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality. Levene test doesn't apply because variance was modeled
shapiro.test(e)

#Evaluate if modeling variances is appropiate by evaluating Akaike information criterion of model with and without VarIdent
aic=AIC(EPM_TimeOpen_TwoWayANOVA, EPM_TimeOpen_TwoWayANOVAVarIdent)
if (aic[1,2]-aic[2,2]>2){
  print(paste('Akaike of VarIdent model is lower than the Akaike of the model without modeling by ',as.character(round((aic[1,2]-aic[2,2]), 2))))
}else{
  print('The mode with variance modeling is not selected')
}

#Evaluate significance
anova(EPM_TimeOpen_TwoWayANOVAVarIdent)


#Relativize percentage_time_on_open to average of each cohort, and perform a t-test with Group as the only factor
#Create 'normalized_time_on_open' variable
for(i in levels(df$Cohort)) {
  df[df$Cohort==i,"EPM_normalized_time_on_open"]=df[df$Cohort==i,"EPM_percentage_time_on_open"]/mean(df[df$Cohort==i,"EPM_percentage_time_on_open"])
}


# EPM: noralized time on open arms, t-test
library(nlme)
EPM_NormTimeOpen<-gls(EPM_normalized_time_on_open~Group, data=df)

#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_NormTimeOpen, type='pearson')
f<-fitted(EPM_NormTimeOpen)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(EPM_NormTimeOpen, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(EPM_normalized_time_on_open~Group, data=df)

#Evaluate significance
summary(EPM_NormTimeOpen)


## EPM: proportion_of_entries_to_open_arms ##

#Given results in percentage of time in open, first try 2 ways model
library(lme4)
EPM_entries_TwoWay<-glm(EPM_proportion_entries_to_open ~ Group*Cohort, df, family=binomial, weights = EPM_total_entries)

#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_entries_TwoWay, type='pearson')
f<-fitted(EPM_entries_TwoWay)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

disp<-sum(e^2)/df.residual(EPM_entries_TwoWay)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


EPM_entries_TwoWay_QuasiBin<-glm(EPM_proportion_entries_to_open ~ Group*Cohort, df, family=quasibinomial, weights = EPM_total_entries)


#Evaluate significance with Wald's test
summary(EPM_entries_TwoWay_QuasiBin)


####      MANOVA  #####
df_M<-subset(df, select = -c(1, 3, 9, 10))
View(df_M)
attach(df_M)

Ctrl<- subset(df_M[,2:7], df_M$Group == "Ctrl")
CeADrd2KO<- subset(df_M[,2:7], df_M$Group == "CeADrd2KO")


## Calculate centroids
(centroid_Ctrl<-colMeans(Ctrl))
(centroid_CeADrd2KO<-colMeans(CeADrd2KO))
(total_centroid<-colMeans(df_M[,-c(1)]))
(centroids<-rbind(centroid_Ctrl,centroid_CeADrd2KO, total_centroid))
centroids<-data.frame(centroids)


#Mahalanobis distance to evaluate multivariate outliers
#For Ctrl group
s_Ctrl<-cov(Ctrl)  #co-variances matrix
D2_Ctrl<-mahalanobis(Ctrl,centroid_Ctrl,s_Ctrl)
(MV_out1<-which(D2_Ctrl>qchisq(.95, df=5)))

#For CeADrd2KO group
s_CeADrd2KO<-cov(CeADrd2KO)   #co-variances matrix
D2_CeADrd2KO<-mahalanobis(CeADrd2KO,centroid_CeADrd2KO,s_CeADrd2KO)
(MV_out2<-which(D2_CeADrd2KO>qchisq(.95, df=5)))


#None multivariate outliers


##Asumptions
#Multivariate normality
library(MVN)
mvn(
  CeADrd2KO,
  subset = NULL,
  mvnTest = "royston",
  covariance = TRUE,
  tol = 1e-25,
  alpha = 0.5,
  scale = FALSE,
  desc = FALSE,
  transform = "none",
  R = 1000,
  univariateTest = "SW",
  univariatePlot = "none",
  multivariatePlot = "none",
  multivariateOutlierMethod = "none",
  bc = FALSE,
  bcType = "rounded",
  showOutliers = FALSE,
  showNewData = FALSE
)

mvn(
  Ctrl,
  subset = NULL,
  mvnTest = "royston",
  covariance = TRUE,
  tol = 1e-25,
  alpha = 0.5,
  scale = FALSE,
  desc = FALSE,
  transform = "none",
  R = 1000,
  univariateTest = "SW",
  univariatePlot = "none",
  multivariatePlot = "none",
  multivariateOutlierMethod = "none",
  bc = FALSE,
  bcType = "rounded",
  showOutliers = FALSE,
  showNewData = FALSE
)

#Data follows multivariate normal distribition

## Covariance matrices equality
#Box's M test
library(biotools)
boxM(df_M[,2:7],as.factor(df_M[,1]))

#Covariance matrices equality is achieved

# data preparation
Group=(as.matrix(df_M[,1]))
dep_vars=as.matrix(df_M[,c(2:7)])

#MANOVA
Manova<-manova(dep_vars~Group)
summary(Manova)




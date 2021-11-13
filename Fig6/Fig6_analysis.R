rm(list=ls())

df <- read.csv("https://raw.githubusercontent.com//casey-e//Casey-et-al-2022//main//Fig6//Fig6_dataframe.csv")
df$Cohort<-factor(df$Cohort,levels=c("1","2"))
df$Group<-factor(df$Group,levels=c("FUG-eGFP","GAD-Cre"))
View(df)



#Evaluate avoidance behavior

#Perform individual univariate tests and verify integrity of the data

#OF: Time in center
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



#DLBT: time in light
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

#DLBT: latency
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

#DLBT: number_of_entries
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


#Quasi-Poisson model since the prior model has sub-dispersion
DLBT_entries_QuasiPoisson<-glm(DLBT_number_of_entries ~ Group, df, family=quasipoisson)


#Evaluate significance with Wald's test
summary(DLBT_entries_QuasiPoisson)



#EPM: percentage of time on open arms
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

#Evaluate significance
summary(EPM_TimeOpen)


#EPM: proportion_of_entries_to_open_arms
library(lme4)
EPM_entries_OneWay<-glm(EPM_proportion_entries_to_open ~ Group, df, family=binomial, weights = EPM_total_entries)

#Check assumptions
#Evaluate residuals 
e<-residuals(EPM_entries_OneWay, type='pearson')
f<-fitted(EPM_entries_OneWay)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

disp<-sum(e^2)/df.residual(EPM_entries_OneWay)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


EPM_entries_OneWay_QuasiBin<-glm(EPM_proportion_entries_to_open ~ Group, df, family=quasibinomial, weights = EPM_total_entries)


#Evaluate significance with Wald's test
summary(EPM_entries_OneWay_QuasiBin)


####      MANOVA  #####
df_M<-subset(df, select = -c(1, 3, 9))
View(df_M)
attach(df_M)

FUG_eGFP<- subset(df_M[,2:7], df_M$Group == "FUG-eGFP")
GAD_Cre<- subset(df_M[,2:7], df_M$Group == "GAD-Cre")


## Calculate centroids
(centroid_FUG_eGFP<-colMeans(FUG_eGFP))
(centroid_GAD_Cre<-colMeans(GAD_Cre))
(total_centroid<-colMeans(df_M[,-c(1)]))
(centroids<-rbind(centroid_FUG_eGFP,centroid_GAD_Cre, total_centroid))
centroids<-data.frame(centroids)


#Mahalanobis distance to evaluate multivariate outliers
#For FUG-eGFP group
s_FUG_eGFP<-cov(FUG_eGFP)  #co-variances matrix
D2_FUG_eGFP<-mahalanobis(FUG_eGFP,centroid_FUG_eGFP,s_FUG_eGFP)
(MV_out1<-which(D2_FUG_eGFP>qchisq(.95, df=5)))

#For GAD-Gre group
s_GAD_Cre<-cov(GAD_Cre) 
D2_GAD_Cre<-mahalanobis(GAD_Cre,centroid_GAD_Cre,s_GAD_Cre)
(MV_out2<-which(D2_GAD_Cre>qchisq(.95, df=5)))


#None multivariate outliers


##Asumptions
#Multivariate normality
library(MVN)
mvn(
  GAD_Cre,
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
  FUG_eGFP,
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

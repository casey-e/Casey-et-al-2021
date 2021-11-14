rm(list=ls())

df <- read.csv("https://raw.githubusercontent.com//casey-e//Casey-et-al-2022//main//Fig7//Fig7_dataframe.csv")
df$Cohort<-factor(df$Cohort,levels=c("1","2"))
df$Group<-factor(df$Group,levels=c("Ctrl","CeADrd2KO"))
View(df)



#### Perform individual univariate tests and verify asumptions ####


## Baseline: binomial model
library(lme4)
baseline<-glm(Training_baseline_proportion_of_freezing ~ Group, df, family=binomial, weights = Training_baseline_number_of_observations)

#Check assumptions
#Evaluate residuals 
e<-residuals(baseline, type='pearson')
f<-fitted(baseline)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

disp<-sum(e^2)/df.residual(baseline)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


baseline_QuasiBin<-glm(Training_baseline_proportion_of_freezing ~ Group, df, family=quasibinomial, weights = Training_baseline_number_of_observations)


#Evaluate significance with Wald's test
summary(baseline_QuasiBin)


## Training: binomial model
#Make dataframe
library(reshape2)
training_df<-df[-c(4,5,7,9,11,12,13)]
names(training_df)[4:6]<-c("Shock1", "Shock2","Shock3")
training_df<-melt(training_df, variable.name="Shock_number", value.name = "Proportion_of_freezing")
total_obs_df<-melt(df[-c(4,5,6,8,10,12,13)], value.name = "Total_observations")
training_df["Total_observations"]<-total_obs_df["Total_observations"]
View(training_df)

library(lme4)
training<-glm(Proportion_of_freezing ~ Group*Shock_number, training_df, family=binomial, weights = Total_observations)

#Check assumptions
#Evaluate residuals 
e<-residuals(training, type='pearson')
f<-fitted(training)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

#Calculate dispersion parameter
disp<-sum(e^2)/df.residual(training)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


training_QuasiBin<-glm(Proportion_of_freezing ~ Group*Shock_number, training_df, family=quasibinomial, weights = Total_observations)


#Evaluate significance with Wald's test
summary(training_QuasiBin)

#Evaluate significance with LRT test
anova(training, test="LRT")




## Testing: binomial model
library(lme4)
testing<-glm(Testing_proportion_of_freezing ~ Group, df, family=binomial, weights = Testing_number_of_observations)

#Check assumptions
#Evaluate residuals 
e<-residuals(baseline, type='pearson')
f<-fitted(baseline)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)

disp<-sum(e^2)/df.residual(baseline)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


baseline_QuasiBin<-glm(Testing_proportion_of_freezing ~ Group, df, family=quasibinomial, weights = Testing_number_of_observations)


#Evaluate significance with Wald's test
summary(baseline_QuasiBin)





####      MANOVA  #####
df_M<-subset(df, select = -c(1, 3, 5,7,9,11,13))
View(df_M)
attach(df_M)

Ctrl<- subset(df_M[,2:6], df_M$Group == "Ctrl")
CeADrd2KO<- subset(df_M[,2:6], df_M$Group == "CeADrd2KO")


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
s_CeADrd2KO<-cov(CeADrd2KO)  #co-variances matrix
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

#Data does not follows multivariate normal distribition, probably because it is a binomial variable


## Covariance matrices equality
#Box's M test
library(biotools)
boxM(df_M[,2:6],as.factor(df_M[,1]))

#Covariance matrices equality is achieved

# data preparation
Group=(as.matrix(df_M[,1]))
dep_vars=as.matrix(df_M[,c(2:6)])

#MANOVA
Manova<-manova(dep_vars~Group)
summary(Manova)




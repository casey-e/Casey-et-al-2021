rm(list=ls())

df <- read.csv("https://raw.githubusercontent.com//casey-e//Casey-et-al-2022//main//Fig5//Fig5_dataframe.csv")
df$Cohort<-factor(df$Cohort,levels=c("1","2"))
df$Group<-factor(df$Group,levels=c("Ctrl","CeADrd2KO"))
df$Day<-factor(df$Day,levels=c("1","2","3"))
View(df)


## Repeated measures model ##
library(nlme)

# Set composed symmetry model and first order autoregressive symmetry model and select the best model using Akaike information criterion

#Model with composed symmetry
distance_rm_CompSym<-gls(OF_distance_cm ~ Group*Day, correlation = corCompSymm(form = ~ 1 |MouseId), data=df)

#Model with first order autoregresive symmetry
distance_rm_AR1<-gls(OF_distance_cm ~ Group*Day, correlation = corAR1(form = ~ 1 |MouseId), data=df)

#Akaike information criterion
aic=AIC(distance_rm_CompSym,distance_rm_AR1)
aic

if (aic[1,2]-aic[2,2]>2){
  print(paste('Akaike of AR1 model is lower than the Akaike of CompSymm model by ',as.character(round((aic[1,2]-aic[2,2]), 2)), '. Ar1 model is selected'))
}else{
  print('AR1 model does not fit the data better than CompSymm model')
}


#Check assumptions
#Evaluate residuals 
e<-residuals(distance_rm_AR1, type='pearson')
f<-fitted(distance_rm_AR1)

windows()
plot(f, e, xlab="Predicted", ylab="Residuals",main="Residuals vs. predicted values" )
abline(0,0)
windows()
qqnorm(distance_rm_AR1, abline = c(0,1)) #QQ plot
#Shapiro test to evaluate normality
shapiro.test(e)
#Levene test to evaluate homocedasticity
library(car)
leveneTest(OF_distance_cm ~ Group*Day, data=df)

#Evaluate significance
anova(distance_rm_AR1)




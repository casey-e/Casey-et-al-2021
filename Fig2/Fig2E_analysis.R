
rm(list=ls())

#Load dataframe
df <- read.csv("https://raw.githubusercontent.com//casey-e/Casey-et-al-2022/main/Fig2/Fig2E_dataframe.csv")
attach(df)
View(df)

#Reorganize levels of factors
df$Treatment<-factor(df$Treatment,levels=c("Vehicle","Quinpirole"))
df$Genotype<-factor(df$Genotype,levels=c("Wild type","Drd2KO"))

# GLMM model with poisson distribution
library(lme4)
TwoWayPoisson<- glmer(Number_of_cFOS ~ Treatment*Genotype + (1|MouseId), df, family=poisson)
summary(TwoWayPoisson)
# Poisson model failed to converge.



#Negative Binomial model
library(glmmTMB)
TwoWayNegBin <- glmmTMB(Number_of_cFOS ~  Treatment*Genotype + (1|MouseId), data=df, family="nbinom2") 

#Check assumptions
#Graphical eveluation of residuals
windows()
e2 <- resid(TwoWayNegBin, type = "pearson")
F2 <- fitted(TwoWayNegBin, type ="response")
plot(x = F2, 
     y = e2, 
     xlab = "Predicted", 
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)


#Dispersion paremeter
disp<-sum(e2^2)/df.residual(TwoWayNegBin)

if (disp>1.2){
  print(paste("dispersion parameter = ",as.character( disp),", the model has overdispersion"))
}else if (disp<0.7){
  print(paste("dispersion parameter = ",as.character( disp ),", the model has subdispersion"))
}else{
  print(paste("dispersion parameter = ",as.character( disp),", the model doesn't have overdispersion nor subdispersion"))
}


#Evaluate significance with likelihood ratio test (LRT):

#Make netled models:
null_model<-glmmTMB(Number_of_cFOS ~  1 + (1|MouseId), data=df, family="nbinom2") 
genotype_model<-glmmTMB(Number_of_cFOS ~  Genotype + (1|MouseId), data=df, family="nbinom2") 
additive_model<-glmmTMB(Number_of_cFOS ~  Genotype+Treatment + (1|MouseId), data=df, family="nbinom2") 

#LRT test
lrt=anova(null_model,genotype_model,additive_model, TwoWayNegBin)
print('Likelihood ratio test:')
print(lrt)

#Post hoc comparissions

#Bolker's patch
recover.data.glmmTMB <- function(object, ...) {
  fcall <- getCall(object)
  recover.data(fcall,delete.response(terms(object)),
               attr(model.frame(object),"na.action"), ...)
}
lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                               mode = "asymptotic", component="cond", ...) {
  if (mode != "asymptotic") stop("only asymptotic mode is available")
  if (component != "cond") stop("only tested for conditional component")
  if (missing(vcov.)) 
    V <- as.matrix(vcov(object)[[component]])
  else V <- as.matrix(.my.vcov(object, vcov.))
  dfargs = misc = list()
  if (mode == "asymptotic") {
    dffun = function(k, dfargs) NA
  }
  ## use this? misc = .std.link.labels(family(object), misc)
  contrasts = attr(model.matrix(object), "contrasts")
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = contrasts)
  bhat = fixef(object)[[component]]
  if (length(bhat) < ncol(X)) {
    kept = match(names(bhat), dimnames(X)[[2]])
    bhat = NA * X[1, ]
    bhat[kept] = fixef(object)[[component]]
    modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
    nbasis = estimability::nonest.basis(modmat)
  }
  else nbasis = estimability::all.estble
  list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
       dfargs = dfargs, misc = misc)
}


#Tukey test for multiple comparissions
library(lsmeans)
lsm<-lsmeans(TwoWayNegBin,pairwise ~ Treatment*Genotype) ## Tukey 

print('Tuckey test for multiple comparissions:')
print(lsm$contrasts)





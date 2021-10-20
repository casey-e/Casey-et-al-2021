
rm(list=ls())

# Load dataframe
df <- read.csv("https://raw.githubusercontent.com//casey-e/Casey-et-al-2022/main/Fig2/Fig2B_cocaine_dataframe.csv")
attach(df)
View(df)

#Reorganize levels of factors
df$Treatment<-factor(df$Treatment,levels=c("Vehicle","Cocaine"))
df$Type_of_neuron<-factor(df$Type_of_neuron,levels=c("PKCd-","PKCd+"))

# GLMM model with poisson distribution
library(lme4)
TwoWayPoisson<- glmer(Number_of_cFOS ~ Treatment*Type_of_neuron + (1|MouseId), data=df, family=poisson)

#Check assumptions
#Graphical eveluation of residuals
windows()
plot(TwoWayPoisson) #Residuals of Pearson vs Predicted


#Overdispersion function
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

disp=overdisp_fun(TwoWayPoisson)

if (disp[4]<0.05){
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model has significant overdispersion (p=", as.character(disp[4]), ")"))
}else if (disp[2]<0.7){
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model has subdispersion"))
}else{
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model doesn't have significant overdispersion (p=", as.character(disp[4]), ")"))
}


#Negative Binomial model

TwoWayNegBin<- glmer.nb(Number_of_cFOS ~ Treatment*Type_of_neuron + (1|MouseId), data=df)

#Check assumptions
#Graphical eveluation of residuals
windows()
plot(TwoWayNegBin) #Residuals of Pearson vs Predicted

disp=overdisp_fun(TwoWayNegBin)

if (disp[4]<0.05){
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model has significant overdispersion (p=", as.character(disp[4]), ")"))
}else if (disp[2]<0.7){
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model has subdispersion"))
}else{
  print(paste("dispesion parameter = ",as.character( disp[2] ),", the model doesn't have significant overdispersion (p=", as.character(disp[4]), ")"))
}


#Evaluate significance with likelihood ratio test (LRT):

#Make netled models:
null_model<- glmer.nb(Number_of_cFOS ~ (1|MouseId), data=df)
neuron_type_model <- glmer.nb(Number_of_cFOS ~ Type_of_neuron + (1|MouseId), data=df)
additive_model<- glmer.nb(Number_of_cFOS ~ Treatment + Type_of_neuron + (1|MouseId), data=df)

#LRT test
lrt=anova(null_model,neuron_type_model,additive_model, TwoWayNegBin)
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
lsm<-lsmeans(TwoWayNegBin,pairwise ~ Treatment*Type_of_neuron) ## Tukey 
print(lsm$contrasts)

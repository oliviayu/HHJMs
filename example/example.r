###############################################
## An example of how to apply the R code. #####
###############################################

library(lme4)
library(tictoc)
library(survival)
## source all the R code
srcpath = "~/repository//HHJMs/R"
setwd(srcpath)
(file.sources = list.files(pattern="*.r$"))
sapply(file.sources,source,.GlobalEnv)

## read data
setwd("~/repository/HHJMs/example/")
long.data <- read.csv("Longdata.csv")
surv.data <- read.csv("Survdata.csv")

################################################
#### Fit joint model using two-step method  ####
##########  to get starting values  ############
################################################

# (1) Model 1: a LME model of Y 
fm1 <- y ~ 1+ year+year2+sindoes+(1|sid)
md1 <- lmer(fm1, data=long.data)

# get the estimated random effect by observation
estBi <- data.frame(row.names(ranef(md1)$sid), 
                    scale(ranef(md1)$sid, center=T,scale=T))
names(estBi) <- c("sid", "estb11")
mydat <- merge(long.data, estBi, by='sid',all=T)


# (2) Model 2: a GLME model of C
# using the estimated random intercept from model 1 (i.e. estb11, scaled) as a covariate
fm2 <- c ~ 1+year+year2+sindoes+ estb11
md2 <- glm(fm2, family=binomial, data=mydat)

# (3) Model 3: a GLME model of Z
# using the estimated random intercept from model 1 (i.e. estb11, scaled) as a covariate
fm3 <- z ~ 1+month+sindoes+doesW+estb11+(month-1|sid)
md3 <- glmer(fm3, family="binomial", data=mydat)

# (4) Model 4: a survival model
# using the estimated random intercept from model 1 (i.e. estb11, scaled) 
# and the estimate random slope from model 3 (i.e. estb21, scaled) as covariates.
Sdata <- surv.data
Sdata$estb11 <- scale(ranef(md1)$sid[,1], center=T, scale=T)
Sdata$estb21 <- scale(ranef(md3)$sid[,1], center=T, scale=T)

# a Cox PH model
fitCOX1 <- coxph(Surv(obs_time, event) ~ base+ estb11+estb21, data = Sdata)   

# a Weibull model
fitCOX2 <- survreg(Surv(obs_time, event) ~ base+ estb11+estb21, data = Sdata, dist='weibull')


###############################
######   Joint modeling   #####
###############################
### (1) Create objects for longitudinal models
LLOQ = 2 # lower limit of quantification of Y

CenObject <- list(
  fm= c ~ 1  + year +year2 +sindoes + b11,
  family='binomial', par='eta', ran.par="b11",
  disp="eta4",
  lower= -Inf, upper=Inf,
  str_val=coef(md2),
  Cregime=1,
  truncated=T, delim_val=LLOQ)

glmeObject1 <- list(
  fm= y ~ 1 + year + year2 + sindoes + b11,
  family='normal', par="beta", ran.par='b11', sigma='sigma',
  disp='beta4',
  lower=0, upper=Inf,   
  str_val=c(fixef(md1), sd(ranef(md1)$sid[,1])),
  CenObject=CenObject)

glmeObject2 <- list(
  fm= z ~ 1 + month + sindoes +doesW + b11 + month*b21,
  family='binomial', par="alpha", ran.par=c('b11','b21'),
  sigma=NULL, 
  disp=c("alpha4", 'alpha5'),
  lower=c(-Inf,0),
  upper=rep(Inf,2),
  str_val=c(fixef(md3), sd(ranef(md3)$sid[,1])),
  CenObject=NULL)

glmeObject <- list(glmeObject1, glmeObject2)


### (2) Create object for survival model
## case 1: if a Cox PH model is fit for survival data
survObject1 <- list(
  fm= obs_time ~ base+b11+b21,
  event="event", 
  par='lambda',
  disp=NULL,
  lower=NULL, upper=NULL,
  distribution=NULL,
  str_val= summary(fitCOX1)$coeff[,1])

## case 2: if a Weibull model is fit for survival data
survObject2 <- list(
  fm= obs_time ~ base+b11+b21,
  event="event", 
  par='lambda',
  disp=NULL,
  lower=NULL, upper=NULL,
#   disp = c("Wlogscale", "Wshape"),
#   lower = c(-Inf, 0), upper = rep(Inf,2),
  distribution='weibull',
  str_val= -summary(fitCOX2)$coeff[-1]/summary(fitCOX2)$scale)


######################################################
#####  Joint modeling using h-likelihood method  #####
######################################################
### case 1: with a Cox PH model ###
tic()
testjm1 <- try(JMfit(glmeObject, survObject1, 
                    long.data, surv.data,
                    idVar="sid", eventTime="obs_time",
                    survFit=fitCOX1,
                    method = "h-likelihood"), silent=T)
# re-estimate SDs of parameter estimates by using the adaptive GH method
new_sd1 = JMsd_aGH(testjm1, ghsize=4, srcpath=srcpath, paralle=T)
ptm <- toc()
(ptm$toc-ptm$tic)/60   # takes 3.33 min
# return coefficient table of Fixed effects 
JMsummary(testjm1)
JMsummary(testjm1, newSD=new_sd1)

### case 2: with a Weibull model ###
tic()
testjm2 <- try(JMfit(glmeObject, survObject2, 
                     long.data, surv.data,
                     idVar="sid", eventTime="obs_time",
                     survFit=fitCOX2,
                     method = "h-likelihood"), silent=T)
# re-estimate SDs of parameter estimates
new_sd2 = JMsd_aGH(testjm2, ghsize=4, srcpath=srcpath, paralle=T)
ptm2 <- toc()
(ptm2$toc-ptm2$tic)/60  # 10.2 min

JMsummary(testjm2)
JMsummary(testjm2, newSD=new_sd2)

######################################################
#####  Joint modeling using adaptive GH method  ######
######################################################
### the survival data must be modelled by a Weibull model 
tic()
testjm3 <- try(JMfit(glmeObject, survObject2, 
                     long.data, surv.data,
                     idVar="sid", eventTime="obs_time",
                     survFit=fitCOX2,
                     method = "aGH", ghsize=3, 
                     srcpath=srcpath, parallel=T), silent=T)
ptm3 <- toc()
(ptm3$toc-ptm3$tic)/60  # 23.6 min
JMsummary(testjm3)

### Example 

In this example, we fit joint models using simulated data. The longitudinal data contain three longitudinal responses, $z$, $y$, and $c$, where $z$ is binary, $y$ is continuous and left-truncated due to a lower limit of quantification (LLOQ), and $c$ is a truncation indicator of $y$ such that $c=1$ if $y$ is truncated and $c=0$ otherwise. The survival data contain the observed event time $obs_time$ and the event indicator $event$. Moreover, $sindoes$, $doesW$, $year$, and $year2$ are some time variables used as explanatory variables in the models, and $sid$ indicates the subject ID in both longitudinal and survival data. See data_description.md for more details. 

Download the source code in 'src' to local computer. Then call the R functions using the following R command. 

```r
## source all the R code
srcpath = "~/Dropbox/HHJMs/src"
setwd(srcpath)
file.sources = list.files(pattern="*.r$")
sapply(file.sources, source, .GlobalEnv)

# You may want to reset the working directory if your data are stored in a different location. 
setwd("~/Dropbox/HHJMs/example/")
## read data
long.data <- read.csv("Longdata.csv")
surv.data <- read.csv("Survdata.csv")

library(lme4)
library(tictoc) # for timing R scripts
```


We first fit the models separately using the two-step method. The resulting estimates will be used as starting values for the joint modeling later. 
```r
############################################################
#### Fit joint model using two-step method  to get starting values  ####
############################################################
# (1) Model 1: a LME model of Y 
fm1 <- y ~ 1+ year+year2+sindoes+(1|sid)
md1 <- lmer(fm1, data=long.data)

# get the estimated random effect
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
# using the estimated random intercept from model 1 (i.e. estb11, scaled) and the estimate random slope from model 3 (i.e. estb21, scaled) as covariates.
Sdata <- surv.data
Sdata$estb11 <- scale(ranef(md1)$sid[,1], center=T, scale=T)
Sdata$estb21 <- scale(ranef(md3)$sid[,1], center=T, scale=T)

# case 1: a Cox PH model
fitCOX1 <- coxph(Surv(obs_time, event) ~ base+ estb11+estb21, data = Sdata)   
# case 2: a Weibull model
fitCOX2 <- survreg(Surv(obs_time, event) ~ base+ estb11+estb21, data = Sdata, dist='weibull')
```

Then we create longitudinal and survival objects.  
```r
### (1) Create objects for longitudinal models
LLOQ = 2  # LLOQ of Y

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
  distribution='weibull',
  str_val= -summary(fitCOX2)$coeff[-1]/summary(fitCOX2)$scale)
```

If a Cox PH model is used for survival data, we can fit a joint model via h-likelihood using the code below.
```r
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
```

If a Weibull model is used for survival data, we can fit a joint model via h-likelihood using the code below.
```r		
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
```

As comparison, we also implement the adaptive GH method for the case when Weibull survival model is used. 
```r
######################################################
#####  Joint modeling using adaptive GH method  ######
######################################################
### the survival data must be modelled by a Weibull model 
tic()
testjm3 <- try(JMfit(glmeObject, survObject2, 
                     long.data, surv.data,
                     idVar="sid", eventTime="obs_time",
                     survFit=fitCOX2,
                     method = "aGH", ghsize=3, srcpath, parallel=T), silent=T)
ptm3 <- toc()
(ptm3$toc-ptm3$tic)/60  
```

The function *JMsummary()* returns the coefficient table for the joint modeling via h-likelihood method or adaptive GH method. For example,
```r
JMsummary(testjm1, newSD=NULL, digits=3)  
# The Std.Error were obtained based on the h-likelihood method.
         Estimate  Std.Error  Zvalue  Pvalue
# beta0      2.079     0.051  40.462  0.000
# beta1      0.941     0.040  23.293  0.000
# beta2     -0.287     0.017 -16.867  0.000
# beta3      1.488     0.022  67.596  0.000
# eta0       0.201     0.232   0.868  0.385
# eta1      -4.285     0.378 -11.326  0.000
# eta2       1.354     0.180   7.523  0.000
# eta3      -6.410     0.374 -17.127  0.000
# alpha0    -1.700     0.099 -17.125  0.000
# alpha1     0.165     0.013  12.892  0.000
# alpha2     1.938     0.110  17.621  0.000
# alpha3    -0.046     0.005  -8.805  0.000
# lambda0   -0.999     0.120  -8.302  0.000
# lambda1   -1.264     0.134  -9.407  0.000
# lambda2    2.147     0.141  15.212  0.000


JMsummary(testjm1, newSD=new_sd1, digits=3)
# The Std.Error=new_sd1 were obtained based on the adaptive GH method.
#       Estimate Std.Error  Zvalue Pvalue
# beta0      2.079     0.062  33.644  0.000
# beta1      0.941     0.051  18.333  0.000
# beta2     -0.287     0.022 -13.065  0.000
# beta3      1.488     0.028  53.832  0.000
# eta0       0.201     0.282   0.714  0.475
# eta1      -4.285     0.477  -8.974  0.000
# eta2       1.354     0.227   5.954  0.000
# eta3      -6.410     0.486 -13.185  0.000
# alpha0    -1.700     0.127 -13.359  0.000
# alpha1     0.165     0.016  10.441  0.000
# alpha2     1.938     0.140  13.834  0.000
# alpha3    -0.046     0.007  -6.944  0.000
# lambda0   -0.999     0.176  -5.688  0.000
# lambda1   -1.264     0.200  -6.308  0.000
# lambda2    2.147     0.222   9.670  0.000


```



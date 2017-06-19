# HHJMs
H-likelihood based Hierarchical Joint Models

### Description

This package fits shared parameter models for the joint modeling of longitudinal data and survival data, where the longitudinal responses may be of mixed types, such as binary and continuous, and may be left censored by lower limit of quantification. For statistical inference, we consider a computationally efficient approximate likelihood method based on h-likelihood method [3]. There is an extensive literature on h-likelihood method (e.g. [1]-[5]). Essentially, the h-likelihood method uses Laplace approximations to the intractable integral in the likelihood. Moreover, it can produce approximate MLEs for the mean parameters and approximate restricted maximum likelihood estimates (REML) for the variance-covariance (dispersion) parameters. 

We also implement the adaptive Gauss-Hermite method to compare with the h-likelihood method. 

A detailed example is given in the **example** folder.

### Usage
```r
# fit joint model
JMfit ( glmeObject = list( ), survObject = list( ), long.data, surv.data, idVar, eventTime, survFit, method, 
          itertol=0.001, Ptol=0.01, epsilon=1e-6, iterMax=10, ghsize=4, Silent = T )

# estimate SEs of h-likelihood based parameter estimates, by using the adaptive GH method
# needed on when method='h-likelihood' in JMfit()
JMsd_aGH( JMobject, ghsize=4, srcpath=NULL, paralle=F )

# print coefficient table 
JMsummary( JMobject, newSD=NULL, digits=3 )

```


##### Arguments

|-----------|-----------|
|glmeObject | A list, indicating the GLME models to be fitted. [Details](../master/others.md) | 
|survObject | A list, indicating the survival model (either Cox PH or Weibull model) to be fitted.  [Details](../master/others.md) |
|long.data  | longitudinal data containing the variables named in formulas in glmeObject |
|surv.data  | survival data containing the variables named in formulas in survObject |
|idVar      | subject id |
|eventTime |    observed event time    |
| survFit | an object returned by the *coxph()* or *survreg()* function to represent a fitted survival model from the two-step method. |
| method | a vector indicating which method to apply. If *method*='h-likelihood' (by default), the h-likelihood method is used; if *method*='aGH', the adaptive Gauss-Hermite method is used.   |
|itertol    | Convergence tolerance on the relative absolute change in log-likelihood function between successive iterations. Convergence is declared when the change is less than itertol. Default is itertol = 0.001. |
|Ptol    | Convergence tolerance on the average relative absolute change in parameter estimates between successive iterations. Convergence is declared when the change is less than Ptol. Default is Ptol = 0.01. |
| epsilon |   A small numerical value, used to calculate the numerical value of the derivative of score function. The default value is 1e-6.|
|iterMax    | The maximum number of iterations. The default value is 10. |
|ghsize | The number of quadrature points used in the adaptive GH method. The default value is 4. |
|Silent     | logical: indicating if messages about convergence success or failure should be suppressed. |


|-----------|-----------|
| JMobject | output of *JMfit()*|
| srcpath |  a character vector of full path names indicating the location of the R code; needed if *parallel*=T;  *srcpath*=NULL by default. |
| parallel | logical: indicating if compute the standard errors of parameter estimates in parallel. By default, *paralle=F*. | 

|-----------|-----------|
| newSD  | If *newSD*=NULL (by default), the p-values of the parameter estimates are calculated based on the SEs in the output of *JMfit()* i.e. *fixedsd*. Otherwise, the p-values are calculated based on the new SEs given by *newSD*, e.g. the output of *JMsd_aGH()*.  | 


##### Outputs
|           |          |
|-----------|-----------|
|fixedest  |  a named vector of estimated coefficients|
|fixedsd  | standard errors of estimated coefficients, calculated based on the given *method*. |
|Bi |  estimated random effects, corresponding to each subject|
|B | estimated random effects, corresponding to each measurement|
|covBi |  estimated covariance matrix of the random effects|
|sigma |  estimates of dispersion parameters |
|convergence| An integer code indicating type of convergence: 0 indicates successful convergence, 1 indicates that the maximum limit for iterations 'iterMax' has been reached without convergence. |
|loglike_value|  value of approximate log likelihood function  |
| ... | ... |



<!--
```r
HHJMsummary( object, digits)
```
##### Arguments
|           |          |
|-----------|-----------|
| object |  an object for which a summary is desired |
| digits |  integer indicating  the number of decimal places to be used|

##### Output
-->

<!--
### Example 
To use this package, first download the source code in 'src' to your local computer. Then call the R functions using the following R command. 

```r
# The 'src' folder has been downloaded to my desktop.
setwd("~/desktop/src")
file.sources = list.files(pattern="*.r$")
sapply(file.sources, source, .GlobalEnv)

# You may reset the working directory if your data are stored in a different location. 
```
In this example, we fit joint models using simulated data. The longitudinal data contain three longitudinal responses, $z$, $y$, and $c$, where $z$ is binary, $y$ is continuous and left-censored due to lower limit of quantification, and $c$ is a censoring indicator of $y$ such that $c=1$ if $y$ is cneosred and $c=0$ otherwise. The survival data contain the observed event time $obs_time$ and the event indicator $event$. Moreover, $sindoes$, $does30$, $t365$, and $t2$ are some time variables used as explanatory variables in the models, and $patientID$ indicates the subject ID in both longitudinal and survival data.

```r
  long.data <- read.csv("LongData.csv", head=T)
  surv.data <- read.csv("SurvData.csv", head=T)

   glmeObject1 <- list(
    fm = z ~ 1 + sindoes + does30 + t365 + (1 | patientID),
    family='binomial',
    par='alpha',
    ran.par='b1',
    sigma=NULL,
    str_val=fixef(fit1),   # the initial values are generated by fitting the model using glmer(), see fit1 below
    CenObject=NULL
  )
  
  CenObject <- list(
    fm=c ~ 1 +sindoes+t365+(1|patientID),
    family='binomial',
    par='eta',
    ran.par='a',
    str_val=fixef(fit3)    # see fit3 below
  )
  
  glmeObject2 <- list(
    fm=y ~ 1 + sindoes + t365 + t2 + (1 | patientID),
    family='normal',
    par='beta',
    ran.par='b2',
    sigma='sigma',
    str_val=fixef(fit2),   # see fit2 below
    CenObject=CenObject
  )
  
  survObject <- list(
    fm = obs_time ~ age,
    event='event',
    par='lambda',
    str_val=summary(fit4)$coeff[,1]    # see fit4 below
  )
  
  set.seed(10) 

  testjm <- HHJMfit(glmeObject=list(glmeObject1, glmeObject2), 
                    survObject,
                    long.data, surv.data, 
                    idVar="patientID",
                    itertol=0.001)
                    
  HHJMsummary(testjm, digits=4) 
  
                  
＃＃ The following models are fitted to generate starting values of the fixed parameters. 
＃  library(lme4)
＃  fit1 <- glmer(z ~ 1 + sindoes + does30 + t365 + (1 | patientID), data=long.data, family='binomial')
＃  fit2 <- lmer(y ~ 1 + sindoes + t365 + t2 + (1 | patientID), data=long.data)
＃  fit3 <- glmer(c ~ 1 + sindoes + t365 + (1 | patientID), data=long.data, family="binomial")
＃ surv.data$nb1 <- ranef(fit1)$patientID
＃ surv.data$nb2 <- ranef(fit2)$patientID 
＃ surv.data$nb3 <- ranef(fit3)$patientID
＃ fit4 <- coxph(Surv(obs_time, event) ~ age + nb1 + nb2 + nb3, data = surv.data)
  
```

The fitting results based on separate analyses and HHJMs are summarized below respectively.

```r
### fitting results based on separate analyses

        Estimate  Std. Error  z value    Pvalue
alpha0   -1.778      0.069   -25.832     0
alpha1    2.017      0.084    24.009     0
alpha2    -0.228     0.016   -13.987     0
alpha3    1.367      0.031    43.710     0
beta0     1.072      0.039    27.323     0
beta1     1.637      0.029    56.262     0
beta2     1.885      0.030    62.964     0
beta3    -0.446      0.012   -36.429     0
eta0      0.585      0.068     8.588     0
eta1     -0.555      0.085    -6.550     0
eta2     -2.018      0.045   -44.710     0
lambda1   0.247      0.027     9.162     0
Asso1    -2.249      0.299    -7.528     0
Asso2    -2.154      0.391    -5.516     0
Asso3    -3.757      0.459    -8.190     0

### fitting results based on HHJMs

        Estimate  Std.Error z.value Pvalue
alpha0    -1.754     0.064 -27.263      0
alpha1     1.985     0.083  23.985      0
alpha2    -0.224     0.016 -13.902      0
alpha3     1.341     0.030  45.052      0
beta0      1.001     0.052  19.265      0
beta1      2.496     0.024 106.203      0
beta2      1.870     0.022  83.917      0
beta3     -0.514     0.009 -60.421      0
eta0       0.558     0.062   8.938      0
eta1      -0.548     0.084  -6.537      0
eta2      -1.988     0.043 -46.194      0
lambda1    0.232     0.007  31.124      0
Asso1     -2.706     0.231 -11.687      0
Asso2     -2.272     0.259  -8.762      0
Asso3     -4.402     0.361 -12.195      0

```
-->



### References
[1] Do Ha, I., Lee, Y., & Song, J. K. (2002). Hierarchical-likelihood approach for mixed linear models with censored data. Lifetime data analysis, 8(2), 163-176.

[2] Ha, I. D., Park, T., & Lee, Y. (2003). Joint modelling of repeated measures and survival time data. Biometrical journal, 45(6), 647-658.

[3] Lee, Y., Nelder, J. A., & Noh, M. (2007). H-likelihood: problems and solutions. Statistics and Computing, 17(1), 49-55.

[4] Lee, Y., & Nelder, J. A. (1996). Hierarchical generalized linear models. Journal of the Royal Statistical Society. Series B (Methodological), 619-678.

[5] Noh, M., & Lee, Y. (2007). REML estimation for binary data in GLMMs. Journal of Multivariate Analysis, 98(5), 896-915.



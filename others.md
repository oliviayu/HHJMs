
### More details on the arguments.

The *glmeObject( )* and *survObject( )* must be in the following format:
```r
  glmeObject <- list(fm, family, par, ran.par, sigma, disp, 
  lower, upper, str_val,  CenObject)

  survObject <- list(fm, event,  par, disp, 
  lower, upper, distribution, str_val)

  CenObject <- list(fm, family="binomial", par, ran.par, 
  disp, lower, upper, str_val, Cregime=1, truncated=T, delim_val)  
```


#### Arguments
|       |        |
|-------|--------| 
|fm     | A two-sided linear formula object with the response on the left of a ~ operator and the terms, separated by + operators, on the right. |
|family | A GLM family.  |
| par | A character string, naming the parameters. Such as, "alpha", "beta", ...  |
| ran.par |  A character string, naming the random effects. Such as, "b1","b2",...  |   
| sigma  |  A character string, naming the standard deviation for normal distribution. Such as, "sigma". |
| disp | A character vector, specifying the dispersion parameters.  |
| lower / upper | lower/upper bounds of the dispersion parameters specified by *disp*.|
| str_val |  A numeric vector of starting values for the fixed parameters in the model. |
| CenObject | A list, indicating the logistic GLME model used to model the censoring mechanism. CenObject=NULL means that the response variable is not censored.| 

<!--For GLME models, random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors. For Cox model, the random effects in GLME models will be automatically incorporated as explanatory variables. So there is no need to include random effects on the right side of the formula.-->

|       |        |
|-------|--------| 
| event |  event indicator  |
| distribution | If *distribution*=NULL, a Cox PH model is fitted. If *distribution*=weibull, a Weibull model is fitted.|
| delim_val |  the lower limit of quantification |
| Cregime | If *Cregime*=1 (by default), we assume that the censored data are from one regime (point mass). If *Cregime*=2, we assume that the censored data are from two regimes, one from normal distribution  and one from point mass.| 
| truncated | logical: if *truncated*=T (by default), we assume the observed response variable follows a **truncated** normal distribution; otherwise, we assume it follows a normal distribution.|


See the example in example.md for details.

<!-- 
|Cvalue | value of the lower limit of quantification|
| Cregime      | number of regimes in the censored data. It defaults to 1.| 
-->


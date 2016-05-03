
Note that *glmeObject( )* and *survObject( )* must be in the following format:
```r
  glmeObject <- list(fm, family, par, ran.par, sigma, str_val,  CenObject)
  survObject <- list(fm, event,  par, str_val)
  CenObject <- list(fm, family="binomial", par, ran.par, str_val)  
```

#### Arguments
|       |        |
|-------|--------| 
|fm     | A two-sided linear formula object with the response on the left of a ~ operator and the terms, separated by + operators, on the right. For GLME models, random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors. |
|family | A GLM family.  |
| event |  event indicator  |
| par | A character string, naming the fixed parameters. Such as, "alpha", "beta", ...  |
| ran.par |  A character string, naming the random effects. Such as, "b1","b2",...  |   
| sigma  |  A character string, naming the standard deviation for normal distribution. Such as, "sigma". |
| str_val |  A numeric vector of starting values for the fixed parameters in the model. |
| CenObject | A list, indicating the logistic GLME model used to model the censoring mechanism. CenObject=NULL means that the response variable is not censored.| 

<!-- 
|Cvalue | value of the lower limit of quantification|
| Cregime      | number of regimes in the censored data. It defaults to 1.| 
-->


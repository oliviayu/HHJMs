
### Data Description

The simulated data have 100 subjects in total, among which 69 experience the event of interest during follow-up and 31 don't. 
The longitudinal and survival data contain the following variables: 

|       |        |
|-------|--------|
| sid     | subject ID |
| time | measurement time (in days)  |
| injectionNO |  The immunization scheduled corresponding to the visit date. Equal to 1, 2, 3, 4, 5, 6, 7 for the Months 1, 2-6, 7-12, 13-18, 19-24, 25-30, 31-36. |
| doest | The time (in days) between the most recent immunization and the measurement time, denoted as $t_{d_{ij}}$. |
| doesW | The time (in weeks) between the most recent immunization and the measurement time, denoted as $t_{d_{ij}}^*$ |
| perd | The time between the previous immunization and the next, denoted as $\Delta_{ij}$. |   
| base | A baseline covariate, simulated from $N(0,1)$. |
| sindoes | $sin(\frac{\pi}{\Delta_{ij}}\times t_{d_{ij}})$, calculated to capture the periodic patterns over time.  |
| month | measurement time (in months), denoted as $t_{ij}$.  |
| month2 |  $t_{ij}^2$ |
| year | measurement time (in years), denoted as $t_{ij}^*$.  |
| year2 | $(t_{ij}^*)^2$  |
|  z |  A binary longitudinal response. |
| y | A continuous longitudinal response, left truncated with a lower limit of quantification of 2.|
| c | A truncation indicator of $y$. If $c=1$, $y$ is left truncated.|
| event | An event indicator. If $event_i=1$, the subject $i$ experiences the event of interest. |
| obs_time | Observed event time (in days). | 


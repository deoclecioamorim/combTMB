
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `combTMB`:Combined Models using Template Model Builder

> [Deoclecio J. Amorim](https://github.com/deoclecioamorim) -
> <deocleciojardim@usp.br>, ESALQ-USP

The `combTMB` is an R package that implements an extended random effects
approach to model repeated, overdispersed binary and count data (called
combined models). Furthermore, fit linear and generalized linear mixed.
The adjustment of the models is done by maximum likelihood with
[TMB](https://github.com/kaskr/adcomp) (Template Model Builder). A
common application is for superdispersed and correlated longitudinal
data. See Molenberghs, G., Verbeke, G. and Demetrio, C.G. (2007)
<doi:10.1007/s10985-007-9064-y> and Molenberghs et al. (2012)
<doi:10.1016/j.jmva.2012.05.005>.

Joint work with [Afrânio Márcio Corrêa
Vieira](https://www.des.ufscar.br/departamento/docentes/afranio-marcio-correa-vieira)
and [Clarice G.B.
Demétrio](http://ce.esalq.usp.br/equipe/clarice-garcia-borges-demetrio).

## Installation

Assuming you have a [C++
compiler](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
installed, you can install the development version of combTMB like so:

``` r
# install.packages("remotes")
remotes::install_github("deoclecioamorim/combTMB", dependencies = TRUE)
```

## Basic use

The syntax of combTMB is very similar to the
[lme4](https://github.com/lme4/lme4/) and
[glmmTMB](https://github.com/glmmTMB/glmmTMB) packages. The main
function is `combTMB(...)`.

``` r
library(combTMB)
#> Package 'combTMB' version 0.0.0.9000
#> Type 'citation("combTMB")' for citing this R package in publications.
# Poisson-Normal model --------------------------------------------------------
m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
summary(m1)
#> Family: poisson 
#> Link function: log 
#> Formula: OT ~ Period + (1 | Donor) 
#> Number of obs: 1148 
#> -2 x logLik         AIC         BIC    df.resid 
#>      7839.5      7845.5      7860.7        1145 
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.69168    0.03307  81.397  < 2e-16 ***
#> PeriodP2     0.06419    0.01859   3.454 0.000553 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#>                   Estimate Std. Error
#> Donor (Intercept)    0.268      0.025
#> 
#> Number of subjects: 318
# Methods --------------------------------------------------------------------
print(m1)
#> 
#> combTMB regression models
#> Call:  combTMB(formula = OT ~ Period + (1 | Donor),
#>            data = embryos,
#>            family = poisson(),
#>            dformula = ~1)
#> 
#> Fixed Effects:
#> (Intercept)     PeriodP2  
#>     2.69168      0.06419  
#> 
#> Residual degrees of freedom: 1145
#> -2 x log-likelihood: 7839.525
#> 
#> For more details, run the summary function
summary(m1)
#> Family: poisson 
#> Link function: log 
#> Formula: OT ~ Period + (1 | Donor) 
#> Number of obs: 1148 
#> -2 x logLik         AIC         BIC    df.resid 
#>      7839.5      7845.5      7860.7        1145 
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.69168    0.03307  81.397  < 2e-16 ***
#> PeriodP2     0.06419    0.01859   3.454 0.000553 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#>                   Estimate Std. Error
#> Donor (Intercept)    0.268      0.025
#> 
#> Number of subjects: 318
# Basic sanity checks on our model--------------------------------------------
sanitycombTMB(m1)
#> ✔ Suggests successful convergence!
#> ✔ Hessian matrix is positive definite!
#> ✔ No extreme or very small eigen values detected!
#> ✔ No fixed-effect standard errors are NA
#> ✔ No fixed-effect standard errors look unreasonably large
# Combined model: Poisson-Gamma-Normal model ---------------------------------
m2 <- combTMB(OT ~ Period+(1|Donor), family=poigamma(), data=embryos)
# Methods --------------------------------------------------------------------
print(m2)
#> 
#> combTMB regression models
#> Call:  combTMB(formula = OT ~ Period + (1 | Donor),
#>            data = embryos,
#>            family = poigamma(),
#>            dformula = ~1)
#> 
#> Fixed Effects:
#> (Intercept)     PeriodP2  
#>     2.69904      0.06813  
#> 
#> Residual degrees of freedom: 1144
#> -2 x log-likelihood: 7650.511
#> 
#> For more details, run the summary function
summary(m2)
#> Family: poigamma 
#> Link function: log 
#> Formula: OT ~ Period + (1 | Donor) 
#> Dformula: ~1 
#> Number of obs: 1148 
#> -2 x logLik         AIC         BIC    df.resid 
#>      7650.5      7658.5      7678.7        1144 
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.69904    0.03475  77.662  < 2e-16 ***
#> PeriodP2     0.06813    0.02520   2.704  0.00686 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>  
#> 
#> Random effects:
#>                   Estimate Std. Error
#> Donor (Intercept)   0.2484      0.024
#> Overdisp.(theta)   21.9242      2.498
#> 
#> Number of subjects: 318
# Basic sanity checks on our model--------------------------------------------
sanitycombTMB(m2)
#> ✔ Suggests successful convergence!
#> ✔ Hessian matrix is positive definite!
#> ✔ No extreme or very small eigen values detected!
#> ✔ No fixed-effect standard errors are NA
#> ✔ No fixed-effect standard errors look unreasonably large
```

Currently, the methods implemented for `"combTMB"` objects are

``` r
methods(class = "combTMB")
#>  [1] anova        coefDisp     df.residual  dispersion   family      
#>  [6] fitted       formula      getME        logLik       model.frame 
#> [11] model.matrix nobs         partvar      predict      print       
#> [16] ranef        residuals    simulate     summary      terms       
#> [21] vcov        
#> see '?methods' for accessing help and source code
```

## License

The `combTMB` package is licensed under the [GNU General Public License,
version 3](https://www.gnu.org/licenses/gpl-3.0.html), see file
`LICENSE.md`, © 2022 Deoclecio J. Amorim.

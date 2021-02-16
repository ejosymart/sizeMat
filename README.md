# sizeMat

[![packageversion](https://img.shields.io/badge/Package%20version-1.1.2-orange.svg?style=flat-square)](commits/master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/sizeMat)](http://cran.r-project.org/package=sizeMat)
[![metacran
downloads](http://cranlogs.r-pkg.org/badges/sizeMat)](http://cran.rstudio.com/web/packages/sizeMat/index.html)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.6-6666ff.svg)](https://cran.r-project.org/)

**Estimate Size at Sexual Maturity**

This package allows to estimate Size at Morphometric and Gonadal
Maturity for organisms, usually fish and invertebrates.

The estimation of morphometric maturity used two allometric variables
and is based on the relative growth.

The estimation of gonadal maturity used one allometric variable and the
stage of sexual maturity (gonad maturation stage).

## Install

Get the released version from CRAN:

``` r
install.packages("sizeMat")
```

Or the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ejosymart/sizeMat")
```

## Examples

This is a basic example which shows you how to estimate Size at
Morphometric and Gonad Maturity:

## Size at Morphometric Maturity

``` r
data(crabdata)

classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_height"), 
                                varSex = "sex_category", selectSex = NULL, method = "ld")
#> all individuals were used in the analysis

print(classify_data)
#> Number in juvenile group = 83 
#> 
#> Number in adult group = 140 
#> 
#> -------------------------------------------------------- 
#> 1) Linear regression for juveniles 
#> 
#> Call:
#> glm(formula = y ~ x, data = juv)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.77010  -0.57399   0.09397   0.56605   1.99008  
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -3.794687   0.497056  -7.634 3.93e-11 ***
#> x            0.161327   0.004701  34.314  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.7320842)
#> 
#>     Null deviance: 921.306  on 82  degrees of freedom
#> Residual deviance:  59.299  on 81  degrees of freedom
#> AIC: 213.63
#> 
#> Number of Fisher Scoring iterations: 2
#> 
#> -------------------------------------------------------- 
#> 2) Linear regression for adults 
#> 
#> Call:
#> glm(formula = y ~ x, data = adt)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -3.3055  -1.0932  -0.0628   1.1178   3.2759  
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -11.246726   1.199496  -9.376   <2e-16 ***
#> x             0.273837   0.008648  31.663   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 2.265729)
#> 
#>     Null deviance: 2584.24  on 139  degrees of freedom
#> Residual deviance:  312.67  on 138  degrees of freedom
#> AIC: 515.79
#> 
#> Number of Fisher Scoring iterations: 2
#> 
#> -------------------------------------------------------- 
#> 3) Difference between slopes (ANCOVA) 
#>               Estimate  Std. Error   t value     Pr(>|t|)
#> (Intercept) -3.7946869 0.757105677 -5.012097 1.109526e-06
#> x            0.1613275 0.007161179 22.528064 6.035478e-59
#> mature      -7.4520389 1.285219562 -5.798261 2.320729e-08
#> x:mature     0.1125093 0.010361046 10.858878 2.956242e-22
#> [1] "slopes are different"


my_ogive = morph_mature(classify_data, method = "fq")

print(my_ogive)
#> formula: Y = 1/1+exp-(A + B*X)
#>     Original Bootstrap (Median)
#> A   -20.753  -20.8576          
#> B   0.1748   0.1757            
#> L50 118.7237 118.6591          
#> R2  -        0.7111
```

## Size at Gonad Maturity

``` r
data(matFish)

my_ogive = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), 
                        inmName = "I", matName = c("II", "III", "IV"), method = "fq", niter = 999)

print(my_ogive)
#> formula: Y = 1/1+exp-(A + B*X)
#>     Original Bootstrap (Median)
#> A   -8.6047  -8.6407           
#> B   0.356    0.3576            
#> L50 24.1694  24.1714           
#> R2  0.5595   -
```

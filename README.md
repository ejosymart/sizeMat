sizeMat
=======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/sizeMat)](http://cran.r-project.org/package=sizeMat) [![](http://cranlogs.r-pkg.org/badges/sizeMat)](http://cran.rstudio.com/web/packages/sizeMat/index.html)

**Estimate Size at Sexual Maturity**

This package allows to estimate Size at Morphometric and Gonadal Maturity for organisms, usually fish and invertebrates.

The estimation of morphometric maturity used two allometric variables and is based on the relative growth.

The estimation of gonadal maturity used one allometric variable and the stage of sexual maturity (gonad maturation stage).

Installation
------------

Get the released version from CRAN:

``` r
install.packages("sizeMat")
```

Or the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ejosymart/sizeMat")
```

Examples
--------

This is a basic example which shows you how to estimate Size at Morphometric and Gonad Maturity:

1.  Size at Morphometric Maturity

``` r
data(crabdata)

classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
                                varSex = "sex_category", selectSex = NULL, method = "ld")
#> all individuals were used in the analysis

print(classify_data)
#> Number in juveline group = 83 
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




# plot(classify_data)

my_ogive = morph_mature(classify_data, method = "fq")

print(my_ogive)
#> formula: Y = 1/1+exp-(A + B*X)
#>        Original Bootstrap (Median)
#> A   -20.7530424        -20.7955022
#> B     0.1748011          0.1751376
#> L50 118.7237428        118.6414514




# plot(my_ogive)
```

2.  Size at Gonad Maturity

``` r
data(matFish)

my_ogive = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), 
                        inmName = "I", matName = c("II", "III", "IV"), method = "fq", niter = 999)

print(my_ogive)
#> formula: Y = 1/1+exp-(A + B*X)
#>       Original Bootstrap (Median)
#> A   -8.6046612         -8.6167068
#> B    0.3560150          0.3568128
#> L50 24.1693811         24.1629202




# plot(my_ogive)
```

sexMat
======

\*\* Estimate Size at Sexual Maturity\*\*

This package allows to estimate Morphometric and Gonadal Size at Sexual maturity for organisms, usually fish and invertebrates.

The estimation of morphometric maturity used two allometric variables and is based in the relative growth.

The estimation of gonadal maturity used one allometric variable and the stage of sexual maturity (gonad maturation stage).

Installation
------------

Get the released version from CRAN:

``` r
install.packages("sexMat")
```

Or the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ejosymart/sexMat")
```

Examples
--------

This is a basic example which shows you how to estimate Morphometric and Gonadal Size at Sexual Maturity:

1.  Morphometric Size at Sexual Maturity

``` r
data(crabdata)

classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
                                varSex = "sex_category", selectSex = NULL, method = "ld")
#> all individuals were used in the analysis

my_ogive = morph_mature(classify_data, method = "fq")

print(my_ogive)
#> formula: Y = 1/1+exp-(A + B*X)
#>        Original Bootstrap (Median)
#> A   -20.7530424        -20.7955022
#> B     0.1748011          0.1751376
#> L50 118.7237428        118.6414514

plot(my_ogive)
```

![](README-unnamed-chunk-2-1.png)![](README-unnamed-chunk-2-2.png)![](README-unnamed-chunk-2-3.png)![](README-unnamed-chunk-2-4.png)

    #> Morphometric size at sexual maturity = 118.6 
    #> Confidence intervals = 115.9 - 121

1.  Gonadal Size at Sexual Maturity

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

plot(my_ogive)
```

![](README-unnamed-chunk-3-1.png)![](README-unnamed-chunk-3-2.png)![](README-unnamed-chunk-3-3.png)![](README-unnamed-chunk-3-4.png)

    #> Morphometric size at sexual maturity = 24.2 
    #> Confidence intervals = 23.8 - 24.6

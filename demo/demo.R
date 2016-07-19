library(sexMat)

#  ------------------------------------------------------------------------
# Morphometric size at sexual maturity ------------------------------------
#  ------------------------------------------------------------------------

# Load data ---------------------------------------------------------------
data(crabdata)

# Classify data -----------------------------------------------------------
classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
varsex = "sex_category", selectSex = NULL, method = "ld")

print(classify_data)

plot(classify_data, col = c(1, 2), pch = c(1, 2), cex = c(1, 1), 
     lty_lines = c(6, 1), lwd_lines =  c(1, 4))

# Estimate Morphometric Size at Sexual Maturity based on Relative Growth --
my_ogive = morph_mature(classify_data, method = "bayes", niter = 1000)

print(my_ogive)

plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"))




#  ------------------------------------------------------------------------
# Gonadal size at sexual maturity -----------------------------------------
#  ------------------------------------------------------------------------

# Load data ---------------------------------------------------------------
data(matFish)

# Estimate Gonadal Size at Sexual Maturity
my_ogive = gonad_mature(matFish, inmName = "I", matName = c("II", "III", "IV"), method = "fq", niter = 1000)

print(my_ogive)

plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"))

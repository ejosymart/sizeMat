
# Load data ---------------------------------------------------------------
data(crabdata)

# Classify data -----------------------------------------------------------
classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
varsex = "sex_category", selectSex = NULL, method = "ld")
plot(classify_data, col = c(1, 2), pch = c(1, 2), cex = c(1, 1), lty_lines = c(6, 1), lwd_lines =  c(1, 4))

# Estimate size at sexual maturity ----------------------------------------
my_ogive = calculate_ogive(classify_data, method = "fq", niter = 1000)
print(my_ogive)
plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"))

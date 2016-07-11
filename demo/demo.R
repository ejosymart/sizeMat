
# Load data ---------------------------------------------------------------
data(crabdata)

# Classify data -----------------------------------------------------------
classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
varsex = "sex_category", useSex = NULL, method = "qd")
plot(classify_data, col = c(1, 2), pch = c(4, 5))

# Estimate size at sexual maturity ----------------------------------------
my_ogive = calculate_ogive(classify_data, method = "bayes")
print(my_ogive)
plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"))

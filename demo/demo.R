
# Load data ---------------------------------------------------------------
my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")

# Read database -----------------------------------------------------------
data    = read_data(my_file)
data

# Classify data --------------------------------------------------------
classify_data = classify_mature(data, method = "ld")
plot(classify_data)


# Calculate ogive ---------------------------------------------------------
my_ogive = calculate_ogive(classify_data, method = "bayes")
plot(my_ogive, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"))





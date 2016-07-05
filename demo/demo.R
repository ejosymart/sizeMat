
# Load data ---------------------------------------------------------------
my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")

# Read database -----------------------------------------------------------
data    = read_data(my_file)
print(data)

# Classify data --------------------------------------------------------
classify_data = classify_mature(data, sex = "all", method = "ld")
plot(classify_data)


# Calculate ogive ---------------------------------------------------------
my_ogive = calculate_ogive(classify_data, method = "fq")
print(my_ogive)
plot(my_ogive, xlab1 = "X", ylab1 = "fitted (Proportion mature)", col = c("blue", "red"))


my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
data    = read_data(my_file)
classify_data = classify_mature(data)
plot(classify_data)
my_ogive = calculate_ogive(classify_data, methodReg = "fq")
plot(my_ogive, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"))

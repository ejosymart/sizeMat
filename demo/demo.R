
#Read Data
file = "crabdat.txt"
crabdat = read_data(file)

#View data
head(crabdat)


#Classify juvelines and adults (PCA + hierarchical clustering + linear or quadratic discriminant analysis)
## linear or quadratic discriminant analysis based on the homogenity covariance matrix
my.mat1 = classify_mature(data = crabdat)
my.mat1


#Plot classify
plot(my.mat1)


#Calculate ogive
my.ogive1 = calculate_ogive(my.mat1, method = "fq")
my.ogive2 = calculate_ogive(my.mat1, method = "bayes")

#Plot
plot(my.ogive1)
plot(my.ogive2)


#install and load this packages
library(UsingR)
library(matrixStats)
library(MCMCpack)
library(biotools)
library(bPCA)

rm(list=ls())
source("calculate_ssm.R")

#Read Data
file = "crabdat"
crabdat = read_data(file, ext = "txt")

#View data
headtail(crabdat)


#Classify juvelines and adults (BayesianPCA + hierarchical clustering + linear or quadratic discriminant analysis)
## linear or quadratic discriminant analysis based on the homogenity covariance matrix
my.mat1 = classify_mature(data = crabdat, method = "bayes")
my.mat1


#Plot classify
plot(my.mat1)


#Calculate ogive
my.ogive1 = calculate_ogive(my.mat1, method = "fq")
my.ogive2 = calculate_ogive(my.mat1, method = "bayes")

#Plot
plot(my.ogive1)
x11()
plot(my.ogive2)


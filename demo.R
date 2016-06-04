#install and load this packages
library(gdata)
library(rJava)
library(xlsx)
library(UsingR)
library(matrixStats)
library(MCMCpack)

rm(list=ls())
source("MatureCode.R")

#Read Data
file = "crabdat"
crabdat = dataMature(file, ext = "txt")

#View data
headtail(crabdat)

#Classify juvelines and adults (minimun distance)
my.mat1 = classifyByDistance(data = crabdat, n.iter = 100)
my.mat1

#Classify juvelines and adults (PCA + clustering)
my.mat2 = classifyCluster(data = crabdat)
my.mat2

#Calculate ogive
my.ogive1  = ogive(my.mat1)
my.ogive2  = ogive(my.mat2)
#Calculate ogive (Bayesian MCMClogit)
my.ogiveB1 = ogiveBayes(my.mat1)
my.ogiveB2 = ogiveBayes(my.mat2)


#Plot
plot(my.ogive1)
plot(my.ogiveB1)

plot(my.ogive2)
plot(my.ogiveB2)


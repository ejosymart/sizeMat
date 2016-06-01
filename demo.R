#install and load this packages
library(gdata)
library(rJava)
library(xlsx)
library(UsingR)

rm(list=ls())
source("MatureCode.R")

#Read Data
file = "crabdat"
crabdat = dataMature(file, ext = "txt")

#View data
headtail(crabdat)

#Calculate and plot
my.mat = matureFit(crabdat, n.iter = 100)
plot(my.mat)


#If you have a diferent extension in your data (csv, csv2, xlsx, delim),
#explore the argument "ext" in the function dataMature 



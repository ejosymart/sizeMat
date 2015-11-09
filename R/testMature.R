source("Mature.code.R")
data = read.table("CrabData.txt")
names(data) = c("x", "y")
matureFit(data) #select bounds




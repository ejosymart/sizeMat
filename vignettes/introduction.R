## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T)
library(ssmRG)

## ----echo=TRUE-----------------------------------------------------------
data(crabdata)

head(crabdata)

names(crabdata)

## ---- echo=TRUE----------------------------------------------------------
#For all the individuals
classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
varsex = "sex_category", selectSex = NULL, method = "ld")

#For males only
classify_data_males = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
varsex = "sex_category", selectSex = "m", method = "ld")


## ---- echo=TRUE, fig.width = 10, fig.height = 10, results='hide', warning=FALSE----
par(mfrow = c(2,2))
plot(classify_data)

plot(classify_data, xlab = "Carapace width (mm.)", ylab = "Chela heigth (mm)")

plot(classify_data, xlab = "Carapace width (mm.)", ylab = "Chela heigth (mm)", 
     col = c(2, 3), pch = c(5, 6))

plot(classify_data, xlab = "Carapace width (mm.)", ylab = "Chela heigth (mm)", 
     col = c(2, 3), pch = c(5, 6), lty_lines = c(1, 2), lwd_lines = c(1, 3), main = "Classification")

## ----echo=TRUE-----------------------------------------------------------
my_ogive = calculate_ogive(classify_data, method = "fq", niter = 1000)



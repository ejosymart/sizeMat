# ssmRG package: Size at Sexual Maturity based on Relative Growth---------------
#' @import MCMCpack
#' @import UsingR
#' @import biotools
#' @import matrixStats
#' @import MASS
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines plot
#' @importFrom stats binomial coef complete.cases cutree dist glm hclust prcomp predict
#' @importFrom utils data read.csv read.csv2 read.table
#'
#' @title Size at Sexual Maturity based on Relative Growth.
#'
#' @name ssmRG-package
#' @description  Package to Estimate size at sexual maturity from morphometic data based on relative growth.
#' @aliases ssmRG-package ssmRG.
#' @docType package
#' @references ssmRG: Size at Sexual Maturity based on Relative Growth (RJournal)
#' @keywords size, sexual - maturity, alometric, relative-growth.


NULL
#' Read data
#' 
#' Read a database with different extesions. The data base have only two variables, first variable is the independet variable and the second is the dependent variable.
#'
#' @param file The filename (extension: txt, csv, csv2)
#'
#' @return Database with x (independent) and y (dependent) variables.
#' @examples 
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data = read_data(my_file)
#' @export
read_data <- function(file){
  data <- strsplit(file, "\\.")  
  ext  <- data[[1]][2]
  
  data <- switch(ext,
                 "txt"   = read.table(file, header = TRUE, stringsAsFactors = FALSE),
                 "csv"   = read.csv(file),
                 "csv2"  = read.csv2(file))
  names(data) <- c("x", "y")
  
  return(data)
}


#' Classify mature
#' 
#' Classify in two groups (juvelines = 0 and adult = 1). The analisys is based on 
#' Principal Components Analisys  with the variables 
#' (x: independent variable, y: dependent variable) in log base, allowing to distinguish 
#' two groups would represent juveniles and adult.
#' The individuals are assigned to each group using a hierarchical classification 
#' procedure (hierarchical cluster). This method is based on establishing a predetermined 
#' number of groups (in this case, two) and assigning individuals to one of the groups according 
#' to their loads on the two axes of the PCA.
#' Using the results of the classification (PCA + cluster), a discriminant analysis (linear or quadratic) 
#' is conducted to obtain a discriminating function that permitted any individuals 
#' to be classified as a juvenile or an adult on the basis of the X and Y variables.
#'
#' @param data The database with two variables (x: independent variable, y: dependent variable)
#'
#' @return Database with x (independent), y (dependent) and mature classification
#' (juveniles = 0, adult = 1) variables.
#' @exportClass classify
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data)
#' classify_data
#' @export
classify_mature <- function(data){
  
  data <- data[complete.cases(data), ]
  data <- log(data)
  
  pca_classify    <- prcomp(data, 2)
  scores          <- pca_classify$x
  clusters        <- hclust(dist(scores, method = 'euclidean'), method = 'ward.D')
  mature_classify <- cutree(clusters, 2) - 1
  
  base            <- data.frame(data, mature_binom = mature_classify)
  
  test_cov_mat <- boxM(base[, c("x", "y")], base[, "mature_binom"])
  if(test_cov_mat$p.value > 0.05){
    dis_reg    <- lda(mature_binom ~ ., data = base)
  }else{
    dis_reg    <- qda(mature_binom ~ ., data = base)
  }

  mature       <- as.numeric(as.character(predict(dis_reg)$class))
  cat("number in juveline and adult group =", as.numeric(table(mature)[1]), ",", as.numeric(table(mature)[2]), "\n")
  
  data         <- data.frame(exp(base[, c("x", "y")]), mature)
  class(data)  <- "classify"
  
  return(data)
}


#' Plot classify data
#'
#' @param x Database (a list) with x (independent), y (dependent) and mature classification (juveniles = 0, adult = 1) variables.
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param col the colors for juveniles and adults group
#' @param pch the character indicating the type of plotting
#' @param \dots Additional arguments to be passed to methods
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data)
#' plot(classify_data, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5))
#' @export
#' @method plot classify
plot.classify <- function(x, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5), ...){
  
  data <- data.frame(do.call("cbind", x))
  
  juv  <- data[data$mature == 0, ]
  adt  <- data[data$mature == 1, ] 
  
  fit_juv <- glm(y ~ x, data = juv)
  fit_adt <- glm(y ~ x, data = adt)
  
  pch <- ifelse (data$mature == 0, pch[1], pch[2])
  col <- ifelse (data$mature == 0, col[1], col[2])
  plot(data$x, data$y, col = col, xlab = xlab, ylab = ylab, pch = pch)
  lines(juv$x, predict(fit_juv), col = col[1], lwd = 2)
  lines(adt$x, predict(fit_adt), col = col[2], lwd = 2)
  eq_juv  <- paste0("Y = ", round(as.numeric(coef(fit_juv)[1]), 2), " + ", round(as.numeric(coef(fit_juv)[2]),2), " *X", sep = "")
  eq_adt  <- paste0("Y = ", round(as.numeric(coef(fit_adt)[1]), 2), " + ", round(as.numeric(coef(fit_adt)[2]),2), " *X", sep = "")
  legend("topleft", c(paste("Juveniles: ", eq_juv), paste("Adults: ", eq_adt)), 
         bty = "n", pch = pch, col = col, cex = 0.8)
  return(invisible(NULL))
}


#' Calculate ogive
#' 
#' Estimate the size at 50\% maturity using a logistic regression
#' relating X variable and maturity stage (classifying individuals into juveniles 
#' or adults depending on their morphometry).
#'
#' @param data The database with the X, Y and mature stages (juvelines = 0, adults = 1)
#' @param methodReg The method to be applied, "fq"frecuentist GLM, or "bayes" bayesian GLM
#' (MCMClogit function).
#' @return List database with the parameters and a data.frame with the X, Y and mature stages
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @exportClass ogive
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data)
#' my_ogive = calculate_ogive(classify_data, methodReg = "fq")
#' @export
calculate_ogive <- function(data, methodReg = "fq"){
  
  input    <- data.frame(do.call("cbind", data))
  
  estimate <- switch(methodReg,
                     fq = .calculate_ogive_fq(input),
                     bayes = .calculate_ogive_bayes(input)) 
  
  out   <- data.frame(x = input$x, y = input$y, mature = input$mature, CIlower = estimate$lower, 
                        fitted = estimate$fitted, CIupper = estimate$upper)
  output     <- list(params = estimate$params, out = out)
  class(output) <- "ogive"
  
  return(out)
}


#' Plot maturity ogive
#'
#' @param x List database with the parameters and a data.frame with the X, Y and mature stages
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @param xlab1 a title for the x axis
#' @param ylab1 a title for the y axis
#' @param col the color for the logistic curve for the 50\% size at sexual maturity 
#' @param \dots Additional arguments to be passed to methods
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data)
#' my_ogive = calculate_ogive(classify_data, methodReg = "fq")
#' plot(my_ogive, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"))
#' @export
#' @method plot ogive
plot.ogive <- function(x, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"), ...){
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  
  plot(sort(unique(x_input)), m_p, xlab = xlab1, ylab = ylab1, pch = 19, col = "darkgrey", axes = F)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  
  wide50 <- - x$params[1] / x$params[2]
  wide95 <- (1 / x$params[2]) * log(1 / 0.05 - 1) - x$params[1] / x$params[2]
  
  lines(sort(x_input), sort(fit$fitted), col = col[1], lwd = 2)
  lines(sort(x_input), sort(fit$CIlower), col = col[1], lwd = 2, lty = 2)
  lines(sort(x_input), sort(fit$CIupper), col = col[1], lwd = 2, lty = 2)
  lines(c(wide50, wide50), c(-1, 0.5), col = col[2], lty = 2, lwd = 2)
  lines(c(-1, wide50), c(0.5, 0.5), col = col[2], lty = 2, lwd = 2)
  legend("topleft", as.expression(bquote(bold(CW[50] == .(round(wide50, 1))))), bty = "n")
  cat("Carapace width of 50% maturity =", round(wide50, 1), "\n")
  cat("Carapace width of 95% maturity =", round(wide95, 1), "\n")
  return(invisible(NULL))
}

# ssmRG package: Size at Sexual Maturity based on Relative Growth---------------
#' @import MCMCpack
#' @import UsingR
#' @import biotools
#' @import matrixStats
#' @import MASS
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines plot points hist par abline
#' @importFrom stats binomial quantile coef complete.cases cutree dist glm hclust prcomp predict
#' @importFrom utils data read.csv read.csv2 read.table
#'
#' @title Size at Sexual Maturity based on Relative Growth.
#'
#' @name ssmRG-package
#' @description  Package to estimate size at sexual maturity from morphometic data based on relative growth. 
#' Principal Components Analysis with two variables (x: independent variable, y: dependent variable), 
#' hierarchical clustering and linear or quadratic discriminant analysis 
#' are used to classify the individuals in two groups (juveniles or adults). 
#' Some basic plotting (classification and maturity ogive) are also provided.
#' @details Package: ssmRG
#' @details Type: Package
#' @aliases ssmRG-package ssmRG.
#' @docType package
#' @author Edgar Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @references Agostinho, C. S. (2000). Use of otoliths to estimate size at sexual maturity in fish. Brazilian Archives of Biology and Technology, 43(4).
#' @references Corgos, A., & Freire, J. (2006). Morphometric and gonad maturity in the spider crab Maja brachydactyla: a comparison of methods for estimating size at maturity in species with determinate growth. ICES Journal of Marine Science: Journal du Conseil, 63(5), 851-859.
#' @references Somerton, D. A. (1980). A computer technique for estimating the size of sexual maturity in crabs. Canadian Journal of Fisheries and Aquatic Sciences, 37(10), 1488-1494.
#' @keywords size, sexual - maturity, alometric, relative-growth.


NULL
#' Read data
#' 
#' Read a database with different extesions. The database has three variables, "x", "y" and "sex".
#' The variables "x" and "y" will be used in the regression. The variable "sex" 
#' is the sex category (\code{m} = male, \code{f} = female, \code{ind} = indeterminate sex.) 
#'
#' @param file the filename (extension: txt, csv, csv2).
#'
#' @return an object of class \code{mature} with x (independent), y (dependent) variables and sex category (\code{m} = male,
#' \code{f} = female, \code{ind} = indeterminate sex).
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
  names(data) <- c("x", "y", "sex")
  class(data) <- "mature"
  
  return(data)
}


#' Print mature database
#'
#' @param x object of class \code{mature}.
#' @param \dots Additional arguments to be passed
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' print(data)
#' @export
#' @method print mature
print.mature <- function(x, ...){
  data <- data.frame(do.call("cbind", x), stringsAsFactors = FALSE)
  cat("number of observations = ", dim(data)[1], "\n")
  cat("number of males = ", length(data$sex[data$sex == "m"]), "\n")
  cat("number of females = ", length(data$sex[data$sex == "f"]), "\n")
  cat("number of indeterminate sex = ", length(data$sex[data$sex == "ind"]), "\n")
}


#' Classify mature
#' 
#' Classify in two groups (juvelines = 0 and adult = 1). The analisys is based on 
#' Principal Components Analisys with the variables 
#' (x: independent variable, y: dependent variable) in log base, allowing to distinguish 
#' two groups that would represent juveniles and adult.
#' The individuals are assigned to each group using a hierarchical classification 
#' procedure (hierarchical cluster). This method is based on establishing a predetermined 
#' number of groups (in this case, two) and assigning individuals to one of the groups according 
#' to their loads on the two axes of the PCA.
#' Using the results of the classification (PCA + cluster), a discriminant analysis (linear or quadratic) 
#' is conducted to obtain a discriminating function that permitted any individuals 
#' to be classified as a juvenile or an adult on the basis of the X and Y variables.
#'
#' @param data object of class \code{mature} with three variables 
#' (x: independent variable, y: dependent variable, sex: sex category)
#' @param sex filter sex category, \code{"m"} = male, \code{"f"} = female, \code{"ind"} = indeterminate sex.
#' If \code{sex = "all"} all the individuals will be used in the analysis. 
#' @param method the discriminant analysis method, linear discriminant analysis \code{"ld"},
#' quadratic discriminant analysis \code{"qd"}. If \code{method} = \code{NULL}, ld or qd will be used 
#' in the discrimination analysis based on the test of homogeneity of covariance matrices.
#'
#' @return object of class \code{classify}, with x (independent), y (dependent) and mature classification
#' (juveniles = 0, adult = 1) variables.
#' @exportClass classify
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data, sex = "all", method = "ld")
#' classify_data
#' @export
classify_mature <- function(data, sex = "all", method = NULL){
  
  input <- data.frame(do.call("cbind", data), stringsAsFactors = FALSE)
  input$x <- as.numeric(input$x)
  input$y <- as.numeric(input$y)
  input <- input[complete.cases(input), ]
  
  if(sex == "all"){
     input <- input[, c("x", "y")]
     cat("all individuals were used in the analysis", "\n")
  }else{
    input <-input[input$sex == sex,]
    cat("only", paste0(unique(input$sex), "-sex", sep =""), "individuals were used in the analysis", "\n\n")
    input <- input[, c("x", "y")]
  }
  
  input <- log(input)
  
  pca_classify    <- prcomp(input, 2)
  scores          <- pca_classify$x
  clusters        <- hclust(dist(scores, method = 'euclidean'), method = 'ward.D')
  mature_classify <- cutree(clusters, 2) - 1
  
  base            <- data.frame(input, mature_binom = mature_classify)

  if(is.null(method)){
    test_cov_mat <- boxM(base[, c("x", "y")], base[, "mature_binom"])
    if(test_cov_mat$p.value > 0.05){
      dis_reg    <- lda(mature_binom ~ ., data = base)
    }else{
      dis_reg    <- qda(mature_binom ~ ., data = base)
    }
  } else {dis_reg  <- switch (method,
                              ld = lda(mature_binom ~ ., data = base),
                              qd = qda(mature_binom ~ ., data = base))
  }
  
  mature       <- as.numeric(as.character(predict(dis_reg)$class))
  cat("number in juveline group =", as.numeric(table(mature)[1]), "\n")
  cat("number in adult group =", as.numeric(table(mature)[2]), "\n")
  data         <- data.frame(exp(base[, c("x", "y")]), mature = mature)
  class(data)  <- "classify"
  
  return(data)
}


#' Plot classify data
#'
#' @param x object of class \code{classify}, with x (independent), y (dependent) and mature classification (juveniles = 0, adult = 1) variables.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col the colors for juveniles and adults group.
#' @param pch the character indicating the type of plotting.
#' @param \dots Additional arguments to be passed 
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data, sex = "all")
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
  plot(data$x, data$y, col = col, xlab = xlab, ylab = ylab, pch = pch, ...)
  lines(juv$x, predict(fit_juv), col = col[1], lwd = 2)
  lines(adt$x, predict(fit_adt), col = col[2], lwd = 2)
  eq_juv <- paste0("Y = ", round(as.numeric(coef(fit_juv)[1]), 2), " + ", round(as.numeric(coef(fit_juv)[2]),2), " *X", sep = "")
  eq_adt <- paste0("Y = ", round(as.numeric(coef(fit_adt)[1]), 2), " + ", round(as.numeric(coef(fit_adt)[2]),2), " *X", sep = "")
  legend("topleft", c(paste("Juveniles: ", eq_juv), paste("Adults: ", eq_adt)), 
         bty = "n", pch = pch, col = col, cex = 0.8)
  return(invisible(NULL))
}


#' Calculate ogive
#' 
#' Estimate the size at maturity using a logistic regression
#' relating X variable and maturity stage (classifying individuals as juveniles 
#' or adults depending on their morphometry).
#'
#' @param data object of class \code{classify} with the X, Y and mature classification (juvelines = 0, adults = 1).
#' @param method the method to be applied, \code{"fq"} frecuentist GLM, or \code{"bayes"} bayesian GLM
#' (MCMClogit function).
#' @param niter number of iterations (bootstrap resampling).
#' @param seed a single value, interpreted as an integer.
#' @return object of class \code{ogive} with the parameters, the L50 and data.frame with the X, Y and mature classification
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @exportClass ogive
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data, sex = "all")
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' my_ogive$A_boot
#' my_ogive$B_boot
#' my_ogive$L50_boot
#' my_ogive$out
#' @export
calculate_ogive <- function(data, method = "fq", niter = 1000, seed = 70387){
  
  input    <- data.frame(do.call("cbind", data))
  
  estimate <- switch(method,
                     fq = .calculate_ogive_fq(input,  niter = niter, seed = seed),
                     bayes = .calculate_ogive_bayes(input,  niter = niter, seed = seed)) 
  
  out   <- data.frame(x = input$x, y = input$y, mature = input$mature, CIlower = estimate$lower, 
                      fitted = estimate$fitted, CIupper = estimate$upper)
  
  output <- list(A_boot   = estimate$parameters_A, 
                 B_boot   = estimate$parameters_B,
                 L50_boot = estimate$L50,
                 out      = out)
  
  class(output) <- "ogive"
  
  return(output)
}


#' Print maturity ogive
#'
#' @param x object of class \code{ogive} with the ogive parameters and a data.frame with the X, Y and mature classification
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @param probs numeric vector of probabilities with values in \code{[0,1]}. 
#' @param \dots Additional arguments to be passed
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data, sex = "all")
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' print(my_ogive)
#' @export
#' @method print ogive
print.ogive <- function(x, probs = c(0.05, 0.5, 0.95), ...){
  L50     <- quantile(x$L50_boot, probs = probs, na.rm = TRUE)
  cat("formula: Y = 1/1+exp-(A + B*X)", "\n\n")
  tab <- matrix(as.numeric(c(L50)), nrow = 1, ncol = 3, byrow = TRUE)
  colnames(tab) <- c("0.05 %", "0.50 %", "0.95 %")
  rownames(tab) <- c("L50")
  tab <- as.table(tab)
  return(tab)
}


#' Plot maturity ogive
#'
#' @param x object of class \code{ogive} with the ogive parameters and a data.frame with the X, Y and mature classification
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @param xlab1 a title for the x axis.
#' @param ylab1 a title for the y axis.
#' @param col the color for the logistic curve and for the L50\% size at sexual maturity.
#' @param \dots Additional arguments to be passed to methods
#' @examples
#' my_file = system.file("extdata", "crabdat.txt", package = "ssmRG")
#' data    = read_data(my_file)
#' classify_data = classify_mature(data, sex = "all")
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' print(my_ogive)
#' plot(my_ogive, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"))
#' @export
#' @method plot ogive
plot.ogive <- function(x, xlab1 = "X", ylab1 = "Proportion mature", col = c("blue", "red"), ...){
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  wide    <- quantile(x$L50_boot, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
  
  hist(x$A_boot, main = "", xlab = "A")
  abline(v = as.numeric(quantile(x$A_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = 2, col = col[1])
  abline(v = c(as.numeric(quantile(x$A_boot, probs = c(0.05, 0.95), na.rm = TRUE))), 
         lty = 2, col = col[1])

  hist(x$B_boot, main = "", xlab = "B")
  abline(v = as.numeric(quantile(x$B_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = 2, col = col[1])
  abline(v = c(as.numeric(quantile(x$B_boot, probs = c(0.05, 0.95), na.rm = TRUE))), 
         lty = 2, col = col[1])

  hist(x$L50_boot, main = "", xlab = "Length of 50% maturity")
  abline(v = as.numeric(quantile(x$L50_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = 2, col = col[1])
  abline(v = c(as.numeric(quantile(x$L50_boot, probs = c(0.05, 0.95), na.rm = TRUE))), 
         lty = 2, col = col[1])

  plot(sort(unique(x_input)), m_p, xlab = xlab1, ylab = ylab1, pch = 19, col = "darkgrey", axes = F, ...)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  lines(sort(x_input), sort(fit$fitted), col = col[1], lwd = 2)
  lines(sort(x_input), sort(fit$CIlower), col = col[1], lwd = 2, lty = 2)
  lines(sort(x_input), sort(fit$CIupper), col = col[1], lwd = 2, lty = 2)
  lines(c(wide[2], wide[2]), c(-1, 0.5), col = col[2], lty = 2, lwd = 2)
  lines(c(-1, wide[2]), c(0.5, 0.5), col = col[2], lty = 2, lwd = 2)
  points(wide[2], 0.5, pch = 19, col = col[2], cex = 1.5)
  legend("topleft", as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))), bty = "n")
  cat("Length of 50% maturity =", round(wide[2], 1), "\n")
  cat("Confidence intervals =", round(wide[1], 1), "-",round(wide[3], 1) ,  "\n")
  
  return(invisible(NULL))
}

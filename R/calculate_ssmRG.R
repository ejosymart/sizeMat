# ssmRG package: Size at Sexual Maturity based on Relative Growth---------------
#' @importFrom MCMCpack MCMClogit
#' @importFrom matrixStats rowQuantiles
#' @importFrom MASS lda qda
#' @importFrom biotools boxM
#' @importFrom grDevices colors
#' @importFrom tools file_ext
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
#' Classify mature
#' 
#' Classify te individuals in two groups (0: juvelines, 1: adults).
#' @param data data.frame with alometric variables and sex category (male, female). 
#' If sex category has NA's, the row will be filter.
#' @param varnames the name of two allometric variables to be used for analysis.
#' @param varsex the name of the variable containing sex information.
#' @param useSex sex category to be used for analysis. If \code{useSex = NULL} all the individuals will be used in the analysis. 
#' @param method the discriminant analysis method, linear discriminant analysis \code{"ld"},
#' quadratic discriminant analysis \code{"qd"}. If \code{method} = \code{NULL}, ld or qd will be used 
#' in the discrimination analysis based on the test of homogeneity of covariance matrices.
#' @return object of class 'classify', with x (independent), y (dependent) and mature classification
#' (juveniles = 0, adult = 1) variables.
#' @details Classify in two groups (juvelines = 0 and adult = 1). The analisys is based on 
#' Principal Components Analisys with the variables (x: independent variable, y: dependent variable) in log base, 
#' allowing to distinguish two groups that would represent juveniles and adult.
#' The individuals are assigned to each group using a hierarchical classification 
#' procedure (hierarchical cluster). This method is based on establishing a predetermined 
#' number of groups (in this case, two) and assigning individuals to one of the groups according 
#' to their loads on the two axes of the PCA.
#' Using the results of the classification (PCA + cluster), a discriminant analysis (linear or quadratic) 
#' is conducted to obtain a discriminating function that permitted any individuals 
#' to be classified as a juvenile or an adult on the basis of the X and Y variables.
#' @exportClass classify
#' @examples
#' data(crabdata)
#' classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
#' varsex = "sex_category", useSex = NULL, method = "ld")
#' classify_data
#' @export
classify_mature <- function(data, varnames = c("x", "y"), varsex = "sex", 
                            useSex = NULL, method = NULL) {
  if(length(varnames) != 2) stop("You must provide two variables only.")
  if(!all(varnames %in% names(data))) stop("'varnames' have not been found in data.")
  input <- data[, c(varnames, varsex)]
  names(input)  <-  c("x", "y", "sex")
  input <- input[complete.cases(input), ] 
  
  if(is.null(useSex)) {
     input <- input[, c("x", "y")]
     cat("all individuals were used in the analysis", "\n\n")
  } else {
    input <-input[which(input$sex == useSex),]
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
    test_cov_mat <- biotools::boxM(base[, c("x", "y")], base[, "mature_binom"])
    if(test_cov_mat$p.value > 0.05){
      dis_reg    <- MASS::lda(mature_binom ~ ., data = base)
    }else{
      dis_reg    <- MASS::qda(mature_binom ~ ., data = base)
    }
  } else {dis_reg  <- switch (method,
                              ld = MASS::lda(mature_binom ~ ., data = base),
                              qd = MASS::qda(mature_binom ~ ., data = base))
  }
  
  mature       <- as.numeric(as.character(predict(dis_reg)$class))
  cat("number in juveline group =", as.numeric(table(mature)[1]), "\n\n")
  cat("number in adult group =", as.numeric(table(mature)[2]), "\n\n")
  data         <- data.frame(exp(base[, c("x", "y")]), mature = mature)
  
  class(data)  <- c("classify", class(data))
  
  return(data)
}


#' Plot classify data
#'
#' @param x an object of class 'classify', with x (independent), y (dependent) and mature classification (juveniles = 0, adult = 1) variables.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col the colors for juveniles and adults group.
#' @param pch the character indicating the type of plotting.
#' @param lty the line type in the regression.
#' @param lwd the line width in the regression.
#' @param cex character expansion in the regression.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' data(crabdata)
#' classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
#' varsex = "sex_category", useSex = NULL, method = "ld")
#' classify_data
#' plot(classify_data, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5))
#' @export
#' @method plot classify
plot.classify <- function(x, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5), 
                          lty = c(1, 1), lwd = c(1, 1), cex = c(1, 1), ...){
  
  data <- data.frame(do.call("cbind", x))
  juv  <- data[data$mature == 0, ]
  adt  <- data[data$mature == 1, ]
  
  fit_juv <- glm(y ~ x, data = juv)
  fit_adt <- glm(y ~ x, data = adt)
  
  PCH <- ifelse (data$mature == 0, pch[1], pch[2])
  COL <- ifelse (data$mature == 0, col[1], col[2])
  
  if(length(col) < 2) stop('col argument must have 2 values. The colors could be the same')
  if(length(col) > 2) warning('col: only the first two colors will be used in the plot')
  if(length(pch) < 2) stop('pch argument must have 2 values. The plotting character could be the same')
  if(length(pch) > 2) warning('pch: only the first two plotting character will be used in the plot')
  if(length(lty) < 2) stop('lty argument must have 2 values. The line type could be the same')
  if(length(lty) > 2) warning('lty: only the first two line type  will be used in the plot')
  if(length(lwd) < 2) stop('lwd argument must have 2 values. The line width could be the same')
  if(length(lwd) > 2) warning('lwd: only the first two line width will be used in the plot')
  if(length(cex) < 2) stop('cex argument must have 2 values. The character expansion could be the same')
  if(length(cex) > 2) warning('cex: only the first two character expansion will be used in the plot')
  
  plot(data$x, data$y, col = COL, xlab = xlab, ylab = ylab, pch = PCH, 
       lty = lty, lwd = lwd, cex = cex, ...)
  lines(juv$x, predict(fit_juv), col = COL[1], lwd = 2)
  lines(adt$x, predict(fit_adt), col = COL[2], lwd = 2)
  eq_juv <- paste0("Y = ", round(as.numeric(coef(fit_juv)[1]), 2), " + ", round(as.numeric(coef(fit_juv)[2]),2), " *X", sep = "")
  eq_adt <- paste0("Y = ", round(as.numeric(coef(fit_adt)[1]), 2), " + ", round(as.numeric(coef(fit_adt)[2]),2), " *X", sep = "")
  legend("topleft", c(paste("Juveniles: ", eq_juv), paste("Adults: ", eq_adt)), 
         bty = "n", pch = PCH, col = COL, cex = 0.8)
  return(invisible(NULL))
}


#' Calculate ogive
#' 
#' Estimate the size at sexual maturity (L50).
#'
#' @param data an object of class 'classify' with the X, Y and mature classification (juvelines = 0, adults = 1).
#' @param method the method to be applied, \code{"fq"} frecuentist GLM, or \code{"bayes"} bayesian GLM (MCMClogit function).
#' @param niter number of iterations (bootstrap resampling).
#' @param seed a single value, interpreted as an integer.
#' @return object of class 'ogive'.
#' 
#' \code{model} the summary model.
#' 
#' \code{A_boot} the 'n iter'values of parameter A.
#' 
#' \code{B_boot} the 'n iter'values of parameter B.
#' 
#' \code{L50} the 'n iter'values of parameter L50 (size at sexual maturity).
#' 
#' \code{out} a dataframe with the X, Y and mature classification variables. 
#' Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @details Estimate the size at sexual maturity using a logit regression with X variable and maturity stage (two categories: juveniles 
#' and adults).
#' @exportClass ogive
#' @examples
#' data(crabdata)
#' classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
#' varsex = "sex_category", useSex = NULL, method = "ld")
#' classify_data
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' my_ogive$A_boot
#' my_ogive$B_boot
#' my_ogive$L50_boot
#' my_ogive$out
#' @export
calculate_ogive <- function(data, method = "fq", niter = 1000, seed = 70387){
  
  estimate <- switch(method,
                     fq = .calculate_ogive_fq(data,  niter = niter, seed = seed),
                     bayes = .calculate_ogive_bayes(data,  niter = niter, seed = seed)) 
  
  out   <- data.frame(x = data$x, y = data$y, mature = data$mature, CIlower = estimate$lower, 
                      fitted = estimate$fitted, CIupper = estimate$upper)
  
  output <- list(model    = estimate$model,
                 A_boot   = estimate$parameters_A, 
                 B_boot   = estimate$parameters_B,
                 L50_boot = estimate$L50,
                 out      = out)
  
  class(output) <- c("ogive", class(output))
  
  return(output)
}


#' Print maturity ogive
#'
#' @param x object of class 'ogive' with the ogive parameters and a data.frame with the X, Y and mature classification
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @param probs numeric vector of probabilities with values in \code{[0,1]}. 
#' @param \dots Additional arguments to be passed
#' @examples
#' data(crabdata)
#' classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
#' varsex = "sex_category", useSex = NULL, method = "ld")
#' classify_data
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' print(my_ogive)
#' @export
#' @method print ogive
print.ogive <- function(x, probs = c(0.025, 0.5, 0.975), ...){
  L50     <- quantile(x$L50_boot, probs = probs, na.rm = TRUE)
  cat("formula: Y = 1/1+exp-(A + B*X)", "\n\n")
  tab <- matrix(as.numeric(c(L50)), nrow = 1, ncol = 3, byrow = TRUE)
  colnames(tab) <- c("0.025 %", "Median", "0.975 %")
  rownames(tab) <- c("L50")
  tab <- as.table(tab)
  return(tab)
}


#' Plot maturity ogive
#'
#' @param x object of class 'ogive' with the ogive parameters and a data.frame with the X, Y and mature classification
#' variables. Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col color for the logistic curve and for the L50\% size at sexual maturity.
#' @param lwd line with for drawing fitted values and confidence intervals.
#' @param lty line type line type for drawing fitted values and confidence intervals
#' @param vline_hist color of the vertival lines in the histogram. The lines represent the 
#' confidence intervals and the median.
#' @param lwd_hist line with for the vertical line in the histogram.
#' @param lty_hist line type for the vertical line in the histogram.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' data(crabdata)
#' classify_data = classify_mature(crabdata, varnames = c("carapace_width", "chela_heigth"), 
#' varsex = "sex_category", useSex = NULL, method = "ld")
#' my_ogive = calculate_ogive(classify_data, method = "fq")
#' plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"))
#' @export
#' @method plot ogive
plot.ogive <- function(x, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"), 
                       lwd = 2, lty = 2, vline_hist = "black", lwd_hist = 2, lty_hist = 2, ...){
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  wide    <- quantile(x$L50_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  # figure 1
  hist(x$A_boot, main = "", xlab = "A")
  abline(v = as.numeric(quantile(x$A_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = lwd_hist, col = vline_hist)
  abline(v = c(as.numeric(quantile(x$A_boot, probs = c(0.025, 0.975), na.rm = TRUE))), 
         lty = lty_hist, col = vline_hist)
  box()

  # figure 2
  hist(x$B_boot, main = "", xlab = "B")
  abline(v = as.numeric(quantile(x$B_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = lwd_hist, col = vline_hist)
  abline(v = c(as.numeric(quantile(x$B_boot, probs = c(0.025, 0.975), na.rm = TRUE))), 
         lty = lty_hist, col = vline_hist)
  box()

  # figure 3
  hist(x$L50_boot, main = "", xlab = "Length of 50% maturity")
  abline(v = as.numeric(quantile(x$L50_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = lwd_hist, col = vline_hist)
  abline(v = c(as.numeric(quantile(x$L50_boot, probs = c(0.025, 0.975), na.rm = TRUE))), 
         lty = lty_hist, col = vline_hist)
  box()
  
  # figure 4
  if(length(col) < 2) stop('col argument must have 2 values. The colors could be the same')
  if(length(col) > 2) warning('col: only the first two colors will be used in the plot')
  plot(sort(unique(x_input)), m_p, xlab = xlab, ylab = ylab, pch = 19, col = "darkgrey", ...)
  lines(sort(x_input), sort(fit$fitted), col = col[1], lwd = lwd)
  lines(sort(x_input), sort(fit$CIlower), col = col[1], lwd = lwd, lty = lty)
  lines(sort(x_input), sort(fit$CIupper), col = col[1], lwd = lwd, lty = lty)
  lines(c(wide[2], wide[2]), c(-1, 0.5), col = col[2], lwd = lwd, lty = lty)
  lines(c(-1, wide[2]), c(0.5, 0.5), col = col[2], lwd = lwd, lty = lty)
  points(wide[2], 0.5, pch = 19, col = col[2], cex = 1.5)
  legend("topleft", as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))), bty = "n")
  cat("Length of 50% maturity =", round(wide[2], 1), "\n")
  cat("Confidence intervals =", round(wide[1], 1), "-",round(wide[3], 1) ,  "\n")
  
  return(invisible(NULL))
}

# sizeMat package: Estimate Size at Sexual Maturity ---------------

#' @importFrom MCMCpack MCMClogit
#' @importFrom matrixStats rowQuantiles
#' @importFrom MASS lda qda
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines plot points hist par abline
#' @importFrom stats binomial quantile coef complete.cases cutree dist glm hclust prcomp predict cov pchisq
#' @importFrom utils data read.csv installed.packages
#'
#' @title Size at Sexual Maturity.
#'
#' @description This package allows to estimate morphometric and gonadal size at sexual maturity for organisms, usually fish and invertebrates. The size at sexual maturity is defined as the length at which a randomly chosen specimen has a 50% chance of being mature. The package includes methods for classification based on relative growth (PCA, hierarchical clustering, discriminant analysis), logit regression (frequentist or bayesian), parameters estimation and some basic plots.
#' @name sizeMat-package
#' @aliases sizeMat-package sizeMat
#' @docType package
#' @author Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @details Package: sizeMat
#' @details Type: Package
#' @details The Morphometric and Gonadal size at sexual maturity are estimating using differents functions (process).
#' 
#' 1) The estimation of the Morphometric Size at Sexual Maturity involves two processes:
#' 
#' 1.1) A Principal Components Analisys is conducted with two allometric variables (x: independent variable, y: dependent variable) in log base, allowing to distinguish 
#' two groups that would represent juveniles and adult. The individuals are assigned to each group using a hierarchical classification procedure (hierarchical cluster). 
#' This method is based on establishing a predetermined number of groups (in this case, two) and assigning individuals to one of the groups according to 
#' their loads on the two axes of the PCA (Corgos & Freire, 2006). Using the results of the classification (PCA + cluster), a discriminant analysis (linear or quadratic) 
#' is carried out to obtain a discriminating function that permitted any individuals to be classified as a juvenile or an adult on the basis of the X and Y 
#' allometric variables.
#' 
#' 1.2) After classification, the logistic approach is used. The morphometric size at 50\% maturity (\eqn{L_50}) is estimated as the length at
#' which a randomly chosen specimen has a 50\% chance of being mature (Somerton  1980, Roa  et al. 1999, Corgos & Freire 2006). 
#' In the regression analysis, \eqn{X} (e.g: carapace width) is considered the explanatory variable and the classification \eqn{CS} 
#' (juvelines: 0, adults: 1) is considered the response variable (binomial). 
#' 
#' The variables are fitted to a logistic function with the form: 
#' 
#' \deqn{P_CS = 1 / [1+e^-(beta_0 + beta_1*X)]}
#' 
#' where:
#' 
#' \eqn{P_CL} is the probability of an individual of being mature at a determinate \eqn{X} length.
#'  
#' \eqn{beta_0} (intercept) and \eqn{beta_1} (slope) are parameters estimated. 
#' 
#' The (\eqn{L_50}) is calculated as: 
#' 
#' \deqn{L_50 = -beta_0 / beta_1}
#' 
#' Some basic plotting (classification, \eqn{beta_0}, \eqn{beta_1} and \eqn{L_50} histogram, and maturity ogive) 
#' are also provided.
#' 
#' 
#' 2) The estimation of Gonadal Size at Sexual Maturity use the logistic approach only. 
#' 
#' To estimate gonadal size at sexual maturity, the database must contains the stage of sexual 
#' maturity and at least one allometric variable (e.g: total length, fork length, carapace width).
#' The stage of sexual maturity is refered to the gonadal maturarion stages (e.g: I, II, III, IV or 0, 1, etc).
#' 
#' So, in the regression analysis, the allometric variable (e.g: total length) is considered the
#' explanatory variable and the stage of sexual maturity (inmature: 0, mature: 1) 
#' is considered the response variable (binomial). The regression  analysis is performed 
#' in the same way as the morphometric size at sexual maturity.
#' 
#' @references Agostinho, C. S. (2000). Use of otoliths to estimate size at sexual maturity in fish. Brazilian Archives of Biology and Technology, 43(4).
#' @references Corgos, A. & Freire, J. (2006). Morphometric and gonad maturity in the spider crab Maja brachydactyla: a comparison of methods for estimating size at maturity in species with determinate growth. ICES Journal of Marine Science: Journal du Conseil, 63(5), 851-859.
#' @references Roa, R., Ernst, B. & Tapia, F. (1999). Estimation of size at sexual maturity: an evaluation of analytical and resampling procedures. Fishery Bulletin, 97(3), 570-580.
#' @references Somerton, D. A. (1980). A computer technique for estimating the size of sexual maturity in crabs. Canadian Journal of Fisheries and Aquatic Sciences, 37(10), 1488-1494.
#' @keywords morphometric, size, sexual-maturity, alometric, relative-growth.
#' @examples
#' #See examples for functions morph_mature() and gonad_mature().

NULL
#' Classify mature
#' 
#' Classify te individuals in two groups (0: juvelines, 1: adults) based on relative growth.
#' @param data data.frame with allometric variables and sex category (male, female). 
#' If sex category contains NA's, that row will be filtered.
#' @param varNames the name of two allometric variables to be used for analysis.
#' @param varSex the name of the variable containing sex information.
#' @param selectSex sex category to be used for analysis. If \code{selectSex = NULL} all the individuals will be used in the analysis. 
#' @param method a character string indicating the discriminant analysis method, linear discriminant analysis \code{"ld"},
#' quadratic discriminant analysis \code{"qd"}.
#' We suggest begin the analysis using the \code{method = "ld"}.
#' @return A data.frame of class 'classify', with x (independent), y (dependent) and classification of maturity
#' (juveniles = 0, adult = 1) variables.
#' @details Classify the individuals in two groups (juvelines = 0 and adult = 1).
#' 
#' A Principal Components Analisys was conducted with two allometric variables (x: independent variable, y: dependent variable) 
#' in log base, allowing to distinguish two groups that would represent juveniles and adult.
#' The individuals are assigned to each group using a hierarchical classification procedure 
#' (hierarchical cluster with agglomeration method: "Ward.D" and the distance measure: "euclidean"). 
#' This method is based on establishing a predetermined number of groups (in this case, two) and assigning individuals 
#' to one of the groups according to their loads on the two axes of the PCA (Corgos & Freire, 2006). 
#' 
#' Using the results of the classification (PCA + cluster), a discriminant analysis (linear or quadratic) is conducted 
#' to obtain a discriminating function that permitted any individuals to be classified as a 
#' juvenile or an adult on the basis of the X and Y allometric variables.
#' @exportClass classify
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#' 
#' classify_data}
#' @export
classify_mature <- function(data, varNames = c("x", "y"), varSex = "sex", 
                            selectSex = NULL, method = "ld") {
  
  if(length(varNames) != 2) stop("You must provide two variables only.")
  if(!all(varNames %in% names(data))) stop("'varNames' have not been found in data.")
  if(varNames[1] == varNames[2]) stop("'varNames' must have different names")
  if(!all(varSex %in% names(data))) stop("'varSex' have not been found in data.")
  
  input <- data[, c(varNames, varSex)]
  names(input)  <-  c("x", "y", "sex")
  input <- input[complete.cases(input), ] 
  
  if(is.null(selectSex)) {
     input <- input[, c("x", "y")]
     cat("all individuals were used in the analysis", "\n\n")
  } else {
    input <-input[which(input$sex == selectSex),]
    cat("only", paste0(unique(input$sex), "-sex", sep =""), "were used in the analysis", "\n\n")
    input <- input[, c("x", "y")]
  }
  
  input <- log(input)
  
  pca_classify    <- prcomp(input, 2)
  scores          <- pca_classify$x
  clusters        <- hclust(dist(scores, method = 'euclidean'), method = 'ward.D')
  mature_classify <- cutree(clusters, 2) - 1
  
  base            <- data.frame(input, mature_binom = mature_classify)

  dis_reg  <- switch (method, 
                      ld = MASS::lda(mature_binom ~ ., data = base), 
                      qd = MASS::qda(mature_binom ~ ., data = base))
  
  mature       <- as.numeric(as.character(predict(dis_reg)$class))
  data         <- data.frame(exp(base[, c("x", "y")]), mature = mature)
    
  class(data)  <- c("classify", class(data))
  
  return(data)
}


#' Print method for classify class
#'
#' @param x an object of class 'classify' with the allometric variables ("X", "Y") and classification of maturity (juvelines = 0, adults = 1).
#' @param \dots Additional arguments to the print method.
#' @return The number of juveniles and adults. Also shows the regression analysis for juveniles and adults.
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#'
#' print(classify_data)
#' }
#' @export
#' @method print classify
print.classify <- function(x, ...){
  if (!inherits(x, "classify"))
    stop("Use only with 'classify' objects")
  data <- x
  cat("Number in juveline group =", as.numeric(table(data$mature)[1]), "\n\n")
  cat("Number in adult group =", as.numeric(table(data$mature)[2]), "\n\n")

  juv  <- data[data$mature == 0, ]
  adt  <- data[data$mature == 1, ]
  
  fit_juv <- glm(y ~ x, data = juv)
  fit_adt <- glm(y ~ x, data = adt)
  cat("--------------------------------------------------------", "\n")
  cat("1) Linear regression for juveniles", "\n")
  print(summary(fit_juv))
  cat("--------------------------------------------------------", "\n")
  cat("2) Linear regression for adults", "\n")
  print(summary(fit_adt))
  
  return(invisible())
}


#' Plot method for classify class
#'
#' @param x an object of class 'classify' with the allometric variables ("X", "Y") and classification of maturity (juvelines = 0, adults = 1).
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col the colors for juveniles and adults group.
#' @param pch the character indicating the type of plotting.
#' @param lty_lines the line type in the regression.
#' @param lwd_lines the line width in the regression.
#' @param cex character expansion in the regression.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#' 
#' ## Showing different plots
#' plot(classify_data, xlab = "X")
#' 
#' plot(classify_data, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5), cex = c(1, 3))
#' 
#' plot(classify_data, xlab = "Carapace width (mm.)", ylab = "Y", col = c(1, 2), 
#' pch = c(4, 5), cex = c(1, 3), lwd_lines = c(1, 3))
#' 
#' plot(classify_data, xlab = "Carapace width (mm.)", ylab = "Y", col = c(1, 2), 
#' pch = c(4, 5), cex = c(1, 3), lwd_lines = c(1, 3), main = "Classification")
#' }
#' @export
#' @method plot classify
plot.classify <- function(x, xlab = "X", ylab = "Y", col = c(1, 2), pch = c(4, 5), 
                          cex = c(1, 1), lty_lines = c(1, 1), lwd_lines = c(1, 1), ...){
  
  if (!inherits(x, "classify"))
    stop("Use only with 'classify' objects")
  
  data <- x
  juv  <- data[data$mature == 0, ]
  adt  <- data[data$mature == 1, ]
  
  fit_juv <- glm(y ~ x, data = juv)
  fit_adt <- glm(y ~ x, data = adt)
  
  if(length(col) < 2) stop('col argument must have 2 values. The colors could be the same')
  if(length(col) > 2) warning('col: only the first two colors will be used in the plot')
  if(length(pch) < 2) stop('pch argument must have 2 values. The plotting character could be the same')
  if(length(pch) > 2) warning('pch: only the first two plotting character will be used in the plot')
  if(length(cex) < 2) stop('cex argument must have 2 values. The character expansion could be the same')
  if(length(cex) > 2) warning('cex: only the first two character expansion will be used in the plot')
  if(length(lty_lines) < 2) stop('lty argument must have 2 values. The line type could be the same')
  if(length(lty_lines) > 2) warning('lty: only the first two line type  will be used in the plot')
  if(length(lwd_lines) < 2) stop('lwd argument must have 2 values. The line width could be the same')
  if(length(lwd_lines) > 2) warning('lwd: only the first two line width will be used in the plot')
  
  PCH <- ifelse (data$mature == 0, pch[1], pch[2])
  COL <- ifelse (data$mature == 0, col[1], col[2])
  CEX <- ifelse (data$mature == 0, cex[1], cex[2])
  LTY <- ifelse (data$mature == 0, lty_lines[1], lty_lines[2])
  LWD <- ifelse (data$mature == 0, lwd_lines[1], lwd_lines[2])
  
  plot(data$x, data$y, type = "p", col = COL, xlab = xlab, ylab = ylab, pch = PCH, cex = CEX, ...)
  lines(juv$x, predict(fit_juv), col = COL[1], lwd = LWD[1], lty = LTY[1])
  lines(adt$x, predict(fit_adt), col = COL[2], lwd = LWD[2], lty = LTY[2])
  eq_juv <- paste0("Y = ", round(as.numeric(coef(fit_juv)[1]), 2), " + ", round(as.numeric(coef(fit_juv)[2]),2), " *X", sep = "")
  eq_adt <- paste0("Y = ", round(as.numeric(coef(fit_adt)[1]), 2), " + ", round(as.numeric(coef(fit_adt)[2]),2), " *X", sep = "")
  legend("topleft", c(paste("Juveniles: ", eq_juv), paste("Adults: ", eq_adt)), 
         bty = "n", pch = PCH, col = COL, cex = 0.8)
  return(invisible(NULL))
}


#' Calculate morphometric mature
#' 
#' Estimate the morphometric size at sexual maturity (\eqn{L50}).
#'
#' @param data an object of class 'classify' with the allometric variables (X", "Y") and classification of maturity (juvelines = 0, adults = 1).
#' @param method a character string indicating the method to be applied, \code{"fq"} frecuentist GLM, or \code{"bayes"} bayesian GLM (MCMClogit function).
#' @param niter number of iterations (bootstrap resampling).
#' @param seed a single value, interpreted as an integer.
#' @return An object of class 'morphMat'.
#' 
#' \code{model} the summary statistics of the model.
#' 
#' \code{A_boot} the 'n iter' values of parameter A.
#' 
#' \code{B_boot} the 'n iter' values of parameter B.
#' 
#' \code{L50} the 'n iter' values of parameter L50 (morphometric size at sexual maturity).
#' 
#' \code{out} a dataframe with the allometric variables "X" and "Y", classification of maturity, the fitted values for  
#' logistic regression and confidence intervals (95\%). Also the summary statistics of the model is provided.
#' @details Estimate the morphometric size at sexual maturity using a logit regression with X variable 
#' and maturity classification (two categories: juveniles and adults). 
#' 
#' The function requires an object of class "classify" with the X, Y (allometric variables) and classification of maturity (juvelines = 0, adults = 1).
#' 
#' The argument `method` requires a character string indicanting which regression will be used for the test.
#' If `method = "fq"` the regression is based on GLM (frequentist) and if `method = "bayes"` a sample from the posterior distribution 
#' of a logistic regression model using a random walk Metropolis algorithm is generated (see MCMClogit function).
#' 
#' The argument `niter` requires a number. For the GLM regression (`method = "fq"`), a non-parametric bootstrap method consists
#' in generate B bootstrap samples, by resampling with replacement the original data. Then all statistics for each parameter 
#' can be calculated from each bootstrap sample (median and confidence intervals). 
#' For the `method = "bayes"`, the argument `niter` is related to the number of Metropolis iterations for the sampler.
#' @exportClass morphMat
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#' 
#' my_mature = morph_mature(classify_data, method = "fq")
#' 
#' ## \eqn{niter} parameters:
#' my_mature$A_boot
#' my_mature$B_boot
#' my_mature$L50_boot
#' my_mature$out
#' }
#' @export
morph_mature <- function(data, method = "fq", niter = 999, seed = 70387){
  
  if (!inherits(data, "classify"))
    stop("Use only with 'classify' objects")
  estimate <- switch(method,
                     fq = .morph_mature_fq(data = data,  niter = niter, seed = seed),
                     bayes = .morph_mature_bayes(data = data,  niter = niter, seed = seed))
  
  out   <- data.frame(x = data$x, y = data$y, mature = data$mature, CIlower = estimate$lower, 
                      fitted = estimate$fitted, CIupper = estimate$upper)
  
  output <- list(model    = estimate$model,
                 A_boot   = estimate$parameters_A, 
                 B_boot   = estimate$parameters_B,
                 L50_boot = estimate$L50,
                 out      = out)
  
  class(output) <- c("morphMat", class(output))
  
  return(output)
}


#' Print method for morphMat (morphometric size at sexual maturity)
#'
#' @param x object of class 'morphMat' with the mature parameters and a data.frame with the allometric variables ("X", "Y") 
#' and classification of maturity. Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @param \dots Additional arguments to the print method.
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#' 
#' my_mature = morph_mature(classify_data, method = "fq")
#' 
#' print(my_mature)
#' }
#' @export
#' @method print morphMat
print.morphMat <- function(x, ...){
  if (!inherits(x, "morphMat"))
    stop("Use only with 'morphMat' objects")
  cat("formula: Y = 1/1+exp-(A + B*X)", "\n\n")
  
  A_b   <- quantile(x$A_boot, probs = 0.5, na.rm = TRUE)
  B_b   <- quantile(x$B_boot, probs = 0.5, na.rm = TRUE)
  L50_b <- quantile(x$L50_boot, probs = 0.5, na.rm = TRUE)
  
  if(is.null(coef(x$model))){
    tab <- matrix(as.numeric(c(A_b, B_b, L50_b)), 
                  nrow = 3, ncol = 1, byrow = TRUE)
    colnames(tab) <- c("Bootstrap (Median)")
    rownames(tab) <- c("A", "B", "L50")
    tab <- as.table(tab)
    return(tab)
  }else{
    A_or   <- coef(x$model)[1]
    B_or   <- coef(x$model)[2]
    L50_or <- -A_or/B_or
    tab <- matrix(as.numeric(c(A_or, A_b, B_or, B_b, L50_or, L50_b)), 
                  nrow = 3, ncol = 2, byrow = TRUE)
    colnames(tab) <- c("Original", "Bootstrap (Median)")
    rownames(tab) <- c("A", "B", "L50")
    tab <- as.table(tab)
    return(tab)
  }
  
  return(invisible())
}


#' Plot method for morphMat (morphometric size at sexual maturity)
#'
#' @param x object of class 'morphMat' with the mature parameters and a data.frame with the allometric variables ("X", "Y") 
#' and classification of maturity. Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col color for the logistic curve and for the L50\% morphometric size at sexual maturity.
#' @param lwd line with for drawing fitted values and confidence intervals.
#' @param lty line type line type for drawing fitted values and confidence intervals
#' @param vline_hist color of the vertival lines in the histogram. The lines represent the 
#' the median and the confidece intervals.
#' @param lwd_hist line with for the vertical line in the histogram.
#' @param lty_hist line type for the vertical line in the histogram.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' \dontrun{
#' data(crabdata)
#' 
#' classify_data = classify_mature(crabdata, varNames = c("carapace_width", "chela_heigth"), 
#' varSex = "sex_category", selectSex = NULL, method = "ld")
#' 
#' my_mature = morph_mature(classify_data, method = "fq")
#' 
#' plot(my_mature, xlab = "Carapace width (mm.)", ylab = "Proportion mature", col = c("blue", "red"))
#' }
#' @export
#' @method plot morphMat
plot.morphMat <- function(x, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"), 
                          lwd = 2, lty = 2, vline_hist = "black", lwd_hist = 2, lty_hist = 2, ...){
  
  if (!inherits(x, "morphMat"))
    stop("Use only with 'morphMat' objects")
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  wide    <- quantile(x$L50_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  # figure 1
  hist(x$A_boot, main = "", xlab = "A", col = "grey90")
  abline(v = as.numeric(quantile(x$A_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = lwd_hist, col = vline_hist)
  abline(v = c(as.numeric(quantile(x$A_boot, probs = c(0.025, 0.975), na.rm = TRUE))), 
         lty = lty_hist, col = vline_hist)
  box()

  # figure 2
  hist(x$B_boot, main = "", xlab = "B", col = "grey90")
  abline(v = as.numeric(quantile(x$B_boot, probs = c(0.5), na.rm = TRUE)), 
         lwd = lwd_hist, col = vline_hist)
  abline(v = c(as.numeric(quantile(x$B_boot, probs = c(0.025, 0.975), na.rm = TRUE))), 
         lty = lty_hist, col = vline_hist)
  box()

  # figure 3
  hist(x$L50_boot, main = "", xlab = "Size at sexual maturity values", col = "grey90")
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
  points(wide[2], 0.5, pch = 19, col = col[2], cex = 1.25)
  legend("topleft", as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))), bty = "n")
  cat("Morphometric size at sexual maturity =", round(wide[2], 1), "\n")
  cat("Confidence intervals =", round(wide[1], 1), "-",round(wide[3], 1) ,  "\n")
  
  return(invisible(NULL))
}

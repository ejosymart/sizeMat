# ssmRG package: Size at Sexual Maturity based on Relative Growth---------------
#' @import MCMCpack
#' @import UsingR
#' @import biotools
#' @import matrixStats
#' @useDynLib ssmRG
#'
#' @title Size at Sexual Maturity based on Relative Growth
#'
#' @name ssmRG-package
#' @description  Package to Estimate size at sexual maturity from morphometic data based on relative growth.
#' @aliases ssmRG-package ssmRG
#' @docType package
#' @references ssmRG: Size at Sexual Maturity based on Relative Growth (RJournal)
#' @keywords size, sexual - maturity, alometric, relative-growth


NULL
#' Read data
#' 
#' Read a database with different extesions. The data base have only two variables, first variable is the independet variable and the second is the dependent variable.
#' @param file The filename
#' @param ext The file extension ("txt", "csv")
#' @return Database with x (independent) and y (dependent) variables.
#' @examples 
#' data = read_data("my_data", ext = "csv")
#' data
#' @export

read_data <- function(file, ext = "txt"){
  file <- paste(file, ".", ext, sep = "")
  
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
#' @param data The database with two variables (x: independent variable, y: dependent variable)
#' @return Database with x (independent), y (dependent) and mature classification
#' (juveniles = 0, adult = 1) variables.
#' @examples
#' data = classify_mature(data, methodPCA = "fq", ...)
#' data
#' plot(data, xlab = "X", ylab = "Y", colors = c(1, 2), pch = c(4, 5))
#' @export

classify_mature <- function(data){
  
  # Removing NA's rows ------------------------------------------------------
  data <- data[complete.cases(data), ]
  
  # Log data ----------------------------------------------------------------  
  data <- log(data)
  
  # Classify PCA ------------------------------------------------------------
  pca_classify    <- prcomp(data, 2)
  scores          <- pca_classify$x
  clusters        <- hclust(dist(scores, method = 'euclidean'), method = 'ward.D')
  mature_classify <- cutree(clusters, 2) - 1
  
  # Temporal base -----------------------------------------------------------
  base           <- data.frame(data, mature_binom = mature_classify)
  
  # linear o quadratic discriminarion analysis ------------------------------
  test_cov_mat <- boxM(base[, c("x", "y")], base[, "mature_binom"])
  if(test_cov_mat$p.value > 0.05){
    dis_reg    <- lda(mature_binom ~ ., data = base)
  }else{
    dis_reg    <- qda(mature_binom ~ ., data = base)
  }
  # Predict mature classification -------------------------------------------
  mature       <- as.numeric(as.character(predict(dis_reg)$class))
  # Print number in juveline and adult group --------------------------------
  cat("number in juveline and adult group =", as.numeric(table(mature)[1]), ",", as.numeric(table(mature)[2]), "\n")
  
  # Ouput -------------------------------------------------------------------
  data         <- data.frame(exp(base[, c("x", "y")]), mature)
  class(data)  <- c("classify")
  
  return(data)
}


plot.classify <- function(data, xlab = "X", ylab = "Y", colors = c(1, 2), pch = c(4, 5)){
  
  data <- data.frame(do.call("cbind", data))
  
  # Data filtering: juvenile and adults ------------------------------------
  juv <- data[data$mature == 0, ]
  adt <- data[data$mature == 1, ] 
  
  # Fit model to data ------------------------------------------------------  
  fit_juv <- glm(y ~ x, data = juv)
  fit_adt <- glm(y ~ x, data = adt)
  
  # Plot and details -------------------------------------------------------
  pch <- ifelse (data$mature == 0, pch[1], pch[2])
  col <- ifelse (data$mature == 0, colors[1], colors[2])
  plot(data[, c("x", "y")], col = col, xlab = xlab, ylab = ylab, pch = pch)
  lines(juv$x, predict(fit_juv), col = colors[1], lwd = 2)
  lines(adt$x, predict(fit_adt), col = colors[2], lwd = 2)
  eq_juv  <- paste0("Y = ", round(as.numeric(coef(fit_juv)[1]), 2), " + ", round(as.numeric(coef(fit_juv)[2]),2), " *X", sep = "")
  eq_adt  <- paste0("Y = ", round(as.numeric(coef(fit_adt)[1]), 2), " + ", round(as.numeric(coef(fit_adt)[2]),2), " *X", sep = "")
  legend("topleft", c(paste("Juveniles: ", eq_juv), paste("Adults: ", eq_adt)), 
         bty = "n", pch = pch, col = colors, cex = 0.8)
  invisible()
}


#' Calculate ogive
#' 
#' Estimate the size at 50\% maturity using a logistic regression
#' relating X variable and maturity stage (classifying individuals into juveniles 
#' or adults depending on their morphometry).
#' @param data The database with the X, Y and mature stages (juvelines = 0, adults = 1)
#' @param methodReg The method to be applied, "fq"frecuentist GLM, or "bayes" bayesian GLM
#' (MCMClogit function).
#' @return List database with the parameters and a data.frame with the X, Y and mature stages
#' variables. Also the fitted values for the logistic regression and confidence intervals.
#' @examples
#' my_ogive = calculate_ogive(data, methodReg = "fq")
#' plot(my_ogive, xlab = "X", ylab = "Proportion mature", col = "blue", col50 = "red")
#' @export

calculate_ogive <- function(data, methodReg = "fq"){
  
  data <- data.frame(do.call("cbind", data))
  
  # Estimate coefficients, Predict data and CI's ----------------------------  
  estimate <- switch(methodReg,
                     fq = .calculate_ogive_fq(data),
                     bayes = .calculate_ogive_bayes(data)) 
  
  # Output ------------------------------------------------------------------
  output  <- data.frame(x = data$x, y = data$y, mature = data$mature, CIlower = estimate$lower, 
                        fitted = estimate$fitted, CIupper = estimate$upper)
  data    <- list(params = estimate$params, output = output)
  class(data) <- "ogive"
  
  return(data)
}


.calculate_ogive_fq <- function(data){
  
  # GLM model ---------------------------------------------------------------  
  model_glm <- glm(mature ~ x, data = data, family = binomial(link = "logit"))
  
  # Coeficients -------------------------------------------------------------  
  A <- as.numeric(model_glm$coef[1])
  B <- as.numeric(model_glm$coef[2])
  params <- c(A, B)
  
  # Predict and Confidence intervals ----------------------------------------
  pred_dat <- predict(model_glm, newdata = data.frame(x = data$x), se.fit = TRUE)
  fitted   <- round(model_glm$fitted, 3)
  lower    <- round(with(pred_dat, exp(fit - 1.96 * se.fit) / (1+exp(fit - 1.96 * se.fit))), 3)
  upper    <- round(with(pred_dat, exp(fit + 1.96 * se.fit) / (1+exp(fit + 1.96 * se.fit))), 3)
  
  # Output ------------------------------------------------------------------  
  estimate <- list(params = params, lower = lower, fitted = fitted,  upper = upper)
  return(estimate)
}


.calculate_ogive_bayes <- function(data){
  
  # Bayesian linear regression model ----------------------------------------  
  model_bayes <- MCMClogit(data$mature ~ data$x, burnin = 1000, mcmc = 100000, thin = 10)
  
  # Coeficients ------------------------------------------------------------- 
  stats  <- summary(model_bayes)
  A      <- stats$quantiles[5]
  B      <- stats$quantiles[6]
  params <- c(A, B)
  
  # Predict and Credible intervals ------------------------------------------  
  create_x    <- cbind(1, data$x)
  x_bayes     <- as.matrix(create_x) %*% t(model_bayes)
  pred_bayes  <- as.data.frame(1 / (1 + exp(-x_bayes)))
  qtl         <- round(rowQuantiles(pred_bayes, probs = c(0.025, 0.5, 0.975)), 3)
  fitted      <- qtl[, 2]
  lower       <- qtl[, 1]
  upper       <- qtl[, 3]
  
  # Output ------------------------------------------------------------------  
  estimate    <- list(params = params, lower = lower, fitted = fitted,  upper = upper)
  return(estimate)
}


plot.ogive <- function(data, xlab = "X", ylab = "Proportion mature", col = "blue", col50 = "red"){
  
  fit     <- data$output
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  
  plot(sort(unique(x_input)), m_p, xlab = xlab, ylab = "Proportion mature", pch = 19, col = "darkgrey", axes = F)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  
  wide50 <- - data$params[1] / data$params[2]
  wide95 <- (1 / data$params[2]) * log(1 / 0.05 - 1) - data$params[1] / data$params[2]
  
  lines(sort(x_input), sort(fit$fitted), col = col, lwd = 2)
  lines(sort(x_input), sort(fit$CIlower), col = col, lwd = 2, lty = 2)
  lines(sort(x_input), sort(fit$CIupper), col = col, lwd = 2, lty = 2)
  lines(c(wide50, wide50), c(-1, 0.5), col = col50, lty = 2, lwd = 2)
  lines(c(-1, wide50), c(0.5, 0.5), col = col50, lty = 2, lwd = 2)
  legend("topleft", as.expression(bquote(bold(CW[50] == .(round(wide50, 1))))), bty = "n")
  cat("Carapace width of 50% maturity =", round(wide50, 1), "\n")
  cat("Carapace width of 95% maturity =", round(wide95, 1), "\n")
  invisible()
}

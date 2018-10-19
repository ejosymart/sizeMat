#' Estimate gonadal maturity
#' 
#' Estimate size at gonad maturity.
#' @param data data.frame with allometric variables and stage of sexual maturity (gonad maturation stages).
#' @param varNames a character string indicating the name of the allometric 
#' and the stage of sexual maturity variables to be used for analysis.
#' @param inmName a character string indicating the name or names of the inmaturity stage.
#' @param matName a character string indicating the name or names of the maturity stage.
#' @param method a character string indicating the method to be applied, \code{"fq"} frecuentist GLM, or \code{"bayes"} Bayes GLM (MCMClogit function). 
#' @param niter number of iterations (bootstrap resampling).
#' @param seed a single value, interpreted as an integer.
#' @return An object of class 'gonadMat'.
#' 
#' \code{model} the summary statistics of the model.
#' 
#' \code{A_boot} the 'n iter' values of parameter A.
#' 
#' \code{B_boot} the 'n iter' values of parameter B.
#' 
#' \code{L50} the 'n iter' values of parameter L50 (size at gonad maturity).
#' 
#' \code{out} a dataframe with the allometric variable "X", stage of sexual maturity, the fitted values for  
#' logistic regression and confidence intervals (95\%). Also the summary statistics of the model is provided.
#' @details Estimate the size at gonad maturity using a logistic regression with X variable and 
#' stages of sexual maturity (two categories: inmature and mature). 
#' 
#' The function requires a data.frame with the X (allometric variable) and 
#' the stage of sexual maturity (gonad maturation stage).
#' 
#' The argument \code{varNames} requires a character string indicating the name of one allometric and the stage
#' of sexual maturity variable to be used for analysis (e.g \code{varNames = c("total_length", "stage_mat")}). 
#' So the argument \code{varNames} must contain two character strings only, the first is the allometric variable 
#' and the second is the stage of sexual maturity.
#' 
#' The arguments \code{inmName} and \code{matName} require a character string indicanting the name 
#' of the stages of sexual maturity in the data.frame. The argument could contain one character string 
#' or could be a vector (e.g \code{inmName = "I"}, \code{matName = c("II", "III", "IV")}).
#' 
#' The argument \code{method} requires a character string indicating which regression will be used for the test.
#' If \code{method = "fq"} the logistic regression is based on GLM (frequentist), if \code{method = "bayes"} a sample from 
#' the posterior distribution of a logistic regression model using a random walk Metropolis algorithm is generated (see MCMClogit function).
#' 
#' The argument \code{niter} requires a number. For the GLM regression (\code{method = "fq"}), a non-parametric bootstrap method consists
#' in generate B bootstrap samples, by resampling with replacement the original data. Then all statistics for each parameter 
#' can be calculated from each bootstrap sample (median and confidence intervals). 
#' For the \code{method = "bayes"}, the argument \code{niter} is related to the number of Metropolis iterations for the sampler.
#' @exportClass gonad_mature
#' @examples
#' data(matFish)
#' 
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), inmName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' # 'niter' parameters:
#' gonad_mat$A_boot
#' gonad_mat$B_boot
#' gonad_mat$L50_boot
#' gonad_mat$out
#' @export
gonad_mature <- function(data, varNames = c("allometric", "stage") , inmName = "inm", matName = "mad", 
                         method = "fq", niter = 999, seed = 070388){
  
  if(length(varNames) != 2) stop("You must provide two variables.")
  if(!all(varNames %in% names(data))) stop("'varNames' have not been found in data.")
  
  data <- data[, varNames]
  names(data) <- c("x", "stage")
  data$stage <- as.factor(data$stage)
  data <- data[complete.cases(data), ] 
  
  if(!all(c(inmName, matName) %in% levels(data$stage))) stop("'inmName' or 'matName' have not been found in data.")
  if(all(inmName %in% matName)) stop("'inmName' and 'matName' must have different stage names")
  
  data$stage <- ifelse(data$stage %in% inmName, 0, 1)
  estimate <- switch(method,
                     fq = .gonad_mature_fq(data = data, niter = niter, seed = seed),
                     bayes = .gonad_mature_bayes(data = data, niter = niter, seed = seed))
  
  out   <- data.frame(x = data$x, mature = data$stage, CIlower = estimate$lower, 
                      fitted = estimate$fitted, CIupper = estimate$upper)
  
  output <- list(model    = estimate$model,
                 A_boot   = estimate$parameters_A, 
                 B_boot   = estimate$parameters_B,
                 L50_boot = estimate$L50,
                 out      = out)
  
  class(output) <- c("gonadMat", class(output))
  
  return(output)
}


#' Print method for gonadMat class (size at gonad maturity)
#'
#' @param x object of class 'gonadMat' with the parameters of the logistic regression and a data.frame with the X and stage of sexual maturity.
#' variables. Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @param \dots Additional arguments to the print method.
#' @return The median of the size at gonad maturity estimation, parameters and the Nagelkerke's R squared.
#' @examples
#' data(matFish)
#' 
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), inmName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' print(gonad_mat)
#' @export
#' @method print gonadMat
print.gonadMat <- function(x, ...){
  if (!inherits(x, "gonadMat"))
    stop("Use only with 'gonadMat' objects")
  cat("formula: Y = 1/1+exp-(A + B*X)", "\n\n")
  
  A_b   <- quantile(x$A_boot, probs = 0.5, na.rm = TRUE)
  B_b   <- quantile(x$B_boot, probs = 0.5, na.rm = TRUE)
  L50_b <- quantile(x$L50_boot, probs = 0.5, na.rm = TRUE)
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature

  model1 <- glm(y_input ~ x_input, family = binomial(link = "logit"))
  R2     <- nagelkerkeR2(model1)
  
    if(is.null(coef(x$model))){
    tab <- matrix(c(round(as.numeric(c(A_b, B_b, L50_b)), 4), round(R2, 4)), 
                  nrow = 4, ncol = 1, byrow = TRUE)
    colnames(tab) <- c("Bootstrap (Median)")
    rownames(tab) <- c("A", "B", "L50", "R2")
    tab <- as.table(tab)
    return(tab)
  }else{
    A_or   <- coef(x$model)[1]
    B_or   <- coef(x$model)[2]
    L50_or <- -A_or/B_or
    tab <- matrix(c(round(as.numeric(c(A_or, A_b, B_or, B_b, L50_or, L50_b)), 4), "-", round(R2, 4)), 
                  nrow = 4, ncol = 2, byrow = TRUE)
    colnames(tab) <- c("Original", "Bootstrap (Median)")
    rownames(tab) <- c("A", "B", "L50", "R2")
    tab <- as.table(tab)
    return(tab)
  }
  
  return(invisible())
}



#' Plot method for gonadMat class (size at gonad maturity)
#'
#' @param x object of class 'gonadMat' with the mature parameters and a data.frame with the X and stage of sexual maturity.
#' variables. Also the fitted values for the logistic regression and confidence intervals (95\%).
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col color for the logistic curve and for the L50\% size at gonad maturity.
#' @param lwd line with for drawing fitted values and confidence intervals.
#' @param lty line type line type for drawing fitted values and confidence intervals
#' @param vline_hist color of the vertival lines in the histogram. The lines represent the 
#' the median and the confidece intervals.
#' @param lwd_hist line with for the vertical line in the histogram.
#' @param lty_hist line type for the vertical line in the histogram.
#' @param onlyOgive plot only the ogive.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' data(matFish)
#' 
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), inmName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' plot(gonad_mat, xlab = "Total length (cm.)", ylab = "Proportion mature", col = c("blue", "red"))
#' @export
#' @method plot gonadMat
plot.gonadMat <- function(x, xlab = "X", ylab = "Proportion mature", col = c("blue", "red"), 
                        lwd = 2, lty = 2, vline_hist = "black", lwd_hist = 2, lty_hist = 2, 
                        onlyOgive = FALSE, ...){
  
  if (!inherits(x, "gonadMat"))
    stop("Use only with 'gonadMat' objects")
  
  fit     <- x$out
  x_input <- fit$x
  y_input <- fit$mature
  m_p     <- tapply(y_input, x_input, mean)
  wide    <- quantile(x$L50_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  # R-square Nagelkerke method
  model1 <- glm(y_input ~ x_input, family = binomial(link = "logit"))
  R2     <- nagelkerkeR2(model1)
  
  if(onlyOgive == TRUE){
    if(length(col) < 2) stop('col argument must have 2 values. The colors could be the same')
    if(length(col) > 2) warning('col: only the first two colors will be used in the plot')
    plot(sort(unique(x_input)), m_p, xlab = xlab, ylab = ylab, pch = 19, col = "darkgrey", ...)
    lines(sort(x_input), sort(fit$fitted), col = col[1], lwd = lwd)
    lines(sort(x_input), sort(fit$CIlower), col = col[1], lwd = lwd, lty = lty)
    lines(sort(x_input), sort(fit$CIupper), col = col[1], lwd = lwd, lty = lty)
    lines(c(wide[2], wide[2]), c(-1, 0.5), col = col[2], lwd = lwd, lty = lty)
    lines(c(-1, wide[2]), c(0.5, 0.5), col = col[2], lwd = lwd, lty = lty)
    points(wide[2], 0.5, pch = 19, col = col[2], cex = 1.25)
  }else{
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
    legend("topleft", c(as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))),
                        as.expression(bquote(bold(R^2 == .(round(R2, 2)))))), 
           bty = "n")
  }
  
  cat("Size at gonad maturity =", round(wide[2], 1), "\n")
  cat("Confidence intervals =", round(wide[1], 1), "-",round(wide[3], 1) ,  "\n")
  cat("Rsquare =", round(R2, 2))
  
  return(invisible(NULL))
}

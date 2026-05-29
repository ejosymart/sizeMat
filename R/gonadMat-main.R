#' Estimate gonadal maturity
#' 
#' Estimate size at gonad maturity.
#' @param data data.frame with allometric variables and stage of sexual maturity (gonad maturation stages).
#' @param varNames a character string indicating the name of the allometric 
#' and the stage of sexual maturity variables to be used for analysis.
#' @param immName a character string indicating the name or names of the immaturity stage.
#' @param matName a character string indicating the name or names of the maturity stage.
#' @param method a character string indicating the method to be applied, \code{"fq"} frequentist GLM, or \code{"bayes"} Bayes GLM (MCMClogit function). 
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
#' stages of sexual maturity (two categories: immature and mature). 
#' 
#' The function requires a data.frame with the X (allometric variable) and 
#' the stage of sexual maturity (gonad maturation stage).
#' 
#' The argument \code{varNames} requires a character string indicating the name of one allometric and the stage
#' of sexual maturity variable to be used for analysis (i.e. \code{varNames = c("total_length", "stage_mat")}). 
#' So the argument \code{varNames} must contain two character strings only, the first is the allometric variable 
#' and the second is the stage of sexual maturity.
#' 
#' The arguments \code{immName} and \code{matName} require a character string indicating the name 
#' of the stages of sexual maturity in the data.frame. The argument could contain one character string 
#' or should be a vector (i.e. \code{immName = "I"}, \code{matName = c("II", "III", "IV")}).
#' 
#' The argument \code{method} requires a character string indicating which regression will be used for the test.
#' If \code{method = "fq"} the logistic regression is based on GLM (frequentist), if \code{method = "bayes"} a sample from 
#' the posterior distribution of a logistic regression model using a random walk Metropolis algorithm is generated (see MCMClogit function).
#' 
#' The argument \code{niter} requires a number. For the GLM regression (\code{method = "fq"}), a non-parametric bootstrap method consists
#' in generate B bootstrap samples, by resampling with replacement the original data. Then all statistics for each parameter 
#' can be calculated from each bootstrap sample (median and confidence intervals). 
#' For the \code{method = "bayes"}, the argument \code{niter} is related to the number of Metropolis iterations for the sampler.
#' @examples
#' data(matFish)
#' 
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), immName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' # 'niter' parameters:
#' gonad_mat$A_boot
#' gonad_mat$B_boot
#' gonad_mat$L50_boot
#' gonad_mat$out
#' @export
gonad_mature <- function(data, varNames = c("allometric", "stage") , immName = "imm", matName = "mad", 
                         method = "fq", niter = 999, seed = 070388){
  
  if(length(varNames) != 2) stop("You must provide two variables.")
  if(!all(varNames %in% names(data))) stop("'varNames' have not been found in data.")
  
  data <- data[, varNames]
  names(data) <- c("x", "stage")
  data$stage <- as.factor(data$stage)
  data <- data[complete.cases(data), ]
  
  if(!all(c(immName, matName) %in% levels(data$stage))) stop("'immName' or 'matName' have not been found in data.")
  if(all(immName %in% matName)) stop("'immName' and 'matName' must have different stage names")
  
  data$stage <- ifelse(data$stage %in% immName, 0, 1)
  estimate <- switch(method,
                     fq = .gonad_mature_fq(data = data, niter = niter, seed = seed),
                     bayes = .gonad_mature_bayes(data = data, niter = niter, seed = seed))
  
  out   <- data.frame(x = data$x, mature = data$stage, CIlower = estimate$lower, 
                      fitted = estimate$fitted, CIupper = estimate$upper)
  
  output <- list(model    = estimate$model,
                 modelglm = estimate$modelglm,
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
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), immName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' print(gonad_mat)
#' @export
#' @method print gonadMat
print.gonadMat <- function(x, ...){
  if (!inherits(x, "gonadMat"))
    stop("Use only with 'gonadMat' objects")
  cat("formula: Y = 1/[1+exp-(A + B*X)]", "\n\n")
  
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
    tab <- matrix(c(round(as.numeric(c(A_or, A_b, B_or, B_b, L50_or, L50_b)), 4), round(R2, 4), "-" ), 
                  nrow = 4, ncol = 2, byrow = TRUE)
    colnames(tab) <- c("Original", "Bootstrap (Median)")
    rownames(tab) <- c("A", "B", "L50", "R2")
    tab <- as.table(tab)
    return(tab)
  }
  
  return(invisible())
}



#' Plot method for gonadMat class (size at gonadal maturity)
#'
#' @param x object of class 'gonadMat' with mature parameters and fitted values
#'   from the logistic regression.
#' @param xlab title for the x axis.
#' @param ylab title for the y axis.
#' @param col colors for the logistic curve and the L50 reference lines.
#' @param lwd line width for fitted values and confidence intervals.
#' @param lty line type for fitted values and confidence intervals.
#' @param vline_hist color of the vertical lines in the histograms.
#' @param lwd_hist line width for vertical lines in the histograms.
#' @param lty_hist line type for vertical lines in the histograms.
#' @param onlyOgive logical. If TRUE, plot only the maturity ogive.
#' @param showLegend logical. If TRUE, show L50 and R2 labels.
#' @param legendPosition position of the legend in base graphics and ggplot2 style. Options are "topleft", "topright", "bottomleft", and "bottomright".
#' @param gg_style logical. If TRUE, return ggplot2-style plots.
#' @param pch plotting character for observed maturity proportions.
#' @param cex point size in base graphics and ggplot2 style.
#' @param point_alpha transparency level for points, used only when gg_style = TRUE.
#' @param base_size base font size, used only when gg_style = TRUE.
#' @param label_size size of L50 and R2 labels, used only when gg_style = TRUE.
#' @param \dots additional arguments passed to base graphics.
#' @examples
#' data(matFish)
#' 
#' gonad_mat = gonad_mature(matFish, varNames = c("total_length", "stage_mat"), immName = "I", 
#' matName = c("II", "III", "IV"), method = "fq", niter = 50)
#' 
#' plot(gonad_mat, xlab = "Total length (cm.)", ylab = "Proportion mature", col = c("blue", "red"))
#' @export
#' @method plot gonadMat
plot.gonadMat <- function(x, xlab = "X", ylab = "Proportion mature",
                          col = c("blue", "red"),
                          lwd = 2, lty = 2,
                          vline_hist = "black",
                          lwd_hist = 2, lty_hist = 2,
                          onlyOgive = FALSE, showLegend = TRUE,
                          legendPosition = "topleft",
                          gg_style = FALSE,
                          pch = 19, cex = 2.5,
                          point_alpha = 0.8, base_size = 13,
                          label_size = 5, ...) {
  
  if (!inherits(x, "gonadMat")) {
    stop("Use only with 'gonadMat' objects")
  }
  
  if (length(col) < 2) {
    stop("'col' argument must have 2 values. The colors could be the same.")
  }
  
  if (length(col) > 2) {
    warning("'col': only the first two colors will be used in the plot.")
  }
  
  col <- col[1:2]
  
  if (!legendPosition %in% c("topleft", "topright", "bottomleft", "bottomright")) {
    stop("'legendPosition' must be one of: 'topleft', 'topright', 'bottomleft', 'bottomright'.")
  }
  
  fit <- x$out
  
  if (!all(c("x", "mature", "fitted", "CIlower", "CIupper") %in% names(fit))) {
    stop("The 'out' object must contain 'x', 'mature', 'fitted', 'CIlower', and 'CIupper'.")
  }
  
  x_input <- fit$x
  y_input <- fit$mature
  
  m_p <- tapply(y_input, x_input, mean)
  
  wide <- stats::quantile(x$L50_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  model1 <- stats::glm(y_input ~ x_input, family = stats::binomial(link = "logit"))
  
  R2 <- nagelkerkeR2(model1)
  
  ogive_points <- data.frame(x = as.numeric(names(m_p)), proportion = as.numeric(m_p))
  
  ogive_points <- ogive_points[order(ogive_points$x), ]
  
  ogive_fit <- data.frame(x = fit$x, fitted = fit$fitted, 
                          CIlower = fit$CIlower, CIupper = fit$CIupper)
  
  ogive_fit <- ogive_fit[order(ogive_fit$x), ]
  
  L50_text <- paste0("L[50] == ", round(wide[2], 1))
  R2_text  <- paste0("R^2 == ", round(R2, 2))
  
  print_summary <- function(){
    cat("Size at gonadal maturity =", round(wide[2], 1), "\n")
    cat("Confidence intervals =", round(wide[1], 1), "-", round(wide[3], 1), "\n")
    cat("Rsquare =", round(R2, 2), "\n")
  }
  
  if(isTRUE(gg_style)){
    
    if(!requireNamespace("ggplot2", quietly = TRUE)){
      stop("Package 'ggplot2' is required when gg_style = TRUE.")
    }
    
    label_x <- switch(legendPosition, 
                      "topleft" = min(ogive_fit$x, na.rm = TRUE), 
                      "topright" = max(ogive_fit$x, na.rm = TRUE), 
                      "bottomleft" = min(ogive_fit$x, na.rm = TRUE), 
                      "bottomright" = max(ogive_fit$x, na.rm = TRUE))
    
    label_y1 <- switch(legendPosition, "topleft" = 0.96, 
                       "topright" = 0.96, "bottomleft" = 0.16, 
                       "bottomright" = 0.16)
    
    label_y2 <- switch(legendPosition, "topleft" = 0.86, 
                       "topright" = 0.86, "bottomleft" = 0.06, 
                       "bottomright" = 0.06)
    
    label_hjust <- switch(legendPosition, "topleft" = 0, 
                          "topright" = 1, "bottomleft" = 0, 
                          "bottomright" = 1)
    
    p_ogive <- ggplot2::ggplot() + 
      ggplot2::geom_point(data = ogive_points, 
                          ggplot2::aes(x = x, y = proportion), 
                          colour = "darkgrey", shape = pch, 
                          size = cex, alpha = point_alpha) +
      ggplot2::geom_line(data = ogive_fit, 
                         ggplot2::aes(x = x, y = fitted), 
                         colour = col[1], 
                         linewidth = lwd) +
      ggplot2::geom_line(data = ogive_fit, 
                         ggplot2::aes(x = x, y = CIlower),
                         colour = col[1], 
                         linewidth = lwd, linetype = lty) +
      ggplot2::geom_line(data = ogive_fit, 
                         ggplot2::aes(x = x, y = CIupper), 
                         colour = col[1], 
                         linewidth = lwd, linetype = lty) +
      ggplot2::geom_segment(ggplot2::aes(x = wide[2], xend = wide[2], 
                                         y = 0, yend = 0.5), 
                            colour = col[2], 
                            linewidth = lwd, linetype = lty) +
      ggplot2::geom_segment(ggplot2::aes(x = min(ogive_fit$x, na.rm = TRUE), 
                                         xend = wide[2], y = 0.5, yend = 0.5), 
                            colour = col[2], 
                            linewidth = lwd, linetype = lty) +
      ggplot2::geom_point(ggplot2::aes(x = wide[2], y = 0.5), 
                          colour = col[2], shape = pch, size = cex * 0.8) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"), 
                     axis.text = ggplot2::element_text(colour = "black"), 
                     plot.margin = ggplot2::margin(8, 8, 8, 8))
    
    if(isTRUE(showLegend)){
      p_ogive <- p_ogive +
        ggplot2::annotate("text", x = label_x, y = label_y1, 
                          label = L50_text, parse = TRUE, 
                          hjust = label_hjust, fontface = "bold", 
                          size = label_size) +
        ggplot2::annotate("text", x = label_x, y = label_y2, 
                          label = R2_text, parse = TRUE, 
                          hjust = label_hjust, 
                          fontface = "bold", 
                          size = label_size)
    }
    
    if(isTRUE(onlyOgive)){
      print_summary()
      return(p_ogive)
    }
    
    make_hist_plot <- function(values, label_x){
      
      values <- values[is.finite(values)]
      
      q <- stats::quantile(values, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      
      ggplot2::ggplot(data.frame(value = values), ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey90", colour = "grey40") +
        ggplot2::geom_vline(xintercept = q[2], colour = vline_hist, linewidth = lwd_hist) +
        ggplot2::geom_vline(xintercept = c(q[1], q[3]), colour = vline_hist, 
                            linewidth = lwd_hist, linetype = lty_hist) +
        ggplot2::labs(x = label_x, y = "Frequency") +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"), 
                       axis.text = ggplot2::element_text(colour = "black"), 
                       plot.margin = ggplot2::margin(8, 8, 8, 8))
    }
    
    p_A <- make_hist_plot(x$A_boot, "A")
    p_B <- make_hist_plot(x$B_boot, "B")
    p_L50 <- make_hist_plot(x$L50_boot, "Size at sexual maturity values")
    
    out <- list(A = p_A, B = p_B, L50 = p_L50, ogive = p_ogive)
    
    print_summary()
    
    return(out)
  }
  
  if(isTRUE(onlyOgive)){
    
    graphics::plot(ogive_points$x, ogive_points$proportion, 
                   xlab = xlab, ylab = ylab, 
                   pch = pch,col = "darkgrey", ...)
    
    graphics::lines(ogive_fit$x, ogive_fit$fitted, col = col[1], lwd = lwd)
    
    graphics::lines(ogive_fit$x, ogive_fit$CIlower, col = col[1], 
                    lwd = lwd, lty = lty)
    
    graphics::lines(ogive_fit$x, ogive_fit$CIupper, col = col[1], 
                    lwd = lwd, lty = lty)
    
    graphics::lines(c(wide[2], wide[2]), c(0, 0.5), 
                    col = col[2], lwd = lwd, lty = lty)
    
    graphics::lines(c(min(ogive_fit$x, na.rm = TRUE), wide[2]), 
                    c(0.5, 0.5), col = col[2], lwd = lwd, lty = lty)
    
    graphics::points(wide[2], 0.5, pch = pch, col = col[2], cex = 1.25)
    
    if (isTRUE(showLegend)){
      graphics::legend(legendPosition, c(as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))),
                                         as.expression(bquote(bold(R^2 == .(round(R2, 2)))))), 
                       bty = "n")
    }
    
    print_summary()
    
    return(invisible(NULL))
  }
  
  graphics::hist(x$A_boot, main = "", xlab = "A", col = "grey90")
  
  graphics::abline(v = as.numeric(stats::quantile(x$A_boot, probs = 0.5, na.rm = TRUE)), 
                   lwd = lwd_hist, col = vline_hist)
  
  graphics::abline(v = as.numeric(stats::quantile(x$A_boot, probs = c(0.025, 0.975), na.rm = TRUE)), 
                   lty = lty_hist, col = vline_hist)
  
  graphics::box()
  
  graphics::hist(x$B_boot, main = "", xlab = "B", col = "grey90")
  
  graphics::abline(v = as.numeric(stats::quantile(x$B_boot, probs = 0.5, na.rm = TRUE)), 
                   lwd = lwd_hist, col = vline_hist)
  
  graphics::abline(v = as.numeric(stats::quantile(x$B_boot, probs = c(0.025, 0.975), na.rm = TRUE)), 
                   lty = lty_hist, col = vline_hist)
  
  graphics::box()
  
  graphics::hist(x$L50_boot, main = "", 
                 xlab = "Size at sexual maturity values", 
                 col = "grey90")
  
  graphics::abline(v = as.numeric(stats::quantile(x$L50_boot, probs = 0.5, na.rm = TRUE)), 
                   lwd = lwd_hist, col = vline_hist)
  
  graphics::abline(v = as.numeric(stats::quantile(x$L50_boot, probs = c(0.025, 0.975), na.rm = TRUE)),
                   lty = lty_hist, col = vline_hist)
  
  graphics::box()
  
  graphics::plot(ogive_points$x, ogive_points$proportion, 
                 xlab = xlab, ylab = ylab, pch = pch, col = "darkgrey", ...)
  
  graphics::lines(ogive_fit$x, ogive_fit$fitted, col = col[1], 
                  lwd = lwd)
  
  graphics::lines(ogive_fit$x, ogive_fit$CIlower, col = col[1], 
                  lwd = lwd, lty = lty)
  
  graphics::lines(ogive_fit$x, ogive_fit$CIupper, col = col[1], 
                  lwd = lwd, lty = lty)
  
  graphics::lines(c(wide[2], wide[2]), c(0, 0.5), col = col[2], 
                  lwd = lwd, lty = lty)
  
  graphics::lines(c(min(ogive_fit$x, na.rm = TRUE), wide[2]), c(0.5, 0.5), 
                  col = col[2], lwd = lwd, lty = lty)
  
  graphics::points(wide[2], 0.5, pch = pch, col = col[2], cex = 1.25)
  
  if (isTRUE(showLegend)){
    graphics::legend(legendPosition,
                     c(as.expression(bquote(bold(L[50] == .(round(wide[2], 1))))), 
                       as.expression(bquote(bold(R^2 == .(round(R2, 2)))))), 
                     bty = "n")
  }
  
  print_summary()
  
  invisible(NULL)
}

#' Nagelkerme method R-square
#' 
#' Estimate Nagelkerke's R squared from the result of glm(). Evaluate the goodness of fit  for logistic regression. 
#' 
#' @param x An object of class 'glm'.
#'
#' @return Rsquare Nagelkerke's R squared.
#' @examples
#' set.seed(7388)
#' n <- 300
#' x <- rnorm(n)
#' a <- 1
#' b <- -2
#' p <- 1/(1+exp(a+b*x))
#' y <- factor(ifelse(runif(n) < p, 1, 0), levels = 0:1)
#' mod1 <- glm(y ~ x, family=binomial)
#' nagelkerkeR2(mod1)
#' @export
nagelkerkeR2 <- function(x){
  if (!inherits(x, "glm") | family(x)$family != "binomial")
    stop("Use with 'glm' class and family 'binomial' only")
  x2 <- update(x, ~ 1)
  M1 <- x$null.deviance - x$deviance
  M2 <- x2$null.deviance - x2$deviance
  M3 <- abs(M1 - M2)
  r2 <- as.numeric((1 - exp(-M3/nrow(x$model)))/(1 - exp(2*as.numeric(logLik(x2)/nrow(x$model)))))  
  return(r2)
}

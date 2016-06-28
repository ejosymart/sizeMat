.calculate_ogive_fq <- function(input){
  
  model_glm <- glm(mature ~ x, data = input, family = binomial(link = "logit"))
  
  A <- as.numeric(model_glm$coef[1])
  B <- as.numeric(model_glm$coef[2])
  params <- c(A, B)
  
  pred_dat <- predict(model_glm, newdata = data.frame(x = input$x), se.fit = TRUE)
  fitted   <- round(model_glm$fitted, 3)
  lower    <- round(with(pred_dat, exp(fit - 1.96 * se.fit) / (1+exp(fit - 1.96 * se.fit))), 3)
  upper    <- round(with(pred_dat, exp(fit + 1.96 * se.fit) / (1+exp(fit + 1.96 * se.fit))), 3)
  
  estimate <- list(params = params, lower = lower, fitted = fitted,  upper = upper)
  return(estimate)
}


.calculate_ogive_bayes <- function(input){
  
  model_bayes <- MCMClogit(input$mature ~ input$x, burnin = 1000, mcmc = 100000, thin = 10)
  
  stats  <- summary(model_bayes)
  A      <- stats$quantiles[5]
  B      <- stats$quantiles[6]
  params <- c(A, B)
  
  create_x    <- cbind(1, input$x)
  x_bayes     <- as.matrix(create_x) %*% t(model_bayes)
  pred_bayes  <- as.data.frame(1 / (1 + exp(-x_bayes)))
  qtl         <- round(rowQuantiles(pred_bayes, probs = c(0.025, 0.5, 0.975)), 3)
  fitted      <- qtl[, 2]
  lower       <- qtl[, 1]
  upper       <- qtl[, 3]
  
  estimate    <- list(params = params, lower = lower, fitted = fitted,  upper = upper)
  return(estimate)
}

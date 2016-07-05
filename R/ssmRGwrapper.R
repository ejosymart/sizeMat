.calculate_ogive_fq <- function(input, niter = niter, seed = seed){
  
  set.seed(seed)
  new_input <- list()
  n_coef <- list()
  for(i in 1:niter){
  new_input[[i]] <- input[sample(nrow(input), nrow(input), replace = T), ]
  model_glm <- glm(mature ~ x, data = new_input[[i]], family = binomial(link = "logit"))
  glm_coef  <- coef(model_glm)
  n_coef    <- rbind(glm_coef, n_coef)
  }
  
  A    <- as.numeric(n_coef[,1])
  B    <- as.numeric(n_coef[,2])
  L50  <- -A/B

  create_x <- cbind(1, input$x)
  x_fq     <- as.matrix(create_x) %*% t(as.matrix(cbind(A,B)))
  pred_fq  <- as.data.frame(1 / (1 + exp(-x_fq)))
  qtl      <- round(rowQuantiles(pred_fq, probs = c(0.05, 0.5, 0.95)), 3)
  fitted   <- qtl[, 2]
  lower    <- qtl[, 1]
  upper    <- qtl[, 3]
  
  estimate <- list(parameters_A = A, parameters_B = B, L50 = L50,
                   lower = lower, fitted = fitted,  upper = upper)
  
  return(estimate)
}


.calculate_ogive_bayes <- function(input, niter = niter, seed = seed){
  
  set.seed(seed)
  model_bayes <- MCMClogit(input$mature ~ input$x, mcmc = niter, thin = 1)
  
  A    <- as.numeric(model_bayes[,1])
  B    <- as.numeric(model_bayes[,2])
  L50  <- -A/B
  
  create_x    <- cbind(1, input$x)
  x_bayes     <- as.matrix(create_x) %*% t(model_bayes)
  pred_bayes  <- as.data.frame(1 / (1 + exp(-x_bayes)))
  qtl         <- round(rowQuantiles(pred_bayes, probs = c(0.05, 0.5, 0.95)), 3)
  fitted      <- qtl[, 2]
  lower       <- qtl[, 1]
  upper       <- qtl[, 3]
  
  estimate    <- list(parameters_A = A, parameters_B = B, L50 = L50, 
                      lower = lower, fitted = fitted,  upper = upper)
  
  return(estimate)
}

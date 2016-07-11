.calculate_ogive_fq <- function(data, niter = niter, seed = seed){
  
  model_glm <- glm(mature ~ x, data = data, family = binomial(link = "logit"))
  smry_model <- summary(model_glm)
  
  set.seed(seed)
  new_data <- list()
  n_coef <- list()
  for(i in 1:niter){
  new_data[[i]] <- data[sample(nrow(data), nrow(data), replace = T), ]
  model_boot <- glm(mature ~ x, data = new_data[[i]], family = binomial(link = "logit"))
  glm_coef  <- coef(model_boot)
  n_coef    <- rbind(glm_coef, n_coef)
  }
  
  A    <- as.numeric(n_coef[,1])
  B    <- as.numeric(n_coef[,2])
  L50  <- -A/B

  create_x <- cbind(1, data$x)
  x_fq     <- as.matrix(create_x) %*% t(as.matrix(cbind(A,B)))
  pred_fq  <- as.data.frame(1 / (1 + exp(-x_fq)))
  qtl      <- round(matrixStats::rowQuantiles(pred_fq, probs = c(0.025, 0.5, 0.975)), 3)
  fitted   <- qtl[, 2]
  lower    <- qtl[, 1]
  upper    <- qtl[, 3]
  
  estimate <- list(model = smry_model,
                   parameters_A = A, 
                   parameters_B = B, 
                   L50 = L50,
                   lower = lower, 
                   fitted = fitted,  
                   upper = upper)
  
  return(estimate)
}


.calculate_ogive_bayes <- function(data, niter = niter, seed = seed){
  
  set.seed(seed)
  model_bayes <- MCMCpack::MCMClogit(data$mature ~ data$x, mcmc = niter, thin = 1)
  smry_model <- summary(model_bayes)
  A    <- as.numeric(model_bayes[,1])
  B    <- as.numeric(model_bayes[,2])
  L50  <- -A/B
  
  create_x    <- cbind(1, data$x)
  x_bayes     <- as.matrix(create_x) %*% t(model_bayes)
  pred_bayes  <- as.data.frame(1 / (1 + exp(-x_bayes)))
  qtl         <- round(matrixStats::rowQuantiles(pred_bayes, probs = c(0.025, 0.5, 0.975)), 3)
  fitted      <- qtl[, 2]
  lower       <- qtl[, 1]
  upper       <- qtl[, 3]
  
  estimate    <- list(model = smry_model, 
                      parameters_A = A, 
                      parameters_B = B, 
                      L50 = L50, 
                      lower = lower, 
                      fitted = fitted, 
                      upper = upper)
  
  return(estimate)
}

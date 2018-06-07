.morph_mature_fq <- function(data, niter, seed){
  
  model_glm <- glm(mature ~ x, data = data, family = binomial(link = "logit"))
  smry_model <- summary(model_glm)
  
  set.seed(seed)
  n_coef   <- list()
  for(i in seq_len(niter)){
    new_data   <- data[sample(nrow(data), replace = TRUE), ]
    model_boot <- glm(mature ~ x, data = new_data, family = binomial(link = "logit"))
    glm_coef   <- coef(model_boot)
    n_coef     <- rbind(n_coef, glm_coef)
  }
  
  A    <- as.numeric(n_coef[,1])
  B    <- as.numeric(n_coef[,2])
  L50  <- -A/B

  create_x <- cbind(1, data$x)
  x_fq     <- as.matrix(create_x) %*% t(as.matrix(cbind(A,B)))
  pred_fq  <- 1 / (1 + exp(-x_fq))
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


.morph_mature_bayes <- function(data, niter, seed){
  
  set.seed(seed)
  model_bayes <- MCMCpack::MCMClogit(data$mature ~ data$x, mcmc = niter, thin = 1)
  smry_model  <- summary(model_bayes)
  A    <- as.numeric(model_bayes[,1])
  B    <- as.numeric(model_bayes[,2])
  L50  <- -A/B
  
  create_x   <- cbind(1, data$x)
  x_bayes    <- as.matrix(create_x) %*% t(model_bayes)
  pred_bayes <- 1 / (1 + exp(-x_bayes))
  qtl        <- round(matrixStats::rowQuantiles(pred_bayes, probs = c(0.025, 0.5, 0.975)), 3)
  fitted     <- qtl[, 2]
  lower      <- qtl[, 1]
  upper      <- qtl[, 3]
  
  estimate   <- list(model = smry_model, 
                     parameters_A = A, 
                     parameters_B = B, 
                     L50 = L50, 
                     lower = lower, 
                     fitted = fitted, 
                     upper = upper)
  
  return(estimate)
}

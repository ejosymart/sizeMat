# Load data ---------------------------------------------------------------

get_data <- function(file, ext = "txt", ...){
  file <- paste(file, ".", ext, sep = "")
  
  data <- switch(ext,
                "txt"   = read.table(file, header = TRUE, stringsAsFactors = FALSE),
                "csv"   = read.csv(file),
                "csv2"  = read.csv2(file),
                "delim" = read.delim(file),
                "xlsx"  = read.xlsx(file, sheetIndex = 1))
  names(data) <- c("x", "y")
  
  return(data)
}



# Classify data -----------------------------------------------------------


# Classify PCA + cluster  -----
.classify_fq <- function(data, ...){
  
  pca_classify   <- prcomp(data, 2)
  print(summary(pca_classify))
  clusters       <- hclust(dist(pca_classify$x, method = 'euclidean'), method = 'ward.D')
  mature.means   <- cutree(clusters, 2) - 1
  
  return(mature.means)
}

# Classify BayesianPCA + cluster  -----
.classify_bayes <- function(data, ...){
  
  bpca_classify  <- sim.bPCA(data, n.chains = 3, n.iter = 100000, n.burnin = 500)
  scores.chain   <- get.scores.chain.bPCA(bpca_classify, log(data))
  scores.bpca    <- summary.scores.bPCA(scores.chain, axes.to.get = 1:2)
  clusters       <- hclust(dist(as.data.frame(scores.bpca[2]), method = 'euclidean'), method = 'ward.D')
  mature.means   <- cutree(clusters, 2) - 1

  return(mature.means)  
}


classify_mature <- function(data, method = "fq", ...){
  
  data <- data[complete.cases(data), ]
  data <- log(data)
  
  mature_classify <- switch(method,
                           "fq"    = .classify_fq(data),
                           "bayes" = .classify_bayes(data))
  
  base           <- data.frame(data, mature.binom = mature_classify)
  # linear o quadratic discriminarion analysis
  test.cov.mat   <- boxM(base[, c("x", "y")], base[, "mature.binom"])
  if(test.cov.mat$p.value > 0.05){
    dis.reg      <- lda(mature.binom ~ ., data = base)
  }else{
    dis.reg      <- qda(mature.binom ~ ., data = base)
  }
  mature         <- as.numeric(as.character(predict(dis.reg)$class))
  output         <- data.frame(exp(base[, c("x", "y")]), mature)
  cat("number in juveline and adult group =", as.numeric(table(mature)[1]), ",", as.numeric(table(mature)[2]), "\n")
  class(output) <- "classify"

  return(output)
}



# Plot classify  ----------------------------------------------------------

plot.classify <- function(data, xlab = "X", ylab = "Y", colors = c(1, 2), pch = c(4, 5), ...){
  
  data <- data.frame(do.call("cbind", data))
  
  juv <- data[data$mature == 0, ]
  adt <- data[data$mature == 1, ] 
  fit.juv <- glm(y ~ x, data = juv)
  fit.adt <- glm(y ~ x, data = adt)
  eq.juv  <- paste0("Y = ", round(as.numeric(coef(fit.juv)[1]), 2), " + ", round(as.numeric(coef(fit.juv)[2]),2), " *X", sep = "")
  eq.adt  <- paste0("Y = ", round(as.numeric(coef(fit.adt)[1]), 2), " + ", round(as.numeric(coef(fit.adt)[2]),2), " *X", sep = "")
  pch <- ifelse (data$mature == 0, pch[1], pch[2])
  col <- ifelse (data$mature == 0, colors[1], colors[2])
  plot(data[, c("x", "y")], col = col, xlab = xlab, ylab = ylab, pch = pch)
  lines(juv$x, predict(fit.juv), col = colors[1], lwd = 2)
  lines(adt$x, predict(fit.adt), col = colors[2], lwd = 2)
  legend("topleft", c(paste("Juveniles: ", eq.juv), paste("Adults: ", eq.adt)), 
         bty = "n", pch = pch, col = colors, cex = 0.8)
  invisible()
}



# Calculate ogive ---------------------------------------------------------

.calculate_ogive_fq <- function(data, ...){
  
  model.glm <- glm(mature ~ x, data = data, family = binomial(link = "logit"))
  A <- as.numeric(model.glm$coef[1])
  B <- as.numeric(model.glm$coef[2])
  params <- c(A, B)
  
  pred.dat <- predict(model.glm, newdata = data.frame(x = data$x), se.fit = TRUE)
  fitted   <- round(model.glm$fitted, 3)
  lower    <- round(with(pred.dat, exp(fit - 1.96 * se.fit) / (1+exp(fit - 1.96 * se.fit))), 3)
  upper    <- round(with(pred.dat, exp(fit + 1.96 * se.fit) / (1+exp(fit + 1.96 * se.fit))), 3)
  
  estimate <- list(params = params, lower = lower, fitted = fitted,  upper = upper)

  return(estimate)
}


.calculate_ogive_bayes <- function(data, ...){
  
  model.bayes <- MCMClogit(data$mature ~ data$x, burnin = 1000, mcmc = 100000, thin = 10)
  stats       <- summary(model.bayes)
  A <- stats$quantiles[5]
  B <- stats$quantiles[6]
  params <- c(A, B)
  
  X.bayes     <- cbind(1, data$x)
  Xb          <- as.matrix(X.bayes) %*% t(model.bayes)
  pred.bayes  <- as.data.frame(1 / (1 + exp(-Xb)))
  qtl         <- round(rowQuantiles(pred.bayes, probs = c(0.025, 0.5, 0.975)), 3)
  fitted      <- qtl[, 2]
  lower       <- qtl[, 1]
  upper       <- qtl[, 3]
  
  estimate    <- list(params = params, lower = lower, fitted = fitted,  upper = upper)
  
  return(estimate)
}


calculate_ogive <- function(data, method = "fq", ...){
  
  data <- data.frame(do.call("cbind", data))
  estimate <- switch(method,
                     "fq" = .calculate_ogive_fq(data),
                     "bayes" = .calculate_ogive_bayes(data)) 
  
  output  <- data.frame(x = data$x, y = data$y, mature = data$mature, CIlower = estimate$lower, 
                        fitted = estimate$fitted, CIupper = estimate$upper)
  data    <- list(params = estimate$params, output = output)
  class(data) <- "ogive"
  
  return(data)
}


# Plot ogive --------------------------------------------------------------

plot.ogive <- function(data, xlab = "X", ylab = "Proportion mature", col1 = "blue", col2 = "red", ...){
  
  fit     <- data$output
  x.input <- fit$x
  y.input <- fit$mature
  m.p     <- tapply(y.input, x.input, mean)
  
  plot(sort(unique(x.input)), m.p, xlab = xlab, ylab = "Proportion mature", pch = 19, col = "darkgrey", axes = F)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  
  wide50 <- - data$params[1] / data$params[2]
  wide95 <- (1 / data$params[2]) * log(1 / 0.05 - 1) - data$params[1] / data$params[2]
  
  lines(sort(x.input), sort(fit$fitted), col = col1, lwd = 2)
  lines(sort(x.input), sort(fit$CIlower), col = col1, lwd = 2, lty = 2)
  lines(sort(x.input), sort(fit$CIupper), col = col1, lwd = 2, lty = 2)
  lines(c(wide50, wide50), c(-1, 0.5), col = col2, lty = 2, lwd = 2)
  lines(c(-1, wide50), c(0.5, 0.5), col = col2, lty = 2, lwd = 2)
  legend("topleft", as.expression(bquote(bold(CW[50] == .(round(wide50, 1))))), bty = "n")
  cat("Carapace width of 50% maturity =", round(wide50, 1), "\n")
  cat("Carapace width of 95% maturity =", round(wide95, 1), "\n")
  invisible()
}

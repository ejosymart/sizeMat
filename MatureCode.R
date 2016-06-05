#Read Data
dataMature = function(file, ext = "txt", ...){
  file = paste(file, ".", ext, sep="")
  
  data = switch(ext,
                "txt"   = read.table(file, header = TRUE, stringsAsFactors = FALSE),
                "csv"   = read.csv(file),
                "csv2"  = read.csv2(file),
                "delim" = read.delim(file),
                "xlsx"  = read.xlsx(file, sheetIndex = 1))
  names(data) = c("carapace", "chela")
  data = data[complete.cases(data), ]
  return(data)
}


#Mature
classifyByDistance = function(data, xlab = "Log(carapace width)", ylab = "Log(chela dimension)", 
                              col.juv = "green", col.adlt = "red", col.indt = "blue", n.iter = 100, ...){
  y  = data$chela
	x  = data$carapace
	ly = log(y)
	lx = log(x)
	plot(lx, ly, xlab = xlab, ylab = ylab, main = "Upper and lower bounds")
	
	# enter upper and lower bounds
	cat("enter lower bound ")
	xlow     = locator(1)$x
	cat(xlow, "\n")
	abline(v = xlow)
	cat("enter upper bound ")
	xup      = locator(1)$x
	cat(xup, "\n")
	abline(v = xup)
	ly.up    = ly[lx > xup]
	lx.up    = lx[lx > xup]
	ly.low   = ly[lx < xlow]
	lx.low   = lx[lx < xlow]
	ly.mid   = ly[lx >= xlow & lx <= xup]
	lx.mid   = lx[lx >= xlow & lx <= xup]
	plot(lx, ly, type = "n", xlab = xlab, ylab = ylab, main = "Determine initial fit to upper and lower data")
	points(lx.low, ly.low, col = col.juv, pch = 4)
	points(lx.up, ly.up, col = col.adlt, pch = 5)
	points(lx.mid, ly.mid, col = col.indt)
	abline(v = xlow)
	abline(v = xup)

	# determine initial fit to upper and lower data
	fit.up = lsfit(lx.up, ly.up)
	
	fit.low = lsfit(lx.low, ly.low)
	
	# iterative reassign mid points to upper or lower group
	# based on absolute residuals (distance)
	plot(lx, ly, xlab = xlab, ylab = ylab, type = "n",
	     main = "Iterative reassign mid points to upper or lower group \n based on absolute residuals (distance)")
	cat("\n", "Reassign points until min residual ss is found \n")
	for(i in seq_len(n.iter)){
		yp.mid.up  = fit.up$coef[1] + fit.up$coef[2] * lx.mid
		res.up     = abs(ly.mid - yp.mid.up)
		
		yp.mid.low = fit.low$coef[1] + fit.low$coef[2] * lx.mid
		res.low    = abs(ly.mid - yp.mid.low)
		
		#Points tentatively assigned to the closest line
		y.up  = c(ly.mid[res.up <= res.low], ly.up) 
		x.up  = c(lx.mid[res.up <= res.low], lx.up)
		y.low = c(ly.mid[res.low < res.up], ly.low)
		x.low = c(lx.mid[res.low < res.up], lx.low)
		cat("number in upper and lower group =", length(y.up), ",", length(y.low), "\n")
		
		fit.up  = glm(y.up ~ x.up)
		rss1    = sum((fit.up$resid)^2)
		fit.low = glm(y.low ~ x.low)
		rss2    = sum((fit.low$resid)^2)
		rss     = rss1 + rss2
		cat("Residual ss of combined fit =", rss, "\n")
		
		points(x.low, y.low, col = col.juv, pch = 4)
		points(x.up, y.up, col = col.adlt, pch = 5)
		lines(x.low, predict(fit.low), col = col.juv, lwd = 2)
		lines(x.up, predict(fit.up), col = col.adlt, lwd = 2)
	}

	output = data.frame(carapace = round(exp(c(x.up, x.low))), 
	                    chela    = round(exp(c(y.up, y.low))), 
	                    mature   = c(rep(1, length(y.up)), rep(0, length(y.low)))) 
# 	class(output) = "classify"
	return(output)
}


classifyCluster = function(data, xlab = "Log(carapace width)", ylab = "Log(chela dimension)", ...){

  #Principal Components Analysis and cluster
  pca.classify   = prcomp(log(data), 2)
  k.means.scores = kmeans(pca.classify$x, centers = 2, nstart = 20)
    
  clusters = hclust(dist(pca.classify$x, method = 'euclidean'), method = 'ward.D')
  clusterCut = cutree(clusters, 2)
  
#   pch=ifelse(k.means.scores$cluster == 1, 4, 5)
#   col=ifelse(k.means.scores$cluster == 1, 3, 2)
  pch=ifelse(clusterCut == 1, 4, 5)
  col=ifelse(clusterCut == 1, 3, 2)
  plot(log(data), col = col, xlab = xlab, ylab = ylab, pch = pch)
  
  # New data (classification)
#   mature = k.means.scores$cluster - 1  
  mature = clusterCut - 1  
  output = data.frame(data, mature = mature)
  output
#   class(output) = "classify"
  return(output)
}


ogive = function(data, ...){
  x.input = data$carapace
  y.input = data$mature
  # Model, params and deviance
  out = glm(y.input ~ x.input, family = binomial(link = "logit"))
  A = as.numeric(out$coef[1]); B = as.numeric(out$coef[2]); params = c(A, B)
  deviance = (out$null.deviance - out$deviance)/out$null.deviance
  # Predict ans Confidence intervals
  pred.dat = predict(out, newdata = data.frame(x.input = x.input), se.fit = TRUE)
  upper    = with(pred.dat, exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)))
  lower    = with(pred.dat, exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)))
  CI       = data.frame(lower, upper)
  # Output
  output  = data.frame(carapace = x.input, mature = y.input, fitted = round(out$fitted, 3))
  data    = list(model = out, deviance = deviance, params = params, output = output, confidence_intervals = CI)
  class(data) = "ogive"
  return(data)
}


ogiveBayes = function(data, ...){
  x.input = data$carapace
  y.input = data$mature
  # Model and params
  model.bayes = MCMClogit(y.input ~ x.input, burnin = 1000, mcmc = 100000, b0 = 0.1, B0 = 0.0001, thin = 10)
  stats       = summary(model.bayes)
  A = stats$quantiles[5]; B = stats$quantiles[6]; params = c(A, B)
  # Predict bayes and credible intervals
  X.bayes     = cbind(1, x.input)
  Xb          = as.matrix(X.bayes) %*% t(model.bayes)
  pred.bayes  = as.data.frame(1/(1 + exp(-Xb)))
  qtl         = rowQuantiles(pred.bayes, probs=c(0.025,0.5,0.975))
  CI          = data.frame(lower= qtl[,1], upper= qtl[,3])
  #Output
  output = data.frame(carapace = x.input, mature = y.input, fitted = round(qtl[, 2], 3))
  data   = list(params = params, output = output, confidence_intervals = round(CI, 3))
  class(data) = "ogive"
  return(data)
}


plot.ogive = function(data, xlab = "Width", ylab = "Proportion mature", col1 = "blue", col2 = "red", ...){
  
  fit     = data$output
  x.input = fit$carapace
  y.input = fit$mature
  m.p     = tapply(y.input, x.input, mean)
  ci      = data$confidence_intervals
  #Plots
  plot(sort(unique(x.input)), m.p, xlab = xlab, ylab = "Proportion mature", pch = 19, col = "darkgrey", axes = F)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  #Length at which the 50% of population is mature
  wide50 = -data$params[1]/data$params[2]
  #Length at which the 95% of population is mature
  wide95 = (1/data$params[2])*log(1/0.05-1)-data$params[1]/data$params[2]
  #Adding lines, intervlas and legend
  lines(sort(x.input), sort(fit$fitted), col = col1, lwd = 2)
  lines(sort(x.input), sort(ci$lower), col = col1, lwd = 2, lty = 2)
  lines(sort(x.input), sort(ci$upper), col = col1, lwd = 2, lty = 2)
  lines(c(wide50, wide50), c(-1, 0.5), col = col2, lty = 2, lwd = 2)
  lines(c(-1, wide50), c(0.5, 0.5), col = col2, lty = 2, lwd = 2)
  legend("topleft", as.expression(bquote(bold(CW[50] == .(round(wide50, 1))))), bty = "n")
  #Results
  cat("Carapace width of 50% maturity =", round(wide50, 1), "\n")
  cat("Carapace width of 95% maturity =", round(wide95, 1), "\n")
  invisible()
}

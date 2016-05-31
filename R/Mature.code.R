#Read Data
dataMature = function(file, ext = "txt", ...){
  file = paste(file, ".", ext, sep="")
  
  data = switch(ext,
                "txt"   = read.table(file, header = TRUE, stringsAsFactors = FALSE),
                "csv"   = read.csv(file),
                "csv2"  = read.csv2(file),
                "delim" = read.delim(file),
                "xlsx"  = read.xlsx(file, sheetIndex = 1))
  names(data) = c("x", "y")
  
  return(data)
}


#Mature
matureFit = function(data, col.juv = "green", col.adlt = "red", col.indt = "blue", n.iter = 100, ...)
{
	y = data$y
	x = data$x
	ly = log(y)
	lx = log(x)
	plot(lx, ly, xlab = "Log(carapace width)", ylab = "Log(chela dimension)",
	     main = "Upper and lower bounds")
	
	# enter upper and lower bounds
	cat("enter lower bound ")
	xlow = locator(1)$x
	cat(xlow, "\n")
	abline(v = xlow)
	cat("enter upper bound ")
	xup = locator(1)$x
	cat(xup, "\n")
	abline(v = xup)
	ly.up = ly[lx > xup]
	lx.up = lx[lx > xup]
	ly.low = ly[lx < xlow]
	lx.low = lx[lx < xlow]
	ly.mid = ly[lx >= xlow & lx <= xup]
	lx.mid = lx[lx >= xlow & lx <= xup]
	plot(lx, ly, type = "n", xlab = "Log(carapace width)", ylab = "Log(chela dimension)",
	     main = "Determine initial fit to upper and lower data")
	points(lx.low, ly.low, col = col.juv, pch = 4)
	points(lx.up, ly.up, col = col.adlt, pch = 5)
	points(lx.mid, ly.mid, col = col.indt)
	abline(v = xlow)
	abline(v = xup)

	# determine initial fit to upper and lower data
	fit.up = lsfit(lx.up, ly.up)
	# abline(fit.up)
	
	fit.low = lsfit(lx.low, ly.low)
	# abline(fit.low)

	# iterative reassign mid points to upper or lower group
	# based on absolute residuals (distance)
	plot(lx, ly, xlab = "Log(carapace width)", ylab = "Log(chela dimension)", type = "n",
	     main = "Iterative reassign mid points to upper or lower group \n based on absolute residuals (distance)")
	cat("\n", "Reassign points until min residual ss is found \n")
	for(i in seq_len(n.iter)){
		yp.mid.up = fit.up$coef[1] + fit.up$coef[2] * lx.mid
		res.up = abs(ly.mid - yp.mid.up)
		
		yp.mid.low = fit.low$coef[1] + fit.low$coef[2] * lx.mid
		res.low = abs(ly.mid - yp.mid.low)
		
		y.up = c(ly.mid[res.up <= res.low], ly.up) #choose the points that have a short distance to the predict line 
		x.up = c(lx.mid[res.up <= res.low], lx.up)
		y.low = c(ly.mid[res.low < res.up], ly.low)
		x.low = c(lx.mid[res.low < res.up], lx.low)
		cat("number in upper and lower group =", length(y.up), length(y.low), "\n")
		
		fit.up = glm(y.up ~ x.up)
		rss1 = sum((fit.up$resid)^2)
		fit.low = glm(y.low ~ x.low)
		rss2 = sum((fit.low$resid)^2)
		points(x.low, y.low, col = col.juv, pch = 4)
		points(x.up, y.up, col = col.adlt, pch = 5)
		lines(x.low, predict(fit.low), col = col.juv, lwd = 2)
		lines(x.up, predict(fit.up), col = col.adlt, lwd = 2)
		cat("Residual ss of combined fit =", rss1 + rss2, "\n")
	}

	output = data.frame(x = round(exp(c(x.up, x.low))), 
	                    y = c(y.up, y.low), 
	                    mature = c(rep(1, length(y.up)), rep(0, length(y.low)))) 
	class(output) = "mature"
	return(output)
}


plot.mature = function(output, xlab = "Width", ylab = "Proportion mature", col1 = "blue", col2 = "red", ...){
  
  x.in = output$x
  p.mat = output$mature
  out = glm(p.mat ~ x.in, family = "binomial")
  A = out$coef[1]; B = out$coef[2]
  m.p = tapply(p.mat, x.in, mean)
  plot(sort(unique(x.in)), m.p, xlab = xlab, ylab = "Proportion mature", pch = 19, col = "darkgrey", axes = F)
  axis(1)
  axis(2, seq(0, 1, 0.1), las = 1)
  box()
  lines(sort(x.in), sort(out$fitted), col = col1, lwd = 2)
  wide50 = -A/B
  cat("Carapace width of 50% maturity =", wide50, "\n")
  lines(c(wide50, wide50), c(-1, 0.5), col = col2, lty = 2, lwd = 2)
  lines(c(-1, wide50), c(0.5, 0.5), col = col2, lty = 2, lwd = 2)
  invisible()
}
  
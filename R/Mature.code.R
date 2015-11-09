#Mature
#To use the following code in R, use the "source" function after removing all of this header material.

matureFit = function(data)
{
	y = data$y
	x = data$x
	ly = log(y)
	lx = log(x)
	plot(lx, ly, xlab = "Log(carapace width)", ylab = "Log(chela dimension)")
	#
	# select out "bad" data (for particular data in example):  TO DO
# 	ycut.up = -2.7 + 1.33 * lx
# 	ycut.low = -3.5 + 1.41 * lx
# 	abline(c(-2.7, 1.33))
# 	abline(c(-3.5, 1.41))
# 	indx = rep(0, length(y))
# 	indx[ly < ycut.up & ly > ycut.low] = 1
# 	lx <- lx[indx == 1]
# 	ly <- ly[indx == 1]
# 	plot(lx, ly, xlab = "Log(carapace width)", ylab = "Log(chela dimension)")
	#
	# enter upper and lower bounds as in Mature 1
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
	#
	#determine initial fit to upper and lower data
	fit.up = lsfit(lx.up, ly.up)
	abline(fit.up)
	fit.low = lsfit(lx.low, ly.low)
	abline(fit.low)
	#
	#iterative reassign mid points to upper or lower group
	# based on absolute residuals (distance)
	cat("\n", "Reassign points until min residual ss is found \n")
	for(i in 1:100) {
		yp.mid.up = fit.up$coef[1] + fit.up$coef[2] * lx.mid
		res.up = abs(yp.mid.up - ly.mid)
		yp.mid.low = fit.low$coef[1] + fit.low$coef[2] * lx.mid
		res.low = abs(yp.mid.low - ly.mid)
		y.up = c(ly.mid[res.up <= res.low], ly.up)
		x.up = c(lx.mid[res.up <= res.low], lx.up)
		y.low = c(ly.mid[res.up > res.low], ly.low)
		x.low = c(lx.mid[res.up > res.low], lx.low)
		cat("number in upper and lower group =", length(y.up), length(y.low), "\n")
		fit.up = lsfit(x.up, y.up)
		rss1 = sum((fit.up$resid)^2)
		fit.low = lsfit(x.low, y.low)
		rss2 = sum((fit.low$resid)^2)
		abline(fit.up, col = 5)
		abline(fit.low, col = 6)
		cat("Residual ss of combined fit =", rss1 + rss2, "\n")
	}
	cat("\n")
	p.mat = c(rep(1, length(y.up)), rep(0, length(y.low)))
	x.in = round(exp(c(x.up, x.low)))
	m.p = tapply(p.mat, x.in, mean)
	plot(sort(unique(x.in)), m.p, ylab = "Proportion mature", xlab = "Width (mm)", pch = 19, col = "darkgrey")
	#assign("p.mat", p.mat, frame = 1)
	#assign("x.in", x.in, frame = 1)
	out = glm(p.mat ~ x.in, family = "binomial")
	lines(x.in[order(x.in)], out$fitted[order(x.in)], col = "blue", lwd = 2)
	wide50 = -out$coef[1]/out$coef[2]
	cat("Carapace width of 50% maturity=", wide50, "\n")
	lines(c(wide50, wide50), c(0, 0.5), col = "red", lty = 2, lwd = 2)
	lines(c(-1, wide50), c(0.5, 0.5), col = "red", lty = 2, lwd = 2)
	invisible()
}

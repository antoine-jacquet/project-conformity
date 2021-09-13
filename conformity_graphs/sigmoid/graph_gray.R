	
png("~/Documents/Cours/Toulouse School of Economics/PhD papers/conformity/images/sigmoid/sigmoid_step_gray.png")	
	
	layout(1)
	par(mgp=c(2.5,.5,0), mar=c(4,4,1,1)+0.1)
	
	# Parameters
	s0 <- .2 ;							# value at 0 and 1
	s <- 3 ;								# slope at .5
	xlow <- (s + 2*s0 -1)/(2*s) ;		# end of the first plateau
	xhigh <- 1 - xlow ;					# beginning of the second plateau
	id <- function(x) {x} ;				# identity function
	
	# Plot frame with random copying as baseline
	plot(id, 0, 1, lwd=3, lty=2,						# plot of the identity function (random copying)
		xaxp=c(0,1,2), yaxp=c(0,1,2), tck=.02,
		cex.lab=1.5,
		xlab="Proportion of demonstrators displaying A = 'Teachers'",
		ylab="Proportion of observers adopting A = 'Pupils'") ;
	abline(v = 1/2, lty="dashed") ;		# vertical line at x = .5
	
	# Fill conformity areas
	polygon(x = c(0, 1/2, 1/2), y = c(0, 0, 1/2), col=adjustcolor("gray", alpha.f=.3), border=NA)
	polygon(x = c(1, 1/2, 1/2), y = c(1, 1, 1/2), col=adjustcolor("gray", alpha.f=.3), border=NA)
	
	# Sigmoid curves
	for (k in 1:5) {
		f1 <- function(x) x + k/5*x*(1-x)*(2*x-1) - s0
		x1 <- uniroot(f1, c(0, 1/2))$root
		f2 <- function(x) x + k/5*x*(1-x)*(2*x-1) - (1 - s0)
		x2 <- uniroot(f2, c(1/2, 1))$root
		if (k==5) width <- 3 else width <- 1
		curve(x + k/5*x*(1-x)*(2*x-1), x1, x2, col="gray50", lwd=width, add=TRUE);
		curve(x + k/5*x*(1-x)*(2*x-1), 0, x1, col="gray50", lwd=width, lty=1, add=TRUE);
		curve(x + k/5*x*(1-x)*(2*x-1), x2, 1, col="gray50", lwd=width, lty=1, add=TRUE);
		curve(x + k/5*x*(1-x)*(2*x-1), 0, x1, col="white", lwd=width, lty=3, add=TRUE);
		curve(x + k/5*x*(1-x)*(2*x-1), x2, 1, col="white", lwd=width, lty=3, add=TRUE);
	}
	
	# Quasi-step function
	curve(1/2 + s*(x-1/2), from = xlow, to= xhigh, col="black", lwd=3, add=TRUE) ;
	curve(s0*x/x, from=0, to=xlow, col="black", lwd=3, add=TRUE) ;
	curve((1-s0)*x/x, from=xhigh, to=1, col="black", lwd=3, add=TRUE) ;
	
	legend(-.02, .95, bty="n",
		legend=c("Random copying", "Sigmoid curves", "Parts unobserved empirically", "Step function"),
		col=c("black", "gray50", "gray50", "black"),
		lwd=c(3,3,3,3),
		lty=c(2,1,1,1)
		)
	legend(-.02, .95, bty="n",
		legend=c("", "", ""),
		col=c(NA,NA,"white"),
		lwd=c(NA,NA,3.1),
		lty=c(NA,NA,3)
		)
		
	legend(-.012, .8, bty="n",		
		legend="",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)
	legend(.010, .8, bty="n",		# x+0.22
		legend="Conformity area",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)

dev.off()
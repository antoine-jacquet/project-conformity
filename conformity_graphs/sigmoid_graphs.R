


###################################
### GRAPHS FOR CONFORMITY PAPER ###
###################################

q <- seq(0, 1, length=5000)

### GRAPH 1a: BOYD & RICHERSON CONFORMITY CURVES

pdf("~/Documents/GitHub/project-conformity/conformity_graphs/fig1a.pdf")

	layout(1)
	par(mgp=c(2.5,.5,0), mar=c(4,4,1,1)+0.1)

	# Plot frame with random copying as baseline
	plot(q, q, type="n",						# random copying
		xaxp=c(0,1,2), yaxp=c(0,1,2), tck=.02,
		cex.lab=1.5,
		xlab="Proportion of demonstrators displaying A = 'Teachers'",
		ylab="Proportion of observers adopting A = 'Pupils'")
	abline(v = 1/2, lty="dashed")		# vertical line at x = .5
	
	# Fill conformity areas
	polygon(x = c(0, 1/2, 1/2), y = c(0, 0, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	polygon(x = c(1, 1/2, 1/2), y = c(1, 1, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	
	# Sigmoid curves
	B <- 1:4/5
	for (β in B) lines(q, q + β*q*(1-q)*(2*q-1), col="red", lwd=1)
	lines(q, q + 0*q*(1-q)*(2*q-1), col="red", lwd=2)
	lines(q, q + 1*q*(1-q)*(2*q-1), col="red", lwd=2)
	
	# Other example of conformist bias
	lines(q, q^5/(q^5 + (1-q)^5), col="blue", lwd=2)
	
	# Legend
	legend(-.02, .95, bty="n",
		legend=c("Boyd-Richerson sigmoid curves", "Other example of conformist bias"),
		col=c("red", "blue"),
		lwd=c(3,3),
		lty=c(1,1)
		)
	legend(-.009, .88, bty="n",		
		legend="",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)
	legend(.013, .88, bty="n",		# x+0.22
		legend="Boyd-Richerson conformity area",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)
	text(.25, .25, "b = 0", col="red", pos=2)
	text(.25, .15, "b = 1", col="red", pos=4)

dev.off()

	
### GRAPH 1b: LOGISTIC CONFORMITY CURVES

pdf("~/Documents/GitHub/project-conformity/conformity_graphs/fig1b.pdf")

	layout(1)
	par(mgp=c(2.5,.5,0), mar=c(4,4,1,1)+0.1)

	# Plot frame with random copying as baseline
	plot(q, q, type="n",						# random copying
		xaxp=c(0,1,2), yaxp=c(0,1,2), tck=.02,
		cex.lab=1.5,
		xlab="Proportion of demonstrators displaying A = 'Teachers'",
		ylab="Proportion of observers adopting A = 'Pupils'")
	abline(v = 1/2, lty="dashed")		# vertical line at x = .5
	abline(h = 1/2, lty="dashed")		# horizontal line at y = .5
	
	# Fill conformity areas
	polygon(x = c(0, 1/2, 1/2), y = c(0, 0, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	polygon(x = c(1, 1/2, 1/2), y = c(1, 1, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	
	# Sigmoid curves
	B <- c(1.5, 3, 5, 10) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="red", lwd=1)
	β  <- 1 ; lines(q, q^β/(q^β + (1-q)^β), col="black", lwd=2)
	B <- 1/c(1.5, 3, 5, 10) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="orange", lwd=1)
	B <- -c(1/c(1.5, 3, 5, 10), c(1.5, 3, 5, 10)) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="pink", lwd=1)
	
	# Legend
	legend(-.02, .95, bg="white", box.col="white",
		legend=c("β > 1", "0 < β < 1", "β < 0"),
		col=c("red", "orange", "pink"),
		lwd=c(2,2,2),
		lty=c(1,1,1)
		)

dev.off()


### FIGURE 1: BOTH GRAPHS SIDE BY SIDE

png("~/Documents/GitHub/project-conformity/conformity_graphs/fig1.png", width=2*480, height=480)

	layout(matrix(1:2, nrow=1, ncol=2, byrow=T))
	par(mgp=c(2.5,.5,0), mar=c(4,4,3,3)+0.1)
	
	# Plot frame with random copying as baseline
	plot(q, q, type="n",						# random copying
		xaxp=c(0,1,2), yaxp=c(0,1,2), tck=.02,
		cex.lab=1.4,
		xlab="Proportion of demonstrators displaying A = 'Teachers'",
		ylab="Proportion of observers adopting A = 'Pupils'")
	abline(v = 1/2, lty="dashed")		# vertical line at x = .5
	
	# Fill conformity areas
	polygon(x = c(0, 1/2, 1/2), y = c(0, 0, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	polygon(x = c(1, 1/2, 1/2), y = c(1, 1, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	
	# Sigmoid curves
	B <- 1:4/5
	for (β in B) lines(q, q + β*q*(1-q)*(2*q-1), col="red", lwd=1)
	lines(q, q + 0*q*(1-q)*(2*q-1), col="red", lwd=2)
	lines(q, q + 1*q*(1-q)*(2*q-1), col="red", lwd=2)
	
	# Other example of conformist bias
	lines(q, q^5/(q^5 + (1-q)^5), col="blue", lwd=2)
	
	# Legend
	legend(-.02, .96, bty="n",
		legend=c("Boyd-Richerson sigmoid curves", "Other example of conformist bias"),
		col=c("red", "blue"),
		lwd=c(3,3),
		lty=c(1,1)
		)
	legend(-.009, .88, bty="n",		
		legend="",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)
	legend(.013, .88, bty="n",		# x+0.22
		legend="Boyd-Richerson conformity area",
		fill=adjustcolor("gray", alpha.f=.5),
		border=NA
		)
	text(.25, .25, "β = 0", col="red", pos=2)
	text(.25, .15, "β = 1", col="red", pos=4)
	title("A", line=2)
	
	# Plot frame with random copying as baseline
	plot(q, q, type="n",						# random copying
		xaxp=c(0,1,2), yaxp=c(0,1,2), tck=.02,
		cex.lab=1.4,
		xlab="",
		ylab="")
	abline(v = 1/2, lty="dashed")		# vertical line at x = .5
	abline(h = 1/2, lty="dashed")		# horizontal line at y = .5
	
	# Fill conformity areas
	polygon(x = c(0, 1/2, 1/2), y = c(0, 0, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	polygon(x = c(1, 1/2, 1/2), y = c(1, 1, 1/2), col=adjustcolor("gray", alpha.f=.5), border=NA)
	
	# Sigmoid curves
	B <- c(1.5, 3, 5, 10) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="red", lwd=1)
	β  <- 1 ; lines(q, q^β/(q^β + (1-q)^β), col="black", lwd=2)
	B <- 1/c(1.5, 3, 5, 10) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="orange", lwd=1)
	B <- -c(1/c(1.5, 3, 5, 10), c(1.5, 3, 5, 10)) ; for (β in B) lines(q, q^β/(q^β + (1-q)^β), col="pink", lwd=1)
	
	# Sigmoid curve with error rate
	α <- .3 ; β  <- 4
	lines(q, α/2 + (1-α)*q^β/(q^β + (1-q)^β), col="dark green", lwd=2, lty=4)
	
	# Legend
	legend(-.02, .95, bg="white", box.col="white",
		legend=c("β > 1", "0 < β < 1", "β < 0", "β > 1 with error rate"),
		col=c("red", "orange", "pink", "dark green"),
		lwd=c(2,2,2,2),
		lty=c(1,1,1,4)
		)
	title("B", line=2)

dev.off()





























	
	
	

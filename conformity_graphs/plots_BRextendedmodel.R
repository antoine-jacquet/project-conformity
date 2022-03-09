




p_prime <- function(p,n) { k <- 0:n ; p^(n+1) * sum( factorial(2*n+1)/factorial(n+1+k)/factorial(n-k)*p^(k)*(1-p)^(n - k) ) }

p_span <- (0:100)/100
f_plot <- p_span

plot(p_span, p_span, type="n")
for (n in c(3, 5, 7, 9, 11)) {
	for (p in p_span) f_plot[p*100+1] <- p_prime(p, n)
	lines(p_span, f_plot, col=n)
}

n <- 3

k <- 0:n

plot( factorial(n)/factorial(k), factorial(2*n+1)/factorial(n+1+k) )
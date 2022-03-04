rm(list=objects())

set.seed(2)

#### Simulation

n = 500
sigma = 0.5

phi = 0.7 # Choice, < 1

eps = rnorm(n, 0, sd=sigma)

x = arima.sim(list(ar = c(phi)), n=n, innov=eps)

plot(x, type="l", ylab=expression(X[t]))

#### Autocovariance estimation

## Autocorrelation

rho_estimated = acf(x, lag.max=10, plot=FALSE)
rho_estimated = rho_estimated$acf # acf returns a list with various attributes, among them:
# acf = values of the autocorrelation
# lags = lags corresponding to the values in $acf
rho_theoretical = phi^c(0:10)

plot(rho_estimated, type='b', pch=20, ylab=expression(rho(h)), xlab="Lag h")
lines(rho_theoretical, col='red')
legend("topright", legend=c(expression(rho[estimated]), expression(rho[theoretical])),
       lty=1, col=c("black","red"), pch=c(20, NA))

## Autocovariance 

gamma_estimated = rho_estimated*var(x)
gamma_theoretical = sigma^2/(1-phi^2)*phi^c(0:10)

plot(gamma_estimated, type='b', pch=20, ylab=expression(gamma(h)), xlab="Lag h")
lines(gamma_theoretical, col='red')
legend("topright", legend=c(expression(gamma[estimated]), expression(gamma[theoretical])),
       lty=1, col=c("black","red"), pch=c(20, NA))

#### Phi estimator

## One-shot

phi_hat = rho_estimated[2]

## Monte-Carlo

estim_phi = function(n,  sigma, phi)
{
  eps = rnorm(n, 0, sd=sigma)
  x = arima.sim(list(ar = c(phi)), n=n, innov=eps)
  rho_estimated = acf(x, lag.max=1, plot =FALSE)$acf
  phi_hat = rho_estimated[2]
  return(phi_hat)
}

N_simu = 500
phi_hat = sapply(rep(n, N_simu), estim_phi, sigma=sigma, phi=phi)

hist(phi_hat, breaks=30)
abline(v=phi, col='red')

alpha = 0.05

confidence_int = quantile(phi_hat, c(alpha/2,1-alpha/2))
abline(v=confidence_int[1], col='blue')
abline(v=confidence_int[2], col='blue')

# handy solution

a = sort(phi_hat)[floor(alpha/2*N_simu)]
b = sort(phi_hat)[floor((1-alpha/2)*N_simu)]
confidence_int = c(a,b)

N_simu = 100 # to go faster
list_n = seq(50, 500, by=50)
ICs = NULL
averaged_phi_hat = c()
for(n in list_n)
{
  phi_h = sapply(rep(n, N_simu), estim_phi, sigma=sigma, phi=phi)
  ICs = rbind(ICs, quantile(phi_h, c(alpha/2,1-alpha/2)))
  averaged_phi_hat = c(averaged_phi_hat, mean(phi_h))
}


plot(list_n, averaged_phi_hat, ylim = c(0.4,0.8), xlab="n", ylab=expression(hat(phi)))
arrows(list_n, ICs[,1], list_n, ICs[,2], length=0.05, angle=90, code=3)

sizes = abs(ICs[,2]-ICs[,1])
plot(list_n, sizes, type='b', pch=20, xlab="n", ylab="Interval length")

conv = 1/sqrt(list_n)
reg = lm(sizes~conv-1)
lines(list_n, reg$coeff/sqrt(list_n), col='red')
legend("topright", legend=c(expression(frac(alpha,sqrt(n)))), col="red", lty=1)

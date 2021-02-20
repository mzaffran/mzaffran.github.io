rm(list=objects())

#### Data simulation

set.seed(4) # ensures reproducible results

n = 100
t = 1:n

## AR(1)

# Generation

eps = rnorm(n,0,1)

ar1_t = function(t, a, epsilon){
  epsilon = head(epsilon, t)
  powers_a = a^((t):1)
  return(powers_a %*% epsilon) # * is done term-by-term. %*% is the matrix multiplication.
}

y1 = sapply(1:n, ar1_t, a=0.1, epsilon=eps)
y2 = sapply(1:n, ar1_t, a=0.7, epsilon=eps)
y3 = sapply(1:n, ar1_t, a=-0.7, epsilon=eps)

plot(t, y1, type='l', xlab="Time", ylab="AR(1)", ylim=range(y1,y2,y3))
lines(t, y2, col='blue')
lines(t, y3, col='red')
legend("bottomleft", legend=c("a = 0.1", "a = 0.7", "a = -0.7"),
       lty=1, col=c("black","blue","red"), cex=0.8)

# a=0.1 => less variance.
# a=-0.7 => more round trips around 0.

#### ACF and PACF functions

autoCorr = function(x,h)
{
  x.lag = lag.xts(x, k=h, na.pad=T)
  return(cor(x.lag, x, use="pairwise.complete.obs"))
}

PartialAutoCorr<-function(x,h)
{
  x.lag = lapply(c(1:h), lag.xts, x=x, na.pad=T)
  x.lag = matrix(unlist(x.lag), ncol=length(x.lag))
  reg = lm(x~x.lag-1)
  return(as.numeric(tail(reg$coef,1)))
}

#### Results on our data

par(mfrow=c(1,2)) # In the plot window, 1 line and 2 columns: 2 plots next to each other.

## a = 0.1

acf1 = c(1, sapply(c(1:30), autoCorr, x=y1))
plot(0:30, acf1, type='h', ylim=c(-0.3,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),0.1^c(0:30),col='red') 
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)
# ρ(h) =a^|h| for an AR of parameter a. We obtain coherent results with this.

pacf1 = sapply(c(1:30), PartialAutoCorr, x=y1)
plot(1:30, pacf1, type='h', ylim=c(-0.3,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y1, lag.max=30)
lines(c(0:30),0.1^c(0:30),col='red')

pacf(y1, lag.max=30)

## a = 0.7

acf2 = c(1, sapply(c(1:30), autoCorr, x=y2))
plot(0:30, acf2, type='h', ylim=c(-0.2,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),0.7^c(0:30),col='red')
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)

pacf2 = sapply(c(1:30), PartialAutoCorr, x=y2)
plot(1:30, pacf2, type='h', ylim=c(-0.2,0.8), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y2, lag.max=30)
lines(c(0:30),0.7^c(0:30),col='red')

pacf(y2, lag.max=30)

## a = -0.7

acf3 = c(1, sapply(c(1:30), autoCorr, x=y3))
plot(0:30, acf3, type='h', ylim=c(-0.7,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),(-0.7)^c(0:30),col='red')
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)

pacf3 = sapply(c(1:30), PartialAutoCorr, x=y3)
plot(1:30, pacf3, type='h', ylim=c(-0.7,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y3, lag.max=30)
lines(c(0:30),(-0.7)^c(0:30),col='red')

pacf(y3, lag.max=30)

# Generally, we have small differences with our handy functions, and the ones from R.
# It comes from the estimation procedure:
# - ACF uses the mean over the n-h observations and not all the observations
# - PACF fits autoregressive models of successively higher orders up to lag.max

# The exponential fit was not always really good (especially for a=0.7), 
# but we have really few data: thus we have more "noise" in comparison with 
# the processus structure itself.
# Increase the size!

par(mfrow=c(1,1))

n = 200
t = 1:n

eps = rnorm(n,0,1)

y1 = sapply(1:n, ar1_t, a=0.1, epsilon=eps)
y2 = sapply(1:n, ar1_t, a=0.7, epsilon=eps)
y3 = sapply(1:n, ar1_t, a=-0.7, epsilon=eps)

plot(t, y1, type='l', xlab="Time", ylab="AR(1)", ylim=range(y1,y2,y3))
lines(t, y2, col='blue')
lines(t, y3, col='red')
legend("bottomleft", legend=c("a = 0.1", "a = 0.7", "a = -0.7"),
       lty=1, col=c("black","blue","red"), cex=0.8)

par(mfrow=c(1,2)) 

acf1 = c(1, sapply(c(1:30), autoCorr, x=y1))
plot(0:30, acf1, type='h', ylim=c(-0.3,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),0.1^c(0:30),col='red') 
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)

pacf1 = sapply(c(1:30), PartialAutoCorr, x=y1)
plot(1:30, pacf1, type='h', ylim=c(-0.3,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y1, lag.max=30)
lines(c(0:30),0.1^c(0:30),col='red')

pacf(y1, lag.max=30)

acf2 = c(1, sapply(c(1:30), autoCorr, x=y2))
plot(0:30, acf2, type='h', ylim=c(-0.2,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),0.7^c(0:30),col='red')
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)

pacf2 = sapply(c(1:30), PartialAutoCorr, x=y2)
plot(1:30, pacf2, type='h', ylim=c(-0.2,0.8), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y2, lag.max=30)
lines(c(0:30),0.7^c(0:30),col='red')

pacf(y2, lag.max=30)

acf3 = c(1, sapply(c(1:30), autoCorr, x=y3))
plot(0:30, acf3, type='h', ylim=c(-0.7,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
lines(c(0:30),(-0.7)^c(0:30),col='red')
legend("topright", legend=c(expression(a^group("|",h,"|") %~~% rho(h))), 
       lty=1, col="red", cex=0.5)

pacf3 = sapply(c(1:30), PartialAutoCorr, x=y3)
plot(1:30, pacf3, type='h', ylim=c(-0.7,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y3, lag.max=30)
lines(c(0:30),(-0.7)^c(0:30),col='red')

pacf(y3, lag.max=30)

# The exponential speed of decrease is really better fitted!

par(mfrow=c(1,1))

#### Data generation: MA(4)

n = 100
t = 1:n

eps = rnorm(n,0,1)

a4 = c(1,rep(2, 3))
y4 = stats::filter(eps, filter=a4, method = c("convolution"), sides=1, circular=FALSE)[4:n]

a5 = c(1,rep(0.5, 3))
y5 = stats::filter(eps, filter=a5, method = c("convolution"), sides=1, circular=FALSE)[4:n]

a6 = c(1,rep(-0.5, 3))
y6 = stats::filter(eps, filter=a6, method = c("convolution"), sides=1, circular=FALSE)[4:n]

new_t = 4:n

plot(new_t, y4, type='l', xlab="Time", ylab="MA(4)", ylim=range(y4,y5,y6, na.rm=TRUE))
lines(new_t, y5, col='blue')
lines(new_t, y6, col='red')

par(mfrow=c(1,2))

acf4 = c(1, sapply(c(1:30), autoCorr, x=y4))
plot(0:30, acf4, type='h', ylim=c(-0.3,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)
# γ(q+h) = 0, for a MA(q). Here, we should have ACF near 0 starting from lag 5. 

pacf4 = sapply(c(1:30), PartialAutoCorr, x=y4)
plot(1:30, pacf4, type='h', ylim=c(-0.3,0.8), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y4, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y4, lag.max=30)

acf5 = c(1, sapply(c(1:30), autoCorr, x=y5))
plot(0:30, acf5, type='h', ylim=c(-0.3,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)

pacf5 = sapply(c(1:30), PartialAutoCorr, x=y5)
plot(1:30, pacf5, type='h', ylim=c(-0.3,0.7), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y5, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y5, lag.max=30)

acf6 = c(1, sapply(c(1:30), autoCorr, x=y6))
plot(0:30, acf6, type='h', ylim=c(-0.4,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)

pacf6 = sapply(c(1:30), PartialAutoCorr, x=y6)
plot(1:30, pacf6, type='h', ylim=c(-0.4,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y6, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y6, lag.max=30)

par(mfrow=c(1,1))

# The 0 coefficients for lags above 5 was clearly not the case.
# But we have a very small dataset, thus we have more "noise" in comparison with the
# processus structure itself.
# Increase the size!

n = 200
t = 1:n

eps = rnorm(n,0,1)

a4 = c(1,rep(2, 3))
y4 = stats::filter(eps, filter=a4, method = c("convolution"), sides=1, circular=FALSE)[4:n]

a5 = c(1,rep(0.5, 3))
y5 = stats::filter(eps, filter=a5, method = c("convolution"), sides=1, circular=FALSE)[4:n]

a6 = c(1,rep(-0.5, 3))
y6 = stats::filter(eps, filter=a6, method = c("convolution"), sides=1, circular=FALSE)[4:n]

new_t = 4:n

plot(new_t, y4, type='l', xlab="Time", ylab="MA(4)", ylim=range(y4,y5,y6, na.rm=TRUE))
lines(new_t, y5, col='blue')
lines(new_t, y6, col='red')

par(mfrow=c(1,2))

acf4 = c(1, sapply(c(1:30), autoCorr, x=y4))
plot(0:30, acf4, type='h', ylim=c(-0.2,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)

pacf4 = sapply(c(1:30), PartialAutoCorr, x=y4)
plot(1:30, pacf4, type='h', ylim=c(-0.3,0.8), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y4, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y4, lag.max=30)

acf5 = c(1, sapply(c(1:30), autoCorr, x=y5))
plot(0:30, acf5, type='h', ylim=c(-0.2,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)

pacf5 = sapply(c(1:30), PartialAutoCorr, x=y5)
plot(1:30, pacf5, type='h', ylim=c(-0.2,0.7), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y5, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y5, lag.max=30)

acf6 = c(1, sapply(c(1:30), autoCorr, x=y6))
plot(0:30, acf6, type='h', ylim=c(-0.4,1), xlab="Lag", ylab=expression(ACF ~~ (rho(h))))
abline(h=0)
abline(v=5, lty=3, col="red")
legend("topright", legend=c(expression(q+1: ~~ above ~~ it ~~ rho(h) %~~% 0)), 
       lty=3, col="red", cex=0.5)

pacf6 = sapply(1:30, PartialAutoCorr, x=y6)
plot(1:30, pacf6, type='h', ylim=c(-0.4,0.2), xlab="Lag", ylab=expression(PACF ~~ (r(h))))
abline(h=0)

acf(y6, lag.max=30)
abline(v=5, lty=3, col="red")

pacf(y6, lag.max=30)

# With a bigger data-set, we indeed have 0 coefficients for h >= 5.

par(mfrow=c(1,1))

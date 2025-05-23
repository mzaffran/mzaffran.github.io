rm(list=objects())

#####################
#### Data simulation
#####################

set.seed(4) # ensures reproducible results

n = 100
t = 1:n

eps = rnorm(n,0,1)

X1 = eps

X2 = t/5 + eps

s = cos(2*pi*t/10)
X3 = t/5 + s + eps

plot(X1, type='l', ylim=range(X1,X2,X3), xlab="Time", ylab=expression(X[t]))
# expression is a function to write maths expressions and greek letters. 
lines(X2, col='blue')
lines(X3, col='red')
legend("topleft", legend=c(expression(X^1),expression(X^2),expression(X^3)),
       lty=1, col=c("black","blue","red"))

#####################
#### Home made
#####################

#### Smoothing methods

# /!\ Important /!\
# Recall that ychap_(t+1) is obtained by taking the smoothed value of time t
# Smoothing at time t is forecasting time t+1

ExpSmooth = function(x,alpha)
{
  n = length(x)
  xsmooth = array(NA, dim=n) 
  
  # Initialization
  # We have to start somewhere: this would be optimized in R, using the Holt-Winters function.
  # Here, we arbitrarily choose a starting point, that is not absurd (at the beginning we do not
  # smooth the time series since we have no past).
  xsmooth[1] = x[1]
  
  for(i in 2:n)
  {
    xsmooth[i] = alpha*x[i] + (1-alpha)*xsmooth[i-1] 
  }
  
  return(xsmooth)
}

DoubleExpSmooth = function(x,alpha)
{
  n = length(x)
  xsmooth = array(NA, dim=n) 

  # Initialization
  # Again, it is here arbitrary but not absurd. l_t represents the level, thus it is in the same
  # rough size. b_t, instead, represents the trend: it is the same rough size and concept than
  # x[2]-x[1]. For l_1, I put x[2] to be at the "same point" than b_1, but I could have chosen
  # x[1] for example.
  # Anyway, usually these parameters are optimized, and we hope here that their influence will 
  # quickly be insignificant, after a warm up period.
  xsmooth[1] = x[1]
  l = array(NA, dim=n) 
  l[1] = x[2]
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  
  for(i in 2:n)
  {
    l[i] = xsmooth[i-1] + (1-(1-alpha)^2)*(x[i]-xsmooth[i-1]) 
    # when h=1, xsmooth[i-1] = l[i-1]+b[i-1]
    # xsmooth[i-1] represents xchap_i/(i-1)
    b[i] = b[i-1] + alpha^2*(x[i]-xsmooth[i-1])
    xsmooth[i] = l[i] + b[i] # same comment than 2 lines before
  }
  
  # the following lines of code allow to return multiple objects, stored in a list
  # return() does not allow to return multiple object (e.g. return(b,l) is not possible)
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  return(res)
}

DoubleExpSmoothHW=function(x,alpha,beta)
{
  n = length(x)
  xsmooth = array(NA, dim=n)
  
  xsmooth[1] = x[1]
  l = array(NA, dim=n)
  l[1] = x[2]
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  
  for(i in 2:n)
  {
    l[i] = alpha*x[i] + (1-alpha)*xsmooth[i-1]
    b[i] = beta*(l[i]-l[i-1]) + (1-beta)*b[i-1]
    xsmooth[i] = l[i] + b[i]
  }
  
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  return(res)
}

SeasonalDoubleExpSmooth = function(x,alpha,beta,delta,period)
{
  n = length(x)
  xsmooth = array(NA, dim=n)

  xsmooth[1] = x[1]
  l = array(NA, dim=n)
  l[1] = x[2] # x[1] ?
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  s = array(NA, dim=n)
  s[1] = x[1]
  
  for(i in 2:n)
  {
    l[i] = alpha*(x[i]-s[max(i-period,1)]) + (1-alpha)*(l[i-1]+b[i-1])
    b[i] = beta*(l[i]-l[i-1]) + (1-beta)*b[i-1]
    s[i] = delta*(x[i]-l[i]) + (1-delta)*s[max(i-period,1)]
    xsmooth[i] = l[i] + b[i] + s[i]
  }
  
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  res$s = s
  return(res)
}

#### Error functions

mse = function(x,y){
  return(mean((x-y)^2))
}

# Since we are going to forecast, we want to optimize our choice of alpha in this framework,
# and not in a smoothing framework. Thus, we build functions to evaluate the MSE of a 
# 1-day ahead forecast, given a smoothing => we have to shift the smoothing vector.

mse_smoothing_simple = function(true_values, smoothed_values){
  n = length(true_values)
  # we compare what we should have forecasted to what we smooth!
  return(mse(tail(true_values,n-1), head(smoothed_values,n-1))) 
}

mse_smoothing_methods = function(true_values, smoothed_object){
  n = length(true_values)
  return(mse(tail(true_values,n-1), head(smoothed_object$smooth,n-1)))
}

#### Application to our time series: exponential

## X1

alpha = 0.2
X1.smooth = ExpSmooth(X1,alpha)
plot(X1, type='l', xlab="Time", ylab=expression(X^1))
lines(X1.smooth, col='red', lwd=2)
legend("bottomleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, lwd=c(1,2), col=c("black","red"), cex=0.6)

mse_smoothing_simple(X1, X1.smooth)

alpha = seq(0.05,0.95,length=100)

forecast = lapply(alpha, ExpSmooth, x=X1)

erreur = sapply(forecast, mse_smoothing_simple, true_values=X1)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X1.smooth = ExpSmooth(X1, alpha[which.min(erreur)])
plot(X1, type='l', xlab="Time", ylab=expression(X^1))
lines(X1.smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)

## X2

forecast = lapply(alpha, ExpSmooth, x=X2)

erreur = sapply(forecast, mse_smoothing_simple, true_values=X2)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X2.smooth = ExpSmooth(X2, alpha[which.min(erreur)])
plot(X2, type='l', xlab="Time", ylab=expression(X^2))
lines(X2.smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)

## X3

forecast = lapply(alpha, ExpSmooth, x=X3)

erreur = sapply(forecast, mse_smoothing_simple, true_values=X3)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X3.smooth = ExpSmooth(X3, alpha[which.min(erreur)])
plot(X3, type='l', xlab="Time", ylab=expression(X^3))
lines(X3.smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)

#### Application to our time series: double exponential

## X1

forecast = lapply(alpha, DoubleExpSmooth, x=X1)

erreur = sapply(forecast, mse_smoothing_methods, true_values=X1)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X1.smooth = DoubleExpSmooth(X1, alpha[which.min(erreur)])
plot(X1, type='l', xlab="Time", ylab=expression(X^1))
lines(X1.smooth$smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)

plot(X1.smooth$l, type='l', ylim=range(X1.smooth$l,X1.smooth$b), col='blue',
     xlab='Time', ylab="Components of the double exponential smoothing")
lines(X1.smooth$b, col='red')
legend("topleft", legend=c(expression(l[t]), expression(b[t])),
       lty=1, col=c("blue","red"))

## X2

forecast = lapply(alpha, DoubleExpSmooth, x=X2)

erreur = sapply(forecast, mse_smoothing_methods, true_values=X2)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X2.smooth = DoubleExpSmooth(X2, alpha[which.min(erreur)])
plot(X2, type='l', xlab="Time", ylab=expression(X^2))
lines(X2.smooth$smooth, col='red')

plot(X2.smooth$l, type='l', ylim=range(X2.smooth$l,X2.smooth$b), col='blue',
     xlab='Time', ylab="Components of the double exponential smoothing")
lines(X2.smooth$b, col='red', type='l')
legend("topleft", legend=c(expression(l[t]), expression(b[t])),
       lty=1, col=c("blue","red"))

## X3

forecast = lapply(alpha, DoubleExpSmooth, x=X3)

erreur = sapply(forecast, mse_smoothing_methods, true_values=X3)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

X3.smooth = DoubleExpSmooth(X3, alpha[which.min(erreur)])
plot(X3, type='l', xlab="Time", ylab=expression(X^3))
lines(X3.smooth$smooth, col='red')

plot(X3.smooth$l, type='l', ylim=range(X3.smooth$l,X3.smooth$b), col='blue',
     xlab='Time', ylab="Components of the double exponential smoothing")
lines(X3.smooth$b, col='red', type='l')
legend("topleft", legend=c(expression(l[t]), expression(b[t])),
       lty=1, col=c("blue","red"))

#### Application to our time series: seasonal double exponential

# First draft of parameters: they should be optimized!
alpha = 0.35
beta = 0.2
delta = 0.2
period = 10
X3.seassmooth = SeasonalDoubleExpSmooth(X3, alpha, beta, delta, period)

plot(X3, type='l', xlab="Time", ylab=expression(X^3))
lines(X3.seassmooth$smooth, col='red')

par(mfrow=c(3,1))
plot(X3.seassmooth$b, type='l', xlab="Time", ylab=expression(b[t]), col='red')
plot(X3.seassmooth$l, type='l', xlab="Time", ylab=expression(l[t]), col='blue')
plot(X3.seassmooth$s, type='l', xlab="Time", ylab=expression(s[t]))
par(mfrow=c(1,1))

#### Forecast: X3

# We forecast at time t the time series at times t+1,...,t+20.
# In practice, we could refine our forecast of time t+20 at time t+19 for example.
# Here we wish to observe the deterioration for mid-term forecasts. 

## Exponential

last_obs = 80
horizon = 20

plot(1:last_obs, X3[1:last_obs], type="l", xlab="Time", ylab=expression(X^3),
     xlim=c(1,n), ylim=range(X3))
lines((last_obs+1):(last_obs+horizon), X3[(last_obs+1):(last_obs+horizon)], col="red")
legend("topleft", legend=c("Training data", "Forecasting data"), lty=1,
       col=c("black", "red"))

alpha = seq(0.05,0.95,length=100)

smoothing = lapply(alpha, ExpSmooth, x=X3[1:last_obs])

erreur = sapply(smoothing, mse_smoothing_simple, true_values=X3[1:last_obs])

smoothing = ExpSmooth(X3[1:last_obs], alpha[which.min(erreur)])
forecast = c(smoothing, rep(smoothing[last_obs], horizon))

plot(X3, pch=20, ylim=range(X3, forecast), xlab="Time", ylab=expression(X^3))
lines(forecast, col='red', lwd=2)
abline(v=last_obs, lty='dashed')
legend("topleft", legend=c("Observed data", "Smoothed and forecast"), lty=c(NA,1), pch=c(20,NA),
       lwd=c(1,2), col=c("black", "red"))

## Double exponential

smoothing = lapply(alpha, DoubleExpSmooth, x=X3[1:last_obs])

erreur = sapply(smoothing, mse_smoothing_methods, true_values=X3[1:last_obs])

smoothing = DoubleExpSmooth(X3[1:last_obs], alpha[which.min(erreur)])
forecast = c(smoothing$smooth, 
             smoothing$l[last_obs] + smoothing$b[last_obs]*(1:horizon))

plot(X3, pch=20, ylim=range(X3, forecast), xlab="Time", ylab=expression(X^3))
lines(forecast, col='red', lwd=2)
abline(v=last_obs, lty='dashed')
legend("topleft", legend=c("Observed data", "Smoothed and forecast"), lty=c(NA,1), pch=c(20,NA),
       lwd=c(1,2), col=c("black", "red"))

## Seasonal double exponential

alpha = 0.35
beta = 0.2
delta = 0.2
period = 10

smoothing = SeasonalDoubleExpSmooth(X3[1:last_obs], alpha, beta, delta, period)
forecast = c(smoothing$smooth, 
             smoothing$l[last_obs] + smoothing$b[last_obs]*(1:horizon) +
             smoothing$s[last_obs-period+1+(1:horizon-1)%%period])

plot(X3, pch=20, ylim=range(X3, forecast), xlab="Time", ylab=expression(X^3))
lines(forecast, col='red', lwd=2)
abline(v=last_obs, lty='dashed')
legend("topleft", legend=c("Observed data", "Smoothed and forecast"), lty=c(NA,1), pch=c(20,NA),
       lwd=c(1,2), col=c("black", "red"))


#####################
#### Forecast package 
#####################

library(forecast)

?ets

ets1 = ets(y=X1, model='ANN')
ets1
ets1$initstate
plot(ets1)

plot(X1, type='l')
lines(ets1$fitted, col='red')

ets1.forecast = array(NA, dim=n-1)
for(i in 1:(n-1)){
  ets1loop = ets(y=X1[1:i], model='ANN')
  ets1.forecast[i] = forecast(ets1loop)$mean[1]
  plot(forecast(ets1loop))
  lines(ets1loop$fitted, col='orange')
  lines(head(ets1$fitted,i), col='red')
}

plot(tail(X1, n-1), type='l')
lines(head(ets1$fitted, n-1), col='red')
lines(ets1.forecast, col='purple')
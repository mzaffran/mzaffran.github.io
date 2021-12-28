rm(list=objects())

setwd("/Users/these/Documents/TDs/STA202/TP6")

# We redefine the functions for the p-values

pvalue = function(model)
{
  (1-pnorm( abs(model$coef) / sqrt(diag(model$var.coef))) )*2
}

pvalue_BP = function(model, K)
{
  rho = acf(model$residuals, lag.max=K, plot=F)$acf[-1]
  n = model$nobs
  pval = (1-pchisq(n*sum(rho^2), df=K-length(model$coef)))
  return(pval)
}

data = read.table(file='TP6_exercice2.txt', header=T, sep=';')

View(data)
# 3 time series
attach(data)

#### X1

plot(x1, type='l', xlab="Time", ylab=expression(X[1]))
# We see a seasonality (small amplitude)

par(mfrow=c(1,2))
acf(x1,lag.max=60) 
pacf(x1,lag.max=60)
par(mfrow=c(1,1))
# Order 12 seasonality
# Quick decrease towards 0, no need for differentiation

s = 12

par(mfrow=c(1,2))
acf(x1, lag.max=3*s)        # qmax = 1, Qmax = 2
pacf(x1, lag.max=3*s)       # pmax = 2, Pmax = 2
par(mfrow=c(1,1))
# Qmax and Pmax = significant repetitions on the ACF/PACF
# qmax and pmax, as usual, on the first cycle
qmax = 1
Qmax = 2
pmax = 2
Pmax = 2

par(mfrow=c(1,2))
acf(x1, lag.max=3*s) 
abline(v=qmax+1, lty=3, col="red")
pacf(x1, lag.max=3*s) 
abline(v=pmax+1, lty=3, col="red")
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:pmax), q=c(0:qmax), P=c(0:Pmax), Q=c(0:Qmax))
ordre = cbind(ordre[,1], 0, ordre[,2], ordre[,3], 0, ordre[,4])

# We need to build a function because the seasonal and "principal" orders are not
# given in the same parameter to the arima function
sarima = function(x, ordre, s)
{
  arima(x, order=ordre[1:3],
        seasonal=list(order=ordre[4:6], period=s),
        include.mean=F)
}

model.sarima = apply(ordre, 1, sarima, x=x1, s=12)

aic = sapply(model.sarima, function(x) {x$aic})
bic = sapply(model.sarima, function(x) {-2*x$loglik+x$nobs*length(x$coef)})
like = sapply(model.sarima, function(x) {-2*x$loglik})

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, ylim=range(aic,like),
     xlab="(p,q,P,Q) SARIMA order", ylab="Criterion value")
points(like[o], col='red', pch=20, type='b')
axis(1, c(1:length(aic)), 
     paste(ordre[o,1], ordre[o,3], ordre[o,4], ordre[o,6]), las=2)
axis(2)
legend("topleft", legend=c("AIC", "-LogLikelihood"), col=c("black", "red"),
       pch=20)

ordre[which.min(aic),]

model.sarima[[which.min(aic)]]$coef
pvalue(model.sarima[[which.min(aic)]]) 
# We can not reject the null hypothesis of phi2 (coefficient of X_t-2, i.e. ar2)
# p-value > 0.05

# Thus, we select the second best model!
model.sarima[[order(aic)[2]]]$coef
pvalue(model.sarima[[order(aic)[2]]]) 
# This time, all the coefficient pass the test

## Residuals study

model.opt = model.sarima[[order(aic)[2]]]

par(mfrow=c(1,2))
acf(model.opt$residuals)
pacf(model.opt$residuals)
par(mfrow=c(1,1))
# I think (personnally) that there is still a component missing and that we should 
# refine our model. Nevertheless, it is not done in this session.

qqnorm(model.opt$residuals)
qqline(model.opt$residuals, col = "steelblue", lwd = 2)

hist(model.opt$residuals, breaks=50, freq=F, xlab="Residuals",
     main="Histogram of the residuals")
x = seq(min(model.opt$residuals),
        max(model.opt$residuals),
        length=50)
lines(x, dnorm(x,mean(model.opt$residuals),model.opt$sigma2), col='red')
# Furthermore, the histogram does not seem to fit correctly, in my point of view.

pvalue_BP(model.opt, K=10) 
# We do not reject the non-correlation of the residuals, since the p-value = 0.48 > 0.05

#### X2

plot(x2, type='l', xlab="Time", ylab=expression(X[2]))
# Again, a seasonality appears

acf(x2, lag.max=30)
# Order 7 seasonality
# Quick decrease towards 0, no differentiation

s = 7

par(mfrow=c(1,2))
acf(x2, lag.max=3*s)  
# Qmax = 2 significant repetitions
# qmax = 3 (deep of the cycle) 
pacf(x2, lag.max=3*s) 
# Pmax = 1 only one significant repetition 
# pmax = 1 (only one significant lag in the cycle)
par(mfrow=c(1,1))

qmax = 3
Qmax = 2
pmax = 1
Pmax = 1

par(mfrow=c(1,2))
acf(x2, lag.max=3*s) 
abline(v=qmax+1, lty=3, col="red")
pacf(x2, lag.max=3*s) 
abline(v=pmax+1, lty=3, col="red")
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:pmax), q=c(0:qmax), P=c(0:Pmax), Q=c(0:Qmax))
ordre = cbind(ordre[,1], 0, ordre[,2], ordre[,3], 0, ordre[,4])

model.sarima = apply(ordre, 1, sarima, x=x2, s=7)
aic = sapply(model.sarima, function(x) x$aic)
bic = sapply(model.sarima, function(x) -2*x$loglik+x$nobs*length(x$coef))
like = sapply(model.sarima, function(x) -2*x$loglik)

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, ylim=range(aic,like),
     xlab="(p,q,P,Q) SARIMA order", ylab="Criterion value")
points(like[o], col='red', pch=20, type='b')
axis(1, c(1:length(aic)), 
     paste(ordre[o,1], ordre[o,3], ordre[o,4], ordre[o,6]), las=2)
axis(2)
legend("topleft", legend=c("AIC", "-LogLikelihood"), col=c("black", "red"),
       pch=20)

ordre[which.min(aic),]

model.sarima[[which.min(aic)]]$coef
pvalue(model.sarima[[which.min(aic)]]) 
# Both coefficients pass the test

## Residuals study

model.opt = model.sarima[[which.min(aic)]]

par(mfrow=c(1,2))
acf(model.opt$residuals)
pacf(model.opt$residuals)
par(mfrow=c(1,1))

qqnorm(model.opt$residuals)
qqline(model.opt$residuals, col = "steelblue", lwd = 2)

hist(model.opt$residuals, breaks=50, freq=F, xlab="Residuals",
     main="Histogram of the residuals")
x = seq(min(model.opt$residuals),
        max(model.opt$residuals),
        length=50)
lines(x, dnorm(x,mean(model.opt$residuals),model.opt$sigma2), col='red')

pvalue_BP(model.opt, K=10) 
# We do not reject the non-correlation of the residuals, since the p-value = 0.78 > 0.05

#### X3

plot(x3,type='l', xlab="Time", ylab=expression(X[3]))
# need to differentiate, no stationary
acf(x3) 
# confirmed by the ACF

x3.diff = diff(x3, lag=1, differences=1)
plot(x3.diff,type='l', xlab="Time", ylab=expression(X[3] ~ differentiated))
acf(x3.diff)  
# Quick decrease towards 0, good
# No seasonality

par(mfrow=c(1,2))
acf(x3.diff, lag.max=20) # qmax = 4
pacf(x3.diff, lag.max=20) # pmax = 5
par(mfrow=c(1,1))

qmax = 4
pmax = 5

par(mfrow=c(1,2))
acf(x3.diff) 
abline(v=qmax+1, lty=3, col="red")
pacf(x3.diff) 
abline(v=pmax+1, lty=3, col="red")
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:pmax), q=c(0:qmax))

ordre = cbind(ordre[,1], 1, ordre[,2]) 
# this time, we add the 1 in order to differentiate: we will fit a complete ARIMA, 
# that will take care of the differentiation

model = apply(ordre, 1, arima, x=x3, method=c("ML"),
             SSinit=c("Rossignol2011"),
             optim.method="BFGS", include.mean=F)

aic = sapply(model, function(x) x$aic)

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, ylim=range(aic,like),
     xlab="(p,q) ARMA order", ylab="Criterion value")
points(like[o], col='red', pch=20, type='b')
axis(1, c(1:length(aic)), paste(ordre[o,1], ordre[o,3]), las=2)
axis(2)
legend("topleft", legend=c("AIC", "-LogLikelihood"), col=c("black", "red"),
       pch=20)

ordre[which.min(aic),]

model[[which.min(aic)]]$coef
pvalue(model[[which.min(aic)]])
# Both coefficients pass the test

## Residuals study

model.opt = model[[which.min(aic)]]

par(mfrow=c(1,2))
acf(model.opt$residuals)
pacf(model.opt$residuals)
par(mfrow=c(1,1))

qqnorm(model.opt$residuals)
qqline(model.opt$residuals, col = "steelblue", lwd = 2)

hist(model.opt$residuals, breaks=50, freq=F, xlab="Residuals",
     main="Histogram of the residuals")
x = seq(min(model.opt$residuals),
        max(model.opt$residuals),
        length=50)
lines(x, dnorm(x,mean(model.opt$residuals),model.opt$sigma2), col='red')

pvalue_BP(model.opt, K=10) 
# We do not reject the non-correlation of the residuals, since the p-value = 0.88 > 0.05

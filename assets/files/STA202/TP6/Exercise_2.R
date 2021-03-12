rm(list=objects())

setwd("/Users/these/Documents/TDs/STA202/TP6")

pvalue = function(model)
{
  (1-pnorm( abs(model$coef) / sqrt(diag(model$var.coef))) )*2
}

pvalue_BP = function(model, K)
{
  rho = acf(model$residuals, lag.max=K, plot=F)$acf[-1] # we delete the first term, always = 1
  n = model$nobs
  pval = (1-pchisq(n*sum(rho^2), df=K-length(model$coef)))
  return(pval)
}

data = read.table(file='TP6_exercice2.txt', header=T, sep=';')

View(data)
# 3 time series
attach(data)

#### X1

plot(x1, type='l')

par(mfrow=c(1,2))
acf(x1,lag.max=60) 
pacf(x1,lag.max=60)
par(mfrow=c(1,1))
# Order 12 seasonality
# Quick decrease towards 0, no differenciation

s = 12

par(mfrow=c(1,2))
acf(x1, lag.max=3*s)        # qmax = 1, Qmax = 2
pacf(x1, lag.max=3*s)       # pmax = 2, Pmax = 2
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:2), q=c(0:1), P=c(0:2), Q=c(0:2))
ordre = cbind(ordre[,1], 0, ordre[,2], ordre[,3], 0, ordre[,4])
dim(ordre)

sarima = function(x, ordre, s)
{
  arima(x, order=ordre[1:3],
        seasonal=list(order=ordre[4:6], period=s),
        include.mean=F)
}

model.sarima = apply(ordre, 1, sarima, x=x1, s=12)
aic = sapply(model.sarima,function(x) {x$aic})
bic = sapply(model.sarima,function(x) {-2*x$loglik+x$nobs*length(x$coef)})
like = sapply(model.sarima,function(x) {-2*x$loglik})

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, xlab='')
axis(1, c(1:length(aic)), 
     paste(ordre[o,1], ordre[o,3], ordre[o,4], ordre[o,6]),
     las=2)
axis(2)

ordre[which.min(aic),]
ordre[which.min(bic),]

model.sarima[[which.min(aic)]]$coef
pvalue(model.sarima[[which.min(aic)]]) 
#####attention phi2 non signi. non-nulle ? 5% on regarde le mod?le suivant en AIC

model.sarima[[order(aic)[2]]]$coef
pvalue(model.sarima[[order(aic)[2]]]) # correct!

# Residuals study

pvalue_BP(model.sarima[[order(aic)[2]]], K=10) #ok

par(mfrow=c(1,2))
acf(model.sarima[[order(aic)[2]]]$residuals)
pacf(model.sarima[[order(aic)[2]]]$residuals)
par(mfrow=c(1,1))

#### X2

plot(x2, type='l')

acf(x2, lag.max=30)
# Order 7 seasonality
# Quick decrease towards 0, no differenciation

s=7

par(mfrow=c(1,2))
acf(x2, lag.max=3*s)  # qmax = 3, Qmax = 2
pacf(x2, lag.max=3*s) # pmax = 1, Pmax = 1
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:1), q=c(0:3), P=c(0:1), Q=c(0:2))
ordre = cbind(ordre[,1], 0, ordre[,2], ordre[,3], 0, ordre[,4])
dim(ordre)

model.sarima = apply(ordre, 1, sarima, x=x2, s=7)
aic = sapply(model.sarima,function(x) x$aic)
bic = sapply(model.sarima,function(x) -2*x$loglik+x$nobs*length(x$coef))
like = sapply(model.sarima,function(x) -2*x$loglik)

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, xlab='')
axis(1, c(1:length(aic)),
     paste(ordre[o,1], ordre[o,3], ordre[o,4], ordre[o,6]),
     las=2)
axis(2)

ordre[which.min(aic),]

model.sarima[[which.min(aic)]]$coef
pvalue(model.sarima[[which.min(aic)]]) 

# Residuals study

pvalue_BP(model.sarima[[which.min(aic)]], K=10) #ok

par(mfrow=c(2,1))
acf(model.sarima[[which.min(aic)]]$residuals)
pacf(model.sarima[[which.min(aic)]]$residuals)
par(mfrow=c(1,1))

#### X3

plot(x3,type='l')

acf(x3) # need to differenciate, no stationary

x3.diff = diff(x3, lag=1, differences=1)
acf(x3.diff)  
# Quick decrease towards 0, good
# No seasonality
plot(x3.diff,type='l')

par(mfrow=c(1,2))
acf(x3.diff, lag.max=20)        # qmax = 4
pacf(x3.diff, lag.max=20)       # pmax = 5
par(mfrow=c(1,1))

ordre = expand.grid(p=c(0:5), q=c(0:4))
dim(ordre)

ordre = cbind(ordre[,1], 1, ordre[,2])

model = apply(ordre, 1, arima, x=x3, method=c("ML"),
             SSinit=c("Rossignol2011"),
             optim.method="BFGS", include.mean=F)

aic = sapply(model, function(x) x$aic)

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, xlab='')
axis(1, c(1:length(aic)), paste(ordre[o,1],ordre[o,3]), las=2)
axis(2)

ordre[which.min(aic),]

model[[which.min(aic)]]$coef
pvalue(model[[which.min(aic)]])

acf(model[[which.min(aic)]]$residuals)
qqnorm(model[[which.min(aic)]]$residuals)
qqline(model[[which.min(aic)]]$residuals, col = "steelblue", lwd = 2)

# portemanteau pierce test

K = 10
pvalue_BP(model[[which.min(aic)]], K)

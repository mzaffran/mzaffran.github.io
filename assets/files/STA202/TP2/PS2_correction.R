rm(list=objects())
library(xts)
library(mgcv)

###### Time series generation

## Date creation

date1 = strptime("01/01/1900", "%m/%d/%Y")
date2 = strptime("01/01/2000", "%m/%d/%Y")
Date2 = seq(date1,date2, by = "year")

n = length(Date)
t = 1:n # temporal index

## Data simulation

Trend = t/20 + 1

w = 2*pi/50
S = cos(w*t) + sin(w*t)

eps = rnorm(n,0,1)

X = Trend + S + eps

X = xts(X, order.by=Date)
Trend = xts(Trend, order.by=Date)
S = xts(S, order.by=Date)

## Data visualization

plot(X)
lines(Trend, col='red')
lines(Trend + S, col='blue')
addLegend("topleft", legend.names = c(expression(X[t]), expression(S[t]+T[t]), expression(T[t])),
          lty=1, col=c("black", "blue", "red"), bty="o")

acf(X, lag.max=50)

###### Trend estimation

## Linear regression

reg = lm(X~t)
summary(reg)
ychap.lm = reg$fitted
plot(X, type='l', main = "Trend estimation by linear regression")
lines(ychap.lm, col='red')
lines(Trend, col='blue')
addLegend("topleft", legend.names = c(expression(T[t]), expression(widehat(T[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Moving average

mb = filter(X, filter = array(1/50, dim=50), method = c("convolution"), sides = 2, circular = F)
mb = xts(mb,order.by=Date)

plot(X, type='l', main = "Trend estimation by moving average")
lines(mb, col='red')
lines(Trend, col='blue')
addLegend("topleft", legend.names = c(expression(T[t]), expression(widehat(T[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Gaussian kernel

h = 20
x = 1:n

noyau = function(x){
          return( dnorm(x-t, 0, sd=h) / sum(dnorm(x-t, 0, sd=h)))
        }

W = matrix(unlist(lapply(x,noyau)), ncol=n, nrow=n, byrow=F) # each column corresponds to one x

plot(W[,50], type='l')

ychap.kernel = colSums(as.numeric(X)*W)
ychap.kernel = xts(ychap.kernel, order.by=Date)
plot(X, type='l', main = "Trend estimation by Gaussian kernel")
lines(ychap.kernel, col='red')
lines(Trend, col='blue')
addLegend("topleft", legend.names = c(expression(T[t]), expression(widehat(T[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Local polynomials

lo = loess(X~t, degree=1, span=0.9) 
# span controls the smoothing: if too small (close to 0) we estimate the "whole curve"
ychap.lo = xts(lo$fitted,order.by=Date)
plot(X, type='l', main = "Trend estimation by local polynomials")
lines(ychap.lo, col='red')
lines(Trend, col='blue')
addLegend("topleft", legend.names = c(expression(T[t]), expression(widehat(T[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Spline basis

g = gam(X~s(t, k=3))
summary(g)
ychap.gam = xts(g$fitted,order.by=Date)
plot(X, type='l', main = "Trend estimation by spline basis")
lines(ychap.gam, col='red')
lines(Trend, col='blue')
addLegend("topleft", legend.names = c(expression(T[t]), expression(widehat(T[t]))),
          lty=1, col=c("blue", "red"), bty="o")

###### Seasonal estimation

## Removing the trend: linear regression

X.detrend = X-ychap.lm
plot(X.detrend)

## Fourier serie

w = 2*pi/50
fourier = cbind(cos(w*t), sin(w*t))
#K = 20
#for(i in 2:K)
#  {
#    fourier = cbind(fourier,cos(i*w*t), sin(i*w*t))
#   }
matplot(fourier, type='l')
dim(fourier)

reg = lm(X.detrend~fourier[,1:2])
ychap.lm.season = xts(as.numeric(reg$fitted), order.by=Date)
plot(X.detrend, type='l',  main = "Seasonal estimation by Fourier")
lines(ychap.lm.season, col='red')
lines(S, col='blue')
addLegend("bottomleft", legend.names = c(expression(S[t]), expression(widehat(S[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Moving average

K = 20
mb.season = filter(X.detrend, filter=array(1/K,dim=K), method = c("convolution"), 
                  sides = 2, circular = TRUE)
mb.season = xts(mb.season,order.by=Date)

plot(X.detrend, type='l',  main = "Seasonal estimation by moving average")
lines(mb.season, col='red')
lines(S, col='blue')
addLegend("bottomleft", legend.names = c(expression(S[t]), expression(widehat(S[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Gaussian kernel

h = 5
x = 1:n

W = matrix(unlist(lapply(x,noyau)),ncol=n,nrow=n,byrow=F)
plot(W[,10])
ychap.kernel.season = colSums(as.numeric(X.detrend)*W)
ychap.kernel.season = xts(ychap.kernel.season,order.by=Date)

plot(X.detrend, type='l',  main = "Seasonal estimation by Gaussian kernel")
lines(ychap.kernel.season, col='red')
lines(S, col='blue')
addLegend("bottomleft", legend.names = c(expression(S[t]), expression(widehat(S[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Local polynomials

lo = loess(X.detrend~t, degree=1, span=0.3)
ychap.lo.season = xts(lo$fitted, order.by=Date)

plot(X.detrend, type='l',  main = "Seasonal estimation by local polynomials")
lines(ychap.lo.season, col='red')
lines(S, col='blue')
addLegend("bottomleft", legend.names = c(expression(S[t]), expression(widehat(S[t]))),
          lty=1, col=c("blue", "red"), bty="o")

## Cyclic splines regression

cycle = c(rep(c(1:50),2),1)
plot(cycle)

plot(cycle, X.detrend, pch=20)

g = gam(X.detrend~s(cycle, k=20, bs='cc'))
summary(g)
ychap.gam.season = xts(g$fitted,order.by=Date)
plot(X.detrend, type='l',  main = "Seasonal estimation by cyclic splines")
lines(ychap.gam.season, col='red')
lines(S, col='blue')
addLegend("bottomleft", legend.names = c(expression(S[t]), expression(widehat(S[t]))),
          lty=1, col=c("blue", "red"), bty="o")


#### Methods comparison

library(RColorBrewer)
cols = brewer.pal(5, "Set1")

plot(X-eps, type='l', ylim=range(X))
lines(X, col='grey')
lines(ychap.lm+ychap.lm.season, col=cols[1])
lines(ychap.lm+ychap.kernel.season, col=cols[2])
lines(ychap.lm+ychap.lo.season, col=cols[3])
lines(ychap.lm+ychap.gam.season, col=cols[4])
lines(ychap.lm+mb.season, col=cols[5])
addLegend("topleft", legend.names = c("Fourier", "Gaussian kernel", "Local polynomials",
                                      "Cyclic splines", "Moving average"),
          lty=1, col=cols, bty="o")

# Root-Mean-Squared-Error
rmse = function(u){
  return(sqrt(mean(u^2)))
}

# Mean-Absolute-Error
mae = function(u){
  return(mean(abs(u)))
}

errors = matrix(nrow = 5, ncol = 2)
colnames(errors) = c("RMSE","MAE")

epschap = ychap.lm+ychap.lm.season
acf(X-epschap, main = "ACF of the residuals of the Fourier fit")
qqnorm(epschap, pch = 1, frame = FALSE, main = "Q-Q plot of the residuals of the Fourier fit")
qqline(epschap, col = "steelblue", lwd = 2)
errors[1,"RMSE"] = rmse(epschap)
errors[1,"MAE"] = mae(epschap)

epschap = ychap.lm+ychap.kernel.season
acf(X-epschap, main = "ACF of the residuals of the Gaussian kernel fit")
qqnorm(epschap, pch = 1, frame = FALSE, main = "Q-Q plot of the residuals of the Gaussian kernel fit")
qqline(epschap, col = "steelblue", lwd = 2)
errors[2,"RMSE"] = rmse(epschap)
errors[2,"MAE"] = mae(epschap)

epschap = ychap.lm+ychap.lo.season
acf(X-epschap, main = "ACF of the residuals of the local polynomials fit")
qqnorm(epschap, pch = 1, frame = FALSE, main = "Q-Q plot of the residuals of the local polynomials fit")
qqline(epschap, col = "steelblue", lwd = 2)
errors[3,"RMSE"] = rmse(epschap)
errors[3,"MAE"] = mae(epschap)

epschap = ychap.lm+ychap.gam.season
acf(X-epschap, main = "ACF of the residuals of the cyclic splines fit")
qqnorm(epschap, pch = 1, frame = FALSE, main = "Q-Q plot of the residuals of the cyclic splines fit")
qqline(epschap, col = "steelblue", lwd = 2)
errors[4,"RMSE"] = rmse(epschap)
errors[4,"MAE"] = mae(epschap)

epschap = ychap.lm+mb.season
acf(X-epschap,na.action = na.omit, main = "ACF of the residuals of the moving average fit")
qqnorm(epschap, pch = 1, frame = FALSE, main = "Q-Q plot of the residuals of the moving average fit")
qqline(epschap, col = "steelblue", lwd = 2)
errors[5,"RMSE"] = rmse(epschap)
errors[5,"MAE"] = mae(epschap)

View(errors)

# ==> They seem quite equivalent, except maybe the local polynomials approach.
# Nevertheless, the Q-Q plot of the Fourier approach seem nicer.

#### Application to beer data-set

setwd("/Users/these/Documents/TDs/STA202/TP1/TP1_data")

beer = read.csv("beer2.csv", header=TRUE, skip=1)

## Date creation

date1 = strptime(c("01/01/91"), "%m/%d/%y")
date2 = strptime(c("08/01/95"), "%m/%d/%y")
Date = seq(date1, date2, by = "month")
Time = c(1:length(Date))
beer = data.frame(Date,beer$BeerProd,Time)
names(beer) = c("Date","BeerProd","Time")
plot(beer$Date, beer$BeerProd, type='l', xlab="Date", ylab="Beer production")

## Spline regression

Month = as.numeric(format(Date,"%m"))
beer = data.frame(beer, Month)

g = gam(BeerProd~s(Time,k=10)+s(Month,k=4,bs='cc'), data=beer)
ychap.gam = g$fitted

plot(beer$Date, beer$BeerProd, type='l')
lines(beer$Date, ychap.gam, col='red')

plot(g)

terms = predict(g, newdata=beer, type="terms")

plot(beer$Date, beer$BeerProd-mean(beer$BeerProd), type='l')
lines(beer$Date, terms[,1], col='blue')
lines(beer$Date, terms[,2], col='red')






































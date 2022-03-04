rm(list=objects())
setwd("~/Documents/TDs/STA202/TP5/TP5_data")

data = read.table('exercice3.txt',header=T,sep=';')

View(data)

attach(data)

horizon = 10

#### x1

x1 = ts(x1)
par(mfrow=c(1,2))
acf(x1)
pacf(x1) # AR order 2
par(mfrow=c(1,1))
mean(x1)

x1.model = arima(x1, order=c(2,0,0), method=c("ML"),
                 SSinit=c("Rossignol2011"),
                 optim.method="BFGS", include.mean=F)
# in "arima", order = (p,d,q) with p the order of the AR component, q of the MA component
# and d the degree of differencing

x1.model

x1.forecast = predict(x1.model, n.ahead=horizon, se.fit=F)
plot(x1, xlim=c(1,nrow(data)+horizon))
lines(nrow(data)+c(1:horizon), x1.forecast, col='red')

names(x1.model)

#### x2

x2 = ts(x2)
par(mfrow=c(1,2))
acf(x2)
pacf(x2) # AR order 1
par(mfrow=c(1,1))

x2.model = arima(x2, order=c(1,0,0), method=c("ML"), SSinit=c("Rossignol2011"),
                 optim.method="BFGS", include.mean=F)

x2.forecast = predict(x2.model, n.ahead=horizon, se.fit=F)
plot(x2, xlim=c(1,nrow(data)+horizon))
lines(nrow(data)+c(1:horizon), x2.forecast, col='red')


rho = ARMAtoMA(ar=x2.model$coef, ma=0, 12) 
rho
0.7^c(1:12)
plot(rho,type='l')
lines(0.7^c(1:12), col='red')

#### x3

x3 = ts(x3)
par(mfrow=c(1,2))
acf(x3) # MA order 6, no exponential decrease on the PACF! Thus, not an AR.
pacf(x3)
par(mfrow=c(1,1))
mean(x3)
sd(x3)

x3.model = arima(x3, order=c(0,0,6), method=c("ML"), SSinit=c("Rossignol2011"),
                optim.method="BFGS", include.mean=T)
x3.model

x3.forecast = predict(x3.model, n.ahead=horizon, se.fit=F)

plot(x3, xlim=c(1,nrow(data)+horizon))
lines(nrow(data)+c(1:horizon), x3.forecast, col='red')

#### x4

x4 = ts(x4)
par(mfrow=c(1,2))
acf(x4,lag.max=100)
pacf(x4,lag.max=100)
par(mfrow=c(1,1))
mean(x4)
x4.model = arima(x4, order=c(0,0,20), method=c("ML"),SSinit=c("Rossignol2011"),
                optim.method="BFGS", include.mean=F)
x4.model

x4.forecast = predict(x4.model, n.ahead=horizon, se.fit=F)

plot(x4, xlim=c(1,nrow(data)+horizon))
lines(nrow(data)+c(1:horizon), x4.forecast, col='red')

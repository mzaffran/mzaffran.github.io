rm(list=objects())
setwd("~/Documents/TDs/STA202/TP5/TP5_data")

#### Data import

data = read.table("exercice2.txt",sep=';',header=T)
View(data)

plot(data$y, type='l', xlab="Time", ylab=expression(Y[t]))

#### Spectrum

s = spectrum(data$y)
plot(s$freq, s$spec, type='l')

a1 = which.max(s$spec)
f1 = s$freq[a1]
a2 = which.max(s$spec[-a1])
f2 = s$freq[a2]
a3 = which.max(s$spec[-c(a1,a2)])
f3 = s$freq[a3]

X1cos = cos(2*pi*f1*data$t)
X1sin = sin(2*pi*f1*data$t)
X2cos = cos(2*pi*f2*data$t)
X2sin = sin(2*pi*f2*data$t)
X3cos = cos(2*pi*f3*data$t)
X3sin = sin(2*pi*f3*data$t)

model = lm(data$y ~ X1cos + X1sin + X2cos + X2sin + X3cos + X3sin - 1) 
# -1 so that we do not estimate an intercept
# by default, lm fits an intercept
summary(model)

plot(data$y, type='l', xlab="Time", ylab=expression(Y[t]))
lines(model$fitted.values, col="red")

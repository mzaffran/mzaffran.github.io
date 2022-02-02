rm(list=objects())

setwd("~/Documents/TDs/STA202/TP1/TP1_data")

########## Exercise 2: electricity consumption

###### Data import

data = read.csv("conso_2015.csv", header=T, sep=';')

###### First visualization

summary(data)
head(data) 

# Each column corresponds to a different time (steps of 30 minutes).
# Not practical, we prefer to have each row corresponding to a precise date and time.

###### Time serie creation

## Date creation

# head displays the first 6 elements (by defaut, set n=k to obtain the k first elements)
head(data$Date) # Begins the 01/01/2015
# tail is the same than head, for the last elements
tail(data$Date) # Ends the 30/11/2015

date1 = strptime("01/01/2015 00:30:00", "%m/%d/%Y %H:%M:%S")
date2 = strptime("11/30/2015 24:00:00", "%m/%d/%Y %H:%M:%S")
Date = seq(date1, date2, by = "30 min")
summary(Date)
head(Date)

## Time serie

X = as.matrix(t(data[,-1])) 
# All column except the first one (the date), 
# and transposed: each column now corresponds to one day
conso = c(X)
# we transform X to a vector
# conso is the concatenation (by columns) of the matrix X: 
# each column has been put after the previous one.

plot(Date[1:(48*7*3)],conso[1:(48*7*3)], type='l', ylab="Electricity consumption", xlab='Date') # Plots only the 3 first weeks of consumption
plot(Date, conso, type='l', ylab="Electricity consumption") # Plots all the data

library(xts)

conso.xts = xts(conso, order.by=Date)
# note again here that we use the Date vector previously created
plot(conso.xts)
# xts makes automatically beautiful plots!

###### Descriptive analysis

## Statistics

mean(conso.xts)

?.indexmon

## Month average

month = as.factor(.indexmon(conso.xts))
mean.month = tapply(conso.xts, month, mean)
# we could also use the code of Exercise 1, 
# notice it is the same with another way of creating the factor.
plot(mean.month, type='b', xlab="Month", ylab="Electricity consumption", 
     main="Averaged electricity consumption per month")

## Profiles

# Hourly profile

dow = as.factor(.indexwday(conso.xts)) # dow = day of the week
conso_day = tapply(conso.xts, dow, mean)
plot(conso_day, type='b', xlab="Day of the week", ylab="Electricity consumption", 
     main="Averaged electricity consumption per type of day")

# Hourly profile for each day

hour = as.factor(.indexhour(conso.xts))
mean.dow.hour = tapply(conso.xts, dow:hour, mean)
plot(mean.dow.hour, type='l', ylab="Electricity consumption", 
     main="Hourly profile per type of day")
abline(v=seq(1, 24*7, by=24), col='red')

## Month boxplots at 8PM

col.pal = colorRampPalette(c("lightblue", "red"))( 12 )
sel = which(.indexhour(conso.xts)==20) # indexes of data at 8PM
boxplot(conso[sel] ~ month[sel], col=col.pal, xlab="Months", ylab="Electricity consumption")
# In the function boxplot, you can specify a formula "A ~ B" to show 
# the boxplots of A conditionnally on the factor B.
# Here, we also specify [sel] to select only the rows of indexes in sel,
# where sel was created before to contain only row indexes of hour 20.
# We only had 11 months: only one data in December in the data set.

###### Autocorrelations

## ACF

# Lag function
lag.test = lag.xts(conso.xts, k=1, na.pad=T)
lag.test[1:3]
conso.xts[1:3]

# Autocorrelation of order h

autoCorr = function(x,h)
{
  x.lag = lag.xts(x, k=h, na.pad=T)
  return(cor(x.lag, x, use="pairwise.complete.obs"))
}

autoCorr2 = function(h,x)
{
  n = length(x)
  x_mean = mean(x)
  autocov = (1/(n-h))*sum(as.numeric(x[1:(n-h)] - x_mean)*as.numeric(x[(h+1):n] - x_mean))
  return(autocov/((1/n)*sum(as.numeric(x - x_mean)^2)))
}

a1 = sapply(c(1:336), autoCorr, x=conso.xts)
a2 = sapply(c(1:336), autoCorr2, x=conso.xts)

plot(a1,type='h',ylim=c(0,1))
lines(a2,col='red')

a3 = acf(conso.xts, lag.max=336, type="correlation")

plot(a3$acf[-1], type='h', ylim=c(0,1))
lines(a2, col='red')

# Difference because the ACF function uses the mean over the n-h observations, and not the whole mean.

## PACF

PartialAutoCorr<-function(x,h)
{
  x.lag = lapply(c(1:h), lag.xts, x=x, na.pad=T)
  x.lag = matrix(unlist(x.lag), ncol=length(x.lag))
  reg = lm(x~x.lag-1)
  return(as.numeric(tail(reg$coef,1)))
}

PartialAutoCorr(conso, h=1)
pa1 = sapply(c(1:48), PartialAutoCorr, x=conso.xts)
plot(pa1, type='b')

pa2 = pacf(conso, lag.max=7*48*4) 
points(pa2$acf)

pa1-pa2$acf[1:48]

rm(list=objects())

setwd("~/Documents/TDs/STA202/TP1/TP1_data")

########## Exercise 1: beer production


###### Data import

beer = read.csv("beer2.csv",header=TRUE, skip=1)
# to know that you have a header AND that you want to skip the first row,
# you should have first tried beer = read.csv("beer2.csv") 
# (or opened the .csv using your favorite sheet app)
# and by doing so you would have seen the title.

###### First visualization

head(beer)
str(beer)
summary(beer)
plot(beer$BeerProd, type='b', pch=1, ylab="Beer production")

####### Date creation

date1 = strptime(c("01/01/91"), "%m/%d/%y")
date2 = strptime(c("08/01/95"), "%m/%d/%y")

Date = seq(date1, date2, by = "1 month")

# Data frame are really practical!
beer = data.frame(Date, beer$BeerProd)
names(beer) = c("Date","BeerProd")

summary(beer)
plot(beer$Date, beer$BeerProd, type='l', xlab="Date", ylab="Beer production")

####### ts class
beer.ts = ts(beer$BeerProd, start=1, frequency=12) 
# each 12 observations, you have completed one period/season 
# (def of the frequency)
# here we assume a yearly seasonality
plot(beer.ts)

# change the marks of the x-axis:

plot(beer.ts, xaxt='n')
xtick = seq(1, 5.7, length=56) 
# 5.7 = 5 + 8/12, number of years (=seasons) in our data set
axis(1, xtick, labels=Date)

# or, proper version:

plot(beer.ts, xaxt='n')
t = time(beer.ts)
xtick = seq(t[1], tail(t, 1), length=56) 
axis(1, xtick, labels=Date)

####### zoo class
library(zoo)
beer.zoo = zoo(beer$BeerProd, order.by=beer$Date)
# note that here we reuse the Date we have created earlier!
plot(beer.zoo)

####### Base statistics

mean(beer$BeerProd)
sd(beer$BeerProd)
summary(beer)

boxplot(beer$BeerProd)

hist(beer$BeerProd, breaks=20, main="Histogram of the monthly Australian beer production",
     xlab="Beer production")
dim(beer)

year = format(beer$Date,"%Y")
mean.year = tapply(beer$BeerProd, as.factor(year), mean)
# tapply: apply a function (third argument) to a vector (first argument)
# grouped by a category (second argument)
# you need to convert your categories to a factor 
# (true type for R categories) using as.factor. 
plot(mean.year, type='b', axes=F, xlab="Date", ylab="Beer production",
     main="Monthly beer production averaged per year")
axis(1, c(1:5), names(mean.year))
axis(2)

####### Monthly average

month = format(beer$Date,"%m")
mean.month = tapply(beer$BeerProd,as.factor(month),mean)
plot(mean.month, type='b', axes=F, xlab="Date", ylab="Beer production",
     main="Monthly beer production averaged per month")
axis(1, c(1:12), names(mean.month))
axis(2)

rm(list=objects())

#### Simulation and representation

set.seed(110) # to ensure reproducible results (fix the seed of the randomness)

n = 1000
sigma = 1
eps = rnorm(n,0,sd=sigma) # innovation of our ARMA process

ar = c(1,-1/4) # AR orders
ma = c(-1) # MA orders
# Be careful: arima.sim reads the orders when the equation is written in the following way:
# X_t = phi_1*X_t-1 + phi_2*X_t-2 + ... + eps_t + theta_1*eps_t-1 + ...
x1 = arima.sim(n=n, model=list(ar=ar,ma=ma), innov=eps)
# arima.sim returns a TS object. 

plot(x1)

#### ACF and PACF

par(mfrow=c(1,2)) # 1 line, 2 columns
acf(x1)
pacf(x1)
par(mfrow=c(1,1)) # go back to only one plot per window

# Be careful: ACF starts at lag 0 (even though ACF at lag 0 is ALWAYS 1)
# and PACF at lags 1
qmax = 5 # ACF = 0 for h > 5
pmax = 6 
# PACF goes down the line at lag 7 and becomes higher only once after
# We prefer to have a smaller bound on p so that we do not try so many models afterwards
# If our models are not satisfactory, we could increase p at the end!

par(mfrow=c(1,2))
acf(x1)
abline(v=qmax+1, lty=3, col="red")
pacf(x1)
abline(v=pmax+1, lty=3, col="red")
par(mfrow=c(1,1))

#### Model selection

## Creation of the different orders without for loop

ordre = expand.grid(p=c(0:pmax), q=c(0:qmax)) 
View(ordre)
# create all pairs of potentials p and q
# matrix of size (length(q)*length(p)) x 2 = 42 x 2
ordre = cbind(ordre[,1], 0, ordre[,2]) 
# add the 0 order for differentiation in the middle, since we did not observe any trend

head(ordre)
dim(ordre)

## Fitting the different models without for loop

# arima function fits an ARIMA model to our data, given the order we gave it
model = apply(ordre, 1, arima, x=x1, method=c("ML"), SSinit=c("Rossignol2011"),
              optim.method="BFGS", include.mean=F) 
# 1 to specifiy that "apply" is by row (2 = by column)
# returns a list of 42 ARIMA models

# We can have a look at one model to understand the output of arima
model[[2]]

## Selection

aic = sapply(model, function(x) {x$aic}) # returns a list containing the 42 AIC
bic = sapply(model, function(x) {-2*x$loglik+log(x$nobs)*length(x$coef)}) # returns a list containing the 42 BIC
like = sapply(model, function(x) {-2*x$loglik}) # returns a list containing the 42 "- log likelihood"

# Understanding in details the following code is not important for the exam. 
# But interesting for those who will continue a bit with stats.
o = order(aic)
# order returns the indices corresponding to sorting the vector in ascending order
# it does not sort the vector! it gives the indices to sort it
# (if you just need to sort, you can use sort(aic) but here we will sort also
# BIC and loglikeklihood with respect to the AIC order)
# aic[o] is the sorted vector of AIC
# it allows us to sort also BIC and LogLikelihood in the AIC order
plot(aic[o], type='b', pch=20, axes=F, ylim=range(aic,like),
     xlab="(p,q) ARMA order", ylab="Criterion value")
points(like[o], col='red', pch=20, type='b')
points(bic[o], col='green', pch=20, type='b')
axis(1,c(1:length(aic)), paste(ordre[o,1], ordre[o,3]), las=2)
# we change the x axis so that the labels correspond to (p,q) and not the index in the list
axis(2)
legend("topleft", legend=c("AIC","BIC", "-LogLikelihood"), col=c("black", "green", "red"),
       pch=20)

# Chosen model by AIC
ordre.opt = ordre[which.min(aic),]
ordre.opt

# Chosen model by BIC: from the plot we know it is the same than AIC
ordre.opt = ordre[which.min(bic),]
ordre.opt

# we retrieve the model associated to the lowest AIC
model.opt = model[[which.min(aic)]] 
#model.opt = arima(x=x1, order=ordre.opt, method=c("ML"),
#                 SSinit=c("Rossignol2011"),
#                 optim.method="BFGS", include.mean=F)

model.opt
names(model.opt)

## Associated parameters

model.opt$coef # coefficients
model.opt$sigma2 # variance estimation of the innovation/noise
model.opt$var.coef # variance estimation of the estimated coefficients

#### Student test of the coefficients nullity

## pvalue of the Gaussian test (approximation of the Student test for n large)
pvalue = function(model)
{
  2*(1 - pnorm( abs(model$coef) / sqrt(diag(model$var.coef))) )
}

pvalue(model.opt)
# for each coefficient, the p-value is really below 0.05: we significantly reject 
# the null hypothesis, the hypothesis that the coefficient is in fact 0

## Higher orders

# We can test higher orders, to see if we are not missing a coefficient
# Be careful: we can not test both AR and MA higher orders, so we proceed in 2 times

model.opt_pp1 = arima(x=x1, order=ordre.opt+c(1,0,0), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_pp1)
# p-value for AR3 coefficient is 0.87: we can not reject the null hypothesis (that the 
# coefficient is null) at 5% (and even higher). It confirms the order 2 for the AR part.

model.opt_qp1 = arima(x=x1, order=ordre.opt+c(0,0,1), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_qp1)
# p-value for AR3 coefficient is 0.86: we can not reject the null hypothesis (that the 
# coefficient is null) at 5% (and even higher). It confirms the order 2 for the AR part.

## Lower orders

# Be careful: this is just an illustration for the practical session
# Our coefficient are significant according to the p-value, 
# so there is no need to look for lower orders

model.opt_pm1 = arima(x=x1, order=ordre.opt-c(1,0,0), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_pm1)

model.opt_qm1 = arima(x=x1, order=ordre.opt-c(0,0,1), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_qm1)

# Their coefficient pass the test, but it was predictable.

#### Residuals analysis

plot(model.opt$residuals, type='l', ylab="Residuals")
# Visually there is no trend and no seasonality.

par(mfrow=c(1,2))
acf(model.opt$residuals)
pacf(model.opt$residuals)
par(mfrow=c(1,1))
# No significant auto-correlations

# If the residuals are truly gaussian distributed, then their empirical quantiles
# should coincide with the theoretical guassian quantiles
qqnorm(model.opt$residuals)
qqline(model.opt$residuals, col = "steelblue", lwd = 2) # to add the x=y line
# Perfect fit of the quantiles

hist(model.opt$residuals, breaks=50, freq=F, xlab="Residuals",
     main="Histogram of the residuals")
x = seq(min(model.opt$residuals),
        max(model.opt$residuals),
        length=50)
lines(x, dnorm(x,mean(model.opt$residuals),model.opt$sigma2), col='red')
# Histogram globally coincides also with the one of the standard Gaussian distribution

#### Box-Pierce test

# In this test we do not want to reject the null hypothesis
# Indeed, the null hypothesis corresponds to the absence of residuals correlation
pvalue_BP = function(model, K)
{
  rho = acf(model$residuals, lag.max=K, plot=F)$acf[-1] # we delete the first term, always = 1
  n = model$nobs
  pval = (1-pchisq(n*sum(rho^2), df=K-length(model$coef)))
  return(pval)
}

pvalue_BP(model.opt, K=10)
pvalue_BP(model.opt, K=20)
# We do not reject H0, with a threshold of 5%
# => Residuals appear to be not correlated, according to this test
# More rigorously, we have no reason to be worried (not rejecting the test does
# not exactly mean that the hypothesis is correct, but that we do not have reason
# to suspect it is not)

# We can test ARMA(2,0) and ARMA(1,1) for the example

pvalue_BP(model.opt_pm1, K=10)
# We reject the hypothesis of non-correlated residuals at order at most 10
pvalue_BP(model.opt_qm1, K=10) 
# We reject the hypothesis of non-correlated residuals at order at most 10

# => missing significant higher orders leads to poor modelisation

#### Forecasting

x1_train = x1[1:900]
x1_test = x1[-c(1:900)]

# we fit the model only on the first 900 observations
model.opt = arima(x=x1_train, order=ordre.opt, method=c("ML"),
                  SSinit=c("Rossignol2011"),
                  optim.method="BFGS", include.mean=F)

# each point X_t is forecasted given the observations until time t-h
forecast = function(model, h)
{
  forecast = array(NA, dim=100)
  for(i in c(900:999)) # for each value to be predicted
  {
    # we apply the formulae with the fitted coefficients (of line 222)
    # but on the observed values! It is just an evaluation
    mod = arima(x=x1[1:(i-h)], order=ordre.opt, fixed=model$coef, include.mean=F)
    # i-899 because forecast is a vector of size 100, from 1 to 100
    forecast[i-899] = tail(predict(mod,n.ahead=h)$pred, 1) 
    # At each time step we predict the value h time later
  }
  return(forecast)
}

prev = forecast(model.opt, h=1) # horizon 1

plot(x1_test, type='l', xlab="Time", ylab=expression(X[1]))
lines(prev, col='red')
legend("topleft", legend=c("Observed values", "Forecasting"), lty=1,
       col=c("black","red"))

prev = lapply(c(1:10), forecast, model=model.opt)
# we make 10 forecasts: horizon 1, 2, ..., 10
# we can compare them and see the influence of the horizon

erreur = sapply(prev, function(x){mean((x-x1_test)^2)})
plot(erreur,type='b', xlab="Horizon", ylab="Mean squared error")
# the error increases with h
# it is consistent since the method is local, highly dependent on previous values
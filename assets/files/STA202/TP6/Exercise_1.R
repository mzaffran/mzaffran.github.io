rm(list=objects())

set.seed(110)
n = 1000
sigma = 1
eps = rnorm(n,0,sd=sigma)


ar = c(1,-1/4)
ma = c(-1)
x1 = arima.sim(n=n, list(ar=ar,ma=ma), innov=eps)

plot(x1)
plot(head(x1,100))

par(mfrow=c(1,2))
acf(x1)
pacf(x1)
par(mfrow=c(1,1))

qmax = 5 # ACF = 0 for h > 5
pmax = 6 # PACF


ordre = expand.grid(p=c(0:pmax), q=c(0:qmax)) # create all pairs of potentials p and q
ordre = cbind(ordre[,1], 0, ordre[,2]) # add the 0 order for differentiation in the middle

head(ordre)
dim(ordre)

?arima


model = apply(ordre, 1, arima, x=x1, method=c("ML"), SSinit=c("Rossignol2011"),
              optim.method="BFGS", include.mean=F) # 1 to specifiy that apply is by row

aic = sapply(model, function(x) {x$aic})
bic = sapply(model, function(x) {-2*x$loglik+log(x$nobs)*length(x$coef)})
like = sapply(model, function(x) {-2*x$loglik})

o = order(aic)
plot(aic[o], type='b', pch=20, axes=F, ylim=range(aic,like))
points(like[o], col='red', pch=20, type='b')
axis(1,c(1:length(aic)), paste(ordre[o,1], ordre[o,3]), las=2)
axis(2)
?axis

par(mfrow=c(1,1))
o = order(aic)
plot(aic[o[1:10]], type='b', pch=20, axes=F)
axis(1, c(1:10), paste(ordre[o[1:10],1], ordre[o[1:10],3]), las=2)
axis(2)


o = order(bic)
plot(bic[o], type='b', pch=20, axes=F)
axis(1, c(1:length(aic)), paste(ordre[o,1],ordre[o,3]), las=2)
axis(2)
ordre[o,]

## Chosen model by AIC
ordre.opt = ordre[which.min(aic),]
ordre.opt

## Chosen model by BIC
ordre.opt = ordre[which.min(bic),]
ordre.opt

model.opt = model[[which.min(bic)]]
#model.opt = arima(x=x1, order=ordre.opt, method=c("ML"),
#                 SSinit=c("Rossignol2011"),
#                 optim.method="BFGS", include.mean=F)

model.opt
names(model.opt)

## Associated parameters

model.opt$coef
model.opt$sigma2
model.opt$var.coef

## Coefficients diagnosis

## pvalue du test de student
pvalue = function(model)
{
  (1-pnorm( abs(model$coef) / sqrt(diag(model$var.coef))) )*2
}

pvalue(model.opt)

# Higher orders
model.opt_pm1 = arima(x=x1, order=ordre.opt+c(1,0,0), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_pm1)

model.opt_qm1 = arima(x=x1, order=ordre.opt+c(0,0,1), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_qm1)

# Lower orders
model.opt_pm1 = arima(x=x1, order=ordre.opt-c(1,0,0), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_pm1)

model.opt_qm1 = arima(x=x1, order=ordre.opt-c(0,0,1), method=c("ML"),
                      SSinit=c("Rossignol2011"),
                      optim.method="BFGS", include.mean=F)
pvalue(model.opt_qm1)

## Residuals analysis

plot(model.opt$residuals, type='l')

par(mfrow=c(1,2))
acf(model.opt$residuals)
pacf(model.opt$residuals)
par(mfrow=c(1,1))

qqnorm(model.opt$residuals)
qqline(model.opt$residuals, col = "steelblue", lwd = 2)

hist(model.opt$residuals, breaks=50, freq=F)
x = seq(min(model.opt$residuals),
        max(model.opt$residuals),
        length=50)
lines(x, dnorm(x,mean(model.opt$residuals),model.opt$sigma2), col='red')

## Box-Pierce test

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

# We can test ARMA(2,0) and ARMA(1,1)

pvalue_BP(model.opt_pm1, K=10)
# We reject the hypothesis of non-correlated residuals at order at most 10
pvalue_BP(model.opt_qm1, K=10) 
# We reject the hypothesis of non-correlated residuals at order at most 10

#### Forecasting

x1_train = x1[1:900]
x1_test = x1[-c(1:900)]

model.opt = arima(x=x1_train, order=ordre.opt, method=c("ML"),
                  SSinit=c("Rossignol2011"),
                  optim.method="BFGS", include.mean=F)

forecast = function(model, h)
{
  forecast = array(0, dim=100)
  for(i in c(900:999))
  {
    mod = arima(x=x1[1:(i-h)], order=ordre.opt, fixed=model$coef, include.mean=F)
    forecast[i-899] = tail(predict(mod,n.ahead=h)$pred, 1) # At each time step we predict the value h time later
  }
  return(forecast)
}

prev = forecast(model.opt, h=1)

par(mfrow=c(1,1))
plot(x1_test, type='l')
lines(prev, col='red')

prev = lapply(c(1:10), forecast, model=model.opt)

x1_test-prev

erreur = sapply(prev, function(x){mean((x-x1_test)^2)})
plot(erreur,type='b')

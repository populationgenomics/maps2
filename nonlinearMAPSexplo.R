# LT 30/9/2021


# run maps.model.2.R to set up the environment
#source()

#rational function fit to MAPS data

library(pracma)
m = lm(ps ~ mu, weights = variant_count, data=syn)



x = syn$mu
y = syn$ps


## try rationalfit
m = rationalfit(x, y, d1=3, d2=1)
ys <- polyval(m$p1,x) / polyval(m$p2,x)
plot(x, y, type="p", col="blue", ylim=c(0, 1))
points(x, Re(ys), col="red")
#grid()

# really doesnt work at all



# use 
# rational function 
f = function(par, fp=2, fq=2, x, y, bpred=FALSE)
# returns residuals unless bpred => returns predicted values
  {
  p = fp
  q = fq
  if (length(par) != p+q+1) stop('wrong parameter vector\n')
  num = par[1:(p+1)]
  # by convention firt term of denominator is 1
  den = c(1,par[(p+2):length(par)])
  
  X = do.call(cbind, lapply(0:p, function(i) x^i))
  Y = do.call(cbind, lapply(0:q, function(i) x^i))
  
  
  pred = (X %*% num) / (Y %*% den)
  y-pred
  #ifelse(bpred, pred, y-pred)
}

## try lm
m.lm33 = lm(y ~ x + I(x^2) + I(x^3) + I(x*y) + I (x^2*y) + I (x^3*y))
lm.par = c(coef(m.lm33)[1:4], -coef(m.lm33)[5:7])

plot(x, y, main='R(3,3) with lm', log='', xlab='mutability', ylab='Proportion of singletons')
lines(x, y-f(lm.par, 3, 3, x, y), col=2, lwd=2)

#plot linear fit properly (on a equispaced grid)
x.grid=seq(min(x), max(x), length.out=100)
lines(x.grid, y-f(lm.par, 3, 3, x.grid, y), col=3)
# it is terribly overfitting. should I logit/probit transform ? should I log the x-axis ?

# using Levenberg Marquardt algorithm
library(minpack.lm)

m.33 = nls.lm(c(1,0,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x, y=y, control=list(maxiter=1000))
plot(x, y, main='R(1,1) with Levenberg-Marquart', log='', xlab='mutability', ylab='Proportion of singletons')
lines(x, y-residuals(m.33), col=2, lwd=2)
# 4 parameters are 0, it is R(1,1)
lines(x.grid, y-f(coef(m.33), 3, 3, x.grid, y), col=4)

# try using starting par from lm fit
m.33.init = nls.lm(lm.par, lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x, y=y, control=list(maxiter=1000))
plot(x, y, main='R(3,3) with Levenberg-Marquart', log='')
lines(x, y-residuals(m.33.init), col=2, lwd=2)
lines(x.grid, y-f(coef(m.33.init), 3, 3, x.grid, y), col=4)

## try a couple of nls on  transformed data: y'=qnorm(y), x'= log(x)
## very simple polynomial and rational function
y2 = qnorm(y)
x2 = log(x)
x2.grid = seq(min(x2), max(x2), length.out=length(x2))

m.20 = nls(y2 ~ I(a + b * x2 + c * x2^2), start = list(a=0, b=0, c=0))
plot(x2, y2, xlab='log - mutability', ylab='probit - Proportion of singletons')
lines(x2, fitted(m.20), col=2, lwd=2)

m.30 = nls(y2 ~ I(a + b * x2 + c * x2^2 + d * x2^3), start = list(a=1, b=0, c=0, d=0), weights=syn$variant_count)
#plot(x2, y2)
lines(x2, fitted(m.30), col=3, lwd=2)

m.32 = nls.lm(c(1,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=2, x=x2, y=y2, control=list(maxiter=1000))
#lines(x2, y2-residuals(m.32), col=4)
lines(x2.grid, y2-f(coef(m.32), 3, 2, x2.grid, y2), col=4)

m.31 = nls.lm(c(1,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=1, x=x2, y=y2, control=list(maxiter=1000))
#lines(x2, y2-residuals(m.31), col=5)
lines(x2.grid, y2-f(coef(m.31), 3, 1, x2.grid, y2), col=5)

# no benefit of m.33, coefs are different but predicted values are very similar
m.33 = nls.lm(c(1,0,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x2, y=y2, control=list(maxiter=1000))
lines(x2.grid, y2-f(coef(m.33), 3, 3, x2.grid, y2), col=6)

legend('topright', legend=c("R(2,0)","R(3,0)","R(3,1)", "R(3,2)"), col=c(2,3,5,4), lwd=2)


# show R(3,2) fit on original scale
plot(x, y, main='R(3,2) with Levenberg-Marquart fit on probit/log scale', log='x', xlab='mutability', ylab='Proportion of singletons')
# line below needs loading of genome data
#plot(x, y, main='R(3,2) with Levenberg-Marquart fit on probit/log scale genomes', log='x', xlab='mutability', ylab='Proportion of singletons')
lines(exp(x2.grid), pnorm(y2-f(coef(m.32), 3, 2, x2.grid, y2)), col=4, lwd=2)
#lines(exp(x2.grid), pnorm(y2-f(coef(m.33), 3, 3, x2.grid, y2)), col=6)

#residual versus fitted
resid = y - pnorm(y2-f(coef(m.32), 3, 2, x2, y2))
plot(x, resid, log='x', ylab='Residuals', xlab='Mutability (log-scale)')  


# plot dots which size is proportional to # variants
ggplot(syn, aes(x=x, y=y, size=variant_count)) + geom_point(shape=21) #+ geom_line(data=data.frame(x=exp(x2.grid), y=pnorm(y2-f(coef(m.32), 3, 2, x2.grid, y2))))
ggplot(syn, aes(x=x, y=y, size=variant_count)) + geom_point(shape=21) +
  geom_line(data=data.frame(x=exp(x2.grid), y=pnorm(y2-f(coef(m.32), 3, 2, x2.grid, y2)), variant_count=0)) +
  xlab('Mutability') + ylab('Proportion of singletons')


# check whether we can switch to nls for weighted fit
# try on m.32
m.32.nls = nls(y2 ~ a + b * x2, start=c(a=))
# AIC calculations


# final try lasso
library(glmnet)

# a naive fit R(3,3)
# start with lm
m.33.lm = lm(y2 ~ x2 + I(x2^2) + I(x2^3) + I(x2*y2) + I (x2^2*y2) + I (x2^3*y2))
lm.par = c(coef(m.33.lm)[1:4], -coef(m.33.lm)[5:7])
lines(x2, y2-f(lm.par, 3, 3, x2, y2), col=1)
# we are back with strong overfitting
# use lasso
m.33.lasso = glmnet(x2, model.matrix(m.33.lm))
# does not work


# try on polynomial
m.30.lm = lm(y2 ~ x2 + I(x2^2) + I(x2^3))
lines(x2, cbind(1,x2, x2^2, x2^3)%*%coef(m.30.lm), col=5)
# we are back with strong overfitting
# use lasso
m.30.lasso = glmnet(model.matrix(m.30.lm), y2)

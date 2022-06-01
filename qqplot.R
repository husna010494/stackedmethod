library(car)
par(mfrow = c(1, 3), cex = 0.6)

lin.pred <- predict(mod.ln, type = "lp")[dat2$censored == 1]
lin.resid <- log(dat2$censored[dat1$censored == 1]) - lin.pred
weib.pred <- predict(mod.weib, type = "lp")[dat2$censored == 1]
weib.resid <- log(dat2$censored[dat1$censored == 1]) - weib.pred
log.pred <- predict(mod.log, type = "lp")[dat2$censored == 1]
log.resid <- log(dat2$censored[dat2$censored == 1]) - log.pred



qqPlot(exp(weib.resid), dist = "weibull",shape = 1/mod.ln$scale,
main = "Q-Q plot", xlab = "Theor. quantiles", ylab = "Emp. quantiles")

qqPlot(ln.resid, dist = "norm",sd = mod.ln$scale,main = "Q-Q plot", 
xlab = "Theor. quantiles (normal)", ylab = "Emp. quantiles")

qqPlot(log.resid, dist = "norm",sd = mod.ln$scale,main = "Q-Q plot", 
xlab = "Theor. quantiles (normal)", ylab = "Emp. quantiles")
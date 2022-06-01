
ln.predict=predict(mod.ln)
cox.predict=predict(mod.cox,type="lp")
rsf.predict=mod.rsf$predicted
stacked.predict=stacked.est$alphas[1]*ln.predict+
stacked.est$alphas[2]*cox.predict+
stacked.est$alphas[3]*rsf.predict


survivalROC(Stime        = dat2$time,
                status       = dat2$Event,
                marker       = stacked.predict,
                predict.time = 14,
                method       = "NNE",
                span = 0.25 * nrow(dat2)^(-0.20))$AUC

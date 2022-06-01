# I need to create a function that takes a list of survival function estimates and returns 
# a single estimate via stacking


# Change log (for version 2):
#	1. Implemented a continuous smoother for the alpha(t), as opposed to a k-nearest neighbor smoother.
#		This version is probably more appropriate than the nearest neighbor smoother.
#	2. I also started using the library 'alabama' to optimize the constrained least squares problem


library(alabama)


bs.data.t <- function(ys, ds, surv.lst, t){
	
# 'ys' are the observed times for all observations (assumed to be ordered)
# 'ds' are the event indicators for all observations
# 'surv.lst' is a list of all of the potential conditional survival estimates
# 't' is the time at which the Brier Score is being evaluated

	Zi <- 1 * {ys > t}		# This will be the response in the weighted LS problem
	Di <- 1 * {{ds == 1} | {ys >= t}}	# This basically says whether censored observations should be counted at t
	
	modG <- survfit(Surv(ys, I(1 - ds)) ~ 1, type = "kaplan-meier")
	
	tmpSurv <- rep(modG$surv, modG$n.event + modG$n.censor)
	
	Gc <- unlist(sapply(1:length(ys), function(x){ 
			if(mean(ys <= t) != 1){ tmp1 <- {ys[x] <= t} * tmpSurv[x] + {ys[x] > t} * tmpSurv[ys > t][1] } 
			if(mean(ys <= t) == 1){ tmp1 <- tmpSurv[x] }
			return(tmp1) }))
			
	Wi <- Di / Gc			# These are the weights for the weighted LS problem

# We are going to remove any observations with weight 0 or infinity...	
	booRemove <- ! Wi %in% c(0, NaN, Inf) 


	if(mean(ys[ds == 1] < t) != 1){ boo <- ys[ds == 1] == ys[ds == 1][ys[ds == 1] >= t][1] }
	if(mean(ys[ds == 1] < t) == 1){ boo <- c(rep(FALSE, length(ds) - 1), TRUE) }
	
	X <- matrix(unlist(lapply(surv.lst, function(x){ x[, boo] })), ncol = length(surv.lst))	# Design matrix

	rslt <- cbind(W = Wi[booRemove], Y = Zi[booRemove], X[booRemove, ])

	return(rslt)
}



stacking_ensemble <- function(ys, ds, xs, level1.SF, level0.SF, ts){
	
# 'ys' are the observed times
# 'ds' are the event indicators
# 'level1.SF' is the list of survival function estimates that EXCLUDE the obs in the fitting procedure
# 'level2.SF' is the list of survival function estimates that INCLUDE the obs in the fitting procedure
# 'ts' is the vector of times to minimize the Brier Score over

	dat1 <- NULL
	for(k in ts){
		dat1 <- rbind(dat1, bs.data.t(ys, ds, level1.SF, k))		
	}

	fr1 <- function(x){ sum(dat1[,1] * (dat1[,2] - dat1[, 3:ncol(dat1)] %*% t(t(x))) ^ 2) }
	
	gd1 <- function(x){ 
		ind <- 0
		fnl <- NULL
		while(ind < {ncol(dat1) - 2}){
			tmp <- -2 * sum(dat1[, 3 + ind] * dat1[,1] * (dat1[,2] - dat1[, 3:ncol(dat1)] %*% t(t(x)))) 
			fnl <- c(fnl, tmp)
			ind <- ind + 1
		}
		return(fnl)
	}
	
	sm1 <- function(x){
		fnl <- sum(x) - 1
		fnl
	}
	
	in1 <- function(x){
		x
	}
	
	stck1 <- auglag(rep(0, ncol(dat1) - 2), fn = fr1, gr = gd1, hin = in1, heq = sm1, control.outer = list(trace = FALSE))
	
	stacking.prob1 <- round(stck1$par, 6)


	new.ensemble1 <- matrix(0, nrow = nrow(level0.SF[[1]]), ncol = ncol(level0.SF[[1]]))
	new.ensemble2 <- matrix(0, nrow = nrow(level1.SF[[1]]), ncol = ncol(level1.SF[[1]]))
	for(k in 1:length(stacking.prob1)){
			tmp1 <- level0.SF[[k]]
			tmp2 <- level1.SF[[k]]

			new.ensemble1 <- new.ensemble1 + stacking.prob1[k] * tmp1
			new.ensemble2 <- new.ensemble2 + stacking.prob1[k] * tmp2
	}
	result <- list(Surv_predict = new.ensemble1, Surv_predict_oob = new.ensemble2, timeInterest = ys[ds == 1], 
		cens = ds, time = ys, predictors = xs, alphas = stacking.prob1)

	return(result)
}



alpha.determ <- function(Ys, Ds, Level1.SF, tt){
	
# estimates the alpha for a particular value of 't' (or, tt)
	
	dat1 <- bs.data.t(Ys, Ds, Level1.SF, tt)

	fr1 <- function(x){ sum(dat1[,1] * (dat1[,2] - dat1[, 3:ncol(dat1)] %*% t(t(x))) ^ 2) }
	
	gd1 <- function(x){ 
		ind <- 0
		fnl <- NULL
		while(ind < {ncol(dat1) - 2}){
			tmp <- -2 * sum(dat1[, 3 + ind] * dat1[,1] * (dat1[,2] - dat1[, 3:ncol(dat1)] %*% t(t(x)))) 
			fnl <- c(fnl, tmp)
			ind <- ind + 1
		}
		return(fnl)
	}
	
	sm1 <- function(x){
		fnl <- sum(x) - 1
		fnl
	}
	
	in1 <- function(x){
		x
	}
		
	stck1 <- auglag(rep(0, ncol(dat1) - 2), fn = fr1, gr = gd1, hin = in1, heq = sm1, 
		control.outer = list(trace = FALSE))
	
	stacking.prob1 <- round(stck1$par, 6)
		
	return(stacking.prob1)
}


	








#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                     #
#              Age Simulator function                 #         
#		                                      #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Version: Apr, 2022

###-----------------------------------------------------
#		VBGM_fn
###-----------------------------------------------------

	VBGM <- function(age, par, model="t0"){
		if(model=="t0"){
			L = par[1]*(1-exp(-par[2]*(age - par[3])))
		}else if(model=="l0"){
			L = par[1] - (par[1]-par[3])*exp(-par[2]*age)
		}
		
		return(L)
	}

###-----------------------------------------------------
#		Age_Simulator_fn
###-----------------------------------------------------
	
	age_sim <- function(sel_l50 = c(0,1),											
			    slope1 = 0.2, 
			    slope2 = 0.1,
			    max.sel= 1,
			    min.sel = 1,
			    nobs = 1000,
			    par = par, 
	                    model = 't0',
	                    proc.form = 'lnorm',
	                    sig_growth = sig_growth*par[1], 
			    plot=FALSE)
	 {
		#sel_l50 is the % of Linf that the first and second logistic selectivities have their inflection point
		#slope1 and slope2 are respectively the first (ascending) and second (descending) slopes of the logistic curve 
		#nobs is the # of observations
		#par is the growth parameters (Linf, K, t0 or L0)
		#model is either t0 or L0 version of VBGF
		#proc.form is either 'norm' or 'lnorm' to control whether the growth variation is in normal or lognormal space
		#sig_growth is the amount of variation around the mean in the growth curve
		#plot will plot the selectivity curve if TRUE
  	  
		lengths <- seq(0,par[1]*1.5,by=0.1)
  		min.sel <- 1-min.sel
  		logist.1 <- max.sel/(1+exp(-slope1*(lengths-(sel_l50[1]*par[1]))))
  		logist.2 <- 1-min.sel/(1+exp(-slope2*(lengths-(sel_l50[2]*par[1]))))
  		sel <- logist.1*logist.2 
  		selectivity <- data.frame(Length=lengths,selectivity=sel)
  		if(plot){
			plot(lengths,sel,type="l",lwd=2,xlab="Length",ylab="Selectivity")
  		}
  		seq.ages <- seq(0,nage,by=1)
		mean.len <- VBGM(age=seq.ages, par = par, model=model)
		age.prob <- sapply(mean.len, function(x) sel[which.min(abs(x-lengths))])
		age.prob <- pmax(0,age.prob/sum(age.prob))

		dat <- data.frame(ages=sort(apply(rmultinom(nobs, 1, prob=age.prob),2, function(x) which(x==1)))-1)
		if(proc.form=="norm"){
			sigma <- sig_growth
			dat$Len <- rnorm(nobs, VBGM(age=dat$ages, par = par, model = model), sigma)
		}else if(proc.form=="lnorm"){
			m <- VBGM(age=dat$ages, par = par, model = model)
			sd <- sig_growth
			m2 <- log(m^2 / sqrt(sd^2 + m^2))
			sd2 <- sqrt(log(1 + (sd^2 / m^2)))
			r <- c(scale(rnorm(nobs)))
			sigma <- sd/par[1]
			dat$Len <- rlnorm(r, m2, sd2)
		}
  		return(list(df=dat,
		            par = c(par,sigma),
		            selectivity = selectivity,
		            true.len = data.frame(age=seq.ages, len=mean.len)))
	}

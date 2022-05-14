#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                     #
#    		Age-growth Bias Simulator             #         
#		                                      #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Modified: Fabio Caltabellotta
# Version Modified: May, 2022

###-----------------------------------------------------
#		Simulation
###-----------------------------------------------------
	
	library(truncnorm)
	library(mvtnorm)
	library(rstan)
	library(foreach)
	library(parallel)
	library(doParallel)  
	options(mc.cores = parallel::detectCores())
	rstan_options(auto_write = TRUE)

	source("age_simulator_fn.R", local = FALSE)
	source("growth_stan_model.R", local = FALSE)
	source("vbgm_stan_fn.R", local=FALSE)

	###-----------------------------------------------------
	
	# true growth parameters

	Linf = 64.0831
	K = 0.0649209
	t0 = -1.22137

	par = c(Linf, K, t0)

	# number of ages
	
	nage = 118
	ages = 0:(nage-1)

	#----------------------------------------------------
	#		Simulating Scenarios
	#-----------------------------------------------------
	
	sim_unbiased <- foreach(sigma=rep(c(0.02,0.05,0.10),10)) %do% age_sim(sel_l50 = c(0,1),
						  nobs = 200,
						  slope1 = 1, 
					      slope2 = 1,
						  max.sel = 1,
						  min.sel = 1,
					      par = par, 
	                      model = 't0',
	                      proc.form = 'lnorm',
	                      sig_growth = sigma*par[1], 
						  plot = FALSE)

	sim_nsmall <- foreach(sigma=rep(c(0.02,0.05,0.10),10)) %do% age_sim(sel_l50 = c(0.7,1),
						  nobs = 200,
						  slope1 = 1, 
						  slope2 = 1,
						  max.sel = 1,
						  min.sel = 1,
						  par = par, 
	                      model = 't0',
	                      proc.form = 'lnorm',
	                      sig_growth = sigma*par[1], 
						  plot = FALSE)

	sim_nlarge <- foreach(sigma=rep(c(0.02,0.05,0.10),10)) %do% age_sim(sel_l50 = c(0,0.7),
						  nobs = 200,
						  slope1 = 1, 
						  slope2 = 1,
						  max.sel = 1,
						  min.sel = 0,
						  par = par, 
	                      model = 't0',
	                      proc.form = 'lnorm',
	                      sig_growth = sigma*par[1], 
						  plot = FALSE)

	sim_bloated <- foreach(sigma=rep(c(0.02,0.05,0.10),10)) %do% age_sim(sel_l50 = c(0.55,0.85),
						  nobs = 200,
						  slope1 = 1, 
						  slope2 = 1,
						  max.sel = 1,
						  min.sel = 0,
						  par = par, 
	                      model = 't0',
	                      proc.form = 'lnorm',
	                      sig_growth = sigma*par[1], 
						  plot = FALSE)

		unbiased_fit <- vector(mode="list", length = length(sim_unbiased))
		nsmall_fit <- vector(mode="list", length = length(sim_nsmall))
		nlarge_fit <- vector(mode="list", length = length(sim_nlarge))
		bloated_fit <- vector(mode="list", length = length(sim_bloated))

system.time({
			for (i in 1:30){
				unbiased_fit[[i]] <- vbgm_stan_fn(sim_unbiased[[i]], par=sim_unbiased[[i]]$par, model="t0")			
				nsmall_fit[[i]] <- vbgm_stan_fn(sim_nsmall[[i]], par=sim_nsmall[[i]]$par, model="t0")
				nlarge_fit[[i]] <- vbgm_stan_fn(sim_nlarge[[i]], par=sim_nlarge[[i]]$par, model="t0")
				bloated_fit[[i]] <- vbgm_stan_fn(sim_bloated[[i]], par=sim_bloated[[i]]$par, model="t0")
			}	
	})

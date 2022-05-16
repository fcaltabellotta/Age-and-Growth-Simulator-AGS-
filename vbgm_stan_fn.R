#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                     #
#    		   Growth_Stan_Model                  #         
#		                                      #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Modified: Fabio Caltabellotta
# Version Modified: May, 2022

###-----------------------------------------------------
#		VBGM_Stan_Function
###-----------------------------------------------------
	
	vbgm_stan_fn <- function(sim_dat, par, model="t0", nchains=4, nwarmup = 6000, niter = 7500, qprobs = c(0.05,0.5,0.95)){
		
		#fit the par_guess
			vbgm_optim <- function(theta, ages, len, model){
				linf <- exp(theta[1])
				k <- exp(theta[2])
				if(model=="t0"){
					t0 <- theta[3]
				}else if(model=="l0"){
					l0 <- exp(theta[3])
				}
				sigma <- exp(theta[4])

				if(model=="t0"){
					l_est <- linf*(1-exp(-k*(ages - t0)))
				}else if(model=="l0"){
					l_est <- linf-(linf-l0)*exp(-k*ages)
				}		

				#NLL = -1*sum(dlnorm(len, l_est, sigma, log=TRUE))
				NLL = -1*sum(dnorm(len, l_est, sigma, log=TRUE))

				return(NLL)
			}
			if(model=="t0"){
				theta <- c(log(par[1]), log(par[2]), 
				           par[3], log(1))
			}else if(model=="l0"){
				theta <- c(log(par[1]), log(par[2]), 
				           log(par[3]), log(1))
			}
			
			vbgm_optim_fit <- optim(par = theta,
			                        fn = vbgm_optim, 
			                        ages = sim_dat$df$ages, 
			                        len = sim_dat$df$Len, 
			                        model = model)
			if(model=='t0'){
				par_guess <- with(vbgm_optim_fit, c(exp(par[1]), exp(par[2]), par[3], exp(par[4])))
			}else if(model=="l0"){
				par_guess <- with(vbgm_optim_fit, c(exp(par[1]), exp(par[2]), exp(par[3])))
			}
			
		#init functions
			inits.t0 <- function(chain_id){
		
				sig_growth <- rtruncnorm(1,a=0,b=Inf,
				                         mean=par_guess[4],
				                         sd=par_guess[4]*0.5)

				vbgm <- rtruncnorm(3, a=c(0,0,-Inf),
				                   b = c(Inf,Inf,0),
				                   par_guess, abs(par_guess*0.1))

				return(list(Linf = vbgm[1], 
				            K=vbgm[2], 
				            tknot=vbgm[3], 
				            sigma_growth = sig_growth))
			}
			inits.l0 <- function(chain_id){
				
				sig_growth <- rtruncnorm(1,a=0,b=Inf,
				                         mean=par_guess[4],
				                         sd=par_guess[4]*0.1)

				vbgm <- rtruncnorm(3, a=0, par_guess, abs(par_guess*0.1))

				return(list(Linf = vbgm[1], 
				            K = vbgm[2], 
				            lknot = vbgm[3], 
				            sigma_growth = sig_growth))
			}
		#compile the data
			stan_dat <- list(n_obs = nrow(sim_dat$df),
			                 age = sim_dat$df$ages,
			                 l = sim_dat$df$Len,
			                 par_guess = par_guess[1:3],
			                 cv_guess = 0.8,
			                 t_par = sim_dat$par,
			                 n_pred = nrow(sim_dat$true.len),
			                 p_age = sim_dat$true.len$age,
			                 t_len = sim_dat$true.len$len,
					 nseq = length(seq(0.1,max(sim_dat$df$ages),by=0.1)),
					 seq_ages = seq(0.1,max(sim_dat$df$ages),by=0.1))
		#run the model
			if(model=="t0"){
				vbgm_stan_fit <- sampling(object = vbgm_stan,
	                      data = stan_dat,
	                      chains = 4,
	                      warmup = 6000,
	                      iter = 7500,
	                      init = inits.t0,
	                      control = list(adapt_delta = 0.95),
	                      pars = c('Linf','K','tknot','sig_growth','Linf_bias','K_bias','tknot_bias','sig_growth_bias','Lage','p_len'),
	                      open_progress=FALSE)
			}else if(model=="l0"){
				vbgm_stan_fit <- sampling(object = vbgm_stan,
	                      data = stan_dat,
	                      chains = nchains,
	                      warmup = nwarmup,
	                      iter = niter,
	                      init = inits.l0,
	                      control = list(adapt_delta = 0.95),
	                      pars = c('Linf','K','lknot','sig_growth','Linf_bias','K_bias','lknot_bias','sig_growth_bias','Lage','p_len'),
	                      open_progress=FALSE)
			}
			samples <- rstan::extract(vbgm_stan_fit)
			quants_len = apply(samples[[10]], 2, quantile, probs=qprobs, na.rm=TRUE)			
			colnames(quants_len) <- paste0('A',sim_dat$true.len$age)
		#return values
			ret <- list(model = model,
			            data = stan_dat,
			            fit = vbgm_stan_fit,
			            quants = sapply(samples[1:9], quantile, probs=qprobs, na.rm=TRUE),
			            quants_len = quants_len,			            
			            samples = samples)
			return(ret)
	}

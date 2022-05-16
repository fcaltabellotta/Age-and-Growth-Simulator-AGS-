#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                                     #
#    		   Growth_Stan_Model                  #         
#		                                      #
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Modified: Fabio Caltabellotta
# Version Modified: Apr, 2022

###-----------------------------------------------------
#		VBGM_Stan_Models
###-----------------------------------------------------

	vbgm_stan_setup <- function(model='t0'){
		stan_data_sect <- "
			data{
				int<lower=1> n_obs; //number of observations
				vector<lower=0>[n_obs] age; //ages of fish
				vector<lower=0>[n_obs] l; //length of fish
				vector[3] par_guess; //guesses on Linf, K, tknot
				real<lower=0> cv_guess;
				vector[4] t_par; //true parameter values
				int<lower=1> n_pred; //number of predicted ages/lengths
				vector<lower=0>[n_pred] p_age; //ages to predict to
				vector<lower=0>[n_pred] t_len; //true lengths
				int<lower=1> nseq;
				vector<lower=0>[nseq] seq_ages;
			}"
			if(model=="t0"){
				stan_parm_sect <-"
				parameters{
					real<lower=0> Linf; //L infinity
					real<lower=0> K; //vb K
					real<upper=0> tknot; //t knot					
					real<lower=0> sig_growth; //sigma growth
				}"
				stan_mod_sect <- "
				model{
					vector[n_obs] l_pred;
					Linf ~ normal(par_guess[1], par_guess[1]*cv_guess);
					K ~ normal(par_guess[2], par_guess[2]*cv_guess);
					tknot ~ normal(par_guess[3], fabs(par_guess[3]*cv_guess));

					sig_growth ~ normal(0, 1);

					for(i in 1:n_obs){
						l_pred[i] = log(Linf*(1-exp(-K*(age[i]-tknot))));
						l[i] ~ lognormal(l_pred[i], sig_growth);
					}			
				}"
				stan_gq_sect <- "
				generated quantities{
					vector[n_pred] p_len; 
					real Linf_bias;
					real K_bias;
					real tknot_bias;
					real sig_growth_bias;
					vector[nseq] Lage;
					

					//standardized bias for length at age
					for(i in 1:n_pred){
						p_len[i] = (t_len[i]-(lognormal_rng(log(Linf*(1-exp(-K*(p_age[i]-tknot)))), sig_growth)))/t_len[i];
					}

					//generating credible intervals
					for(i in 1:nseq){
						Lage[i] = lognormal_rng(log(Linf*(1-exp(-K*(seq_ages[i]-tknot)))),sig_growth);
					}

					//relative bias
					Linf_bias =  (Linf-t_par[1])/t_par[1];
					K_bias = (K-t_par[2])/t_par[2];
					tknot_bias = (tknot-t_par[3])/t_par[3];
					sig_growth_bias = (sig_growth-t_par[4])/t_par[4];
					
				}"
			}else if(model=="l0"){
				stan_parm_sect <- "
				parameters{
					real<lower=0> Linf; //L infinity
					real<lower=0> K; //vb K
					real<lower=0> lknot; //l knot
					real<lower=0> sig_growth; //sigma growth
					int<lower=1> nseq;
					vector<lower=0>[nseq] seq_ages;
				}"
				stan_mod_sect <- "
				model{
					vector[n_obs] l_pred;
					Linf ~ normal(par_guess[1], par_guess[1]*cv_guess);
					K ~ normal(par_guess[2], par_guess[2]*cv_guess);
					lknot ~ normal(par_guess[3], par_guess[3]*cv_guess);

					sig_growth ~ normal(0, 1);

					for(i in 1:n_obs){
						l_pred[i] = log(Linf-(Linf-lknot)*exp(-K*age[i]));
						l[i] ~ lognormal(l_pred[i], sig_growth);
					}
				}"
				stan_gq_sect <- "
				generated quantities{
					vector[n_pred] p_len;
					real Linf_bias;
					real K_bias;
					real lknot_bias;
					real sig_growth_bias;
					vector[nseq] Lage;

					//standardized bias for length at age
					for(i in 1:n_pred){
						p_len[i] = (t_len[i]-(lognormal_rng(log(Linf-(Linf-lknot)*exp(-K*p_age[i])), 
						sig_growth))])/t_len[i];
					}

					//generating credible intervals
					for(i in 1:nseq){
						Lage[i] = lognormal_rng(log(Linf-(Linf-lknot)*exp(-K*seq_ages[i])),sig_growth);
					}
					
					//relative bias 
					Linf_bias =  (Linf-t_par[1])/t_par[1];
					K_bias = (K-t_par[2])/t_par[2];
					lknot_bias = (lknot-t_par[3])/t_par[3];
					sig_growth_bias = (sig_growth-t_par[4])/t_par[4];
					
				}"
			}
			vbgm_stan <- paste(stan_data_sect,
			                   stan_parm_sect,
			                   stan_mod_sect,
			                   stan_gq_sect)
			mod <- stan_model(model_code=vbgm_stan)
			return(mod)
	}

	vbgm_stan <- vbgm_stan_setup(model='t0')	

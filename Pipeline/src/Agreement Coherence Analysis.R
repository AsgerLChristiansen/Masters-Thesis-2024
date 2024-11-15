# Functions for analysis of gun violence


pacman::p_load(pacman, tidyverse, stringr, LaplacesDemon)

#### PREPROCESSING FUNCTIONS

simulate_pair_agreement= function(max_agreement_change = 20, real = 1){
  
  if (real == 1){
  
  possible_agreement_values = c(40:60)/100
  possible_agreement_change = c(0:max_agreement_change)/100
  
  v1_agreement = sample(possible_agreement_values,1)
  direction = sample(c(1, -1),1)
  v2_agreement = v1_agreement + sample(possible_agreement_change,1)*direction
  }
  
  
  else{
    v1_agreement = sample(c(20:80)/100,1)
    v2_agreement = sample(c(20:80)/100,1)
    
  }
  
  return(c(v1_agreement, v2_agreement))
  
  
}


simulate_condition_and_visits_2 = function(mu_beta_1_channel, 
                                         sigma_beta_1_channel,
                                         mu_beta_2_channel,
                                         sigma_beta_2_channel,
                                         mu_beta_3_channel,
                                         sigma_beta_3_channel,
                                         pair_agreement,
                                         N_conditions = 2, loc =1, scale = 3){
  

  results = c()
  
  visits = c(0,1)
  condition_modifiers = c(1, 1.5)
  for (c in 1:N_conditions){
    # Generate  condition parameters
    beta_1 =  rnorm(n = 1,mean = mu_beta_1_channel, sd =  sigma_beta_1_channel)*condition_modifiers[c]
    beta_2 =  rnorm(n = 1,mean = mu_beta_2_channel, sd =  sigma_beta_2_channel)*condition_modifiers[c]
    beta_3 =  rnorm(n = 1,mean = mu_beta_3_channel, sd =  sigma_beta_3_channel)*condition_modifiers[c]
    sigma = rgamma(n = 1,loc, scale)
    visit_results = c()
    for (index in c(1,2)){
      
      mu =  beta_1 - visits[index]*beta_1 + visits[index]*beta_2 + beta_3*pair_agreement[index]
          
      y = invlogit(rnorm(n = 1, mean = mu, sd = sigma))
          
      visit_results = c(visit_results, y)
          
        }
      
    results = array(c(results, visit_results), dim = c(2, c))
        }
        return(results)
      }
      
  


simulate_channel_data_2 = function(mu_beta_1_pair, sigma_beta_1_pair,
                                 mu_beta_2_pair, sigma_beta_2_pair, 
                                 mu_beta_3_pair, sigma_beta_3_pair, 
                                 pair_agreement,  N_channels, N_conditions = 2, N_visits = 2,
                                 loc =1, scale = 3){
  results = c()
  for (c in 1:N_channels){
      
      mu_beta_1_channel = rnorm(n = 1,mean = mu_beta_1_pair, sd =  sigma_beta_1_pair)
      sigma_beta_1_channel = rgamma(n = 1, loc, scale)
      mu_beta_2_channel =  rnorm(n = 1, mean = mu_beta_2_pair, sd =  sigma_beta_2_pair)
      sigma_beta_2_channel = rgamma(n = 1,loc, scale)
      mu_beta_3_channel =  rnorm(n = 1, mean = mu_beta_3_pair, sd =  sigma_beta_3_pair)
      sigma_beta_3_channel = rgamma(n = 1,loc, scale)
      
      
      y = simulate_condition_and_visits_2(mu_beta_1_channel, sigma_beta_1_channel,
                                        mu_beta_2_channel, sigma_beta_2_channel,
                                        mu_beta_3_channel, sigma_beta_3_channel,
                                        pair_agreement,
                                        loc = loc, scale = scale,N_conditions = 2)
      
      
      dims = dim(y)
      results = array(c(results, y), dim = c(dims, c))
      
    }
    
    return(results)
}

simulate_pairs_data_2 = function(mu_beta_1_real, sigma_beta_1_real, mu_beta_1_surrogate, sigma_beta_1_surrogate,
                         mu_beta_2_real, sigma_beta_2_real, mu_beta_2_surrogate, sigma_beta_2_surrogate,
                         mu_beta_3_real, sigma_beta_3_real, mu_beta_3_surrogate, sigma_beta_3_surrogate,
                         N_pairs,
                         N_pairs_real,
                         loc =1, scale =3, N_channels = 33){

  
    
  
    results = c()
    agreement = c()
    for (p in 1:N_pairs){
      
      if (p <= N_pairs_real){
        # Simulating this pair's agreement
        pair_agreement = simulate_pair_agreement(real = 1)
        
        mu_beta_1_pair = rnorm(n = 1,mean = mu_beta_1_real, sd =  sigma_beta_1_real)
        sigma_beta_1_pair = rgamma(n = 1, loc, scale)
        mu_beta_2_pair =  rnorm(n = 1, mean = mu_beta_2_real, sd =  sigma_beta_2_real)
        sigma_beta_2_pair = rgamma(n = 1, loc, scale)
        mu_beta_3_pair =  rnorm(n = 1, mean = mu_beta_3_real, sd =  sigma_beta_3_real)
        sigma_beta_3_pair = rgamma(n = 1, loc, scale)
        
        
        
        y = simulate_channel_data_2(mu_beta_1_pair, sigma_beta_1_pair,
                                  mu_beta_2_pair, sigma_beta_2_pair, 
                                  mu_beta_3_pair, sigma_beta_3_pair, 
                                  pair_agreement, N_channels,
                                  loc = loc, scale = scale,  N_conditions = 2)
        
        
        
      }
      
      else{
        
        # Simulating this pair's agreement
        pair_agreement = simulate_pair_agreement(real = 0) %>% array(dim = c(1,2))
        
        
        mu_beta_1_pair = rnorm(n = 1,
        mean = mu_beta_1_surrogate, sd = sigma_beta_1_surrogate)
        
        sigma_beta_1_pair = rgamma(n = 1, loc, scale)
        mu_beta_2_pair = rnorm(n = 1,mean = mu_beta_2_surrogate, sd = sigma_beta_2_surrogate)
        sigma_beta_2_pair = rgamma(n = 1, loc, scale)
        
        mu_beta_3_pair = rnorm(n = 1, 
        mean = mu_beta_3_surrogate, sd = sigma_beta_3_surrogate)
        
        sigma_beta_3_pair = rgamma(n = 1, loc, scale)
        
        y = simulate_channel_data_2(mu_beta_1_pair, sigma_beta_1_pair,
                                  mu_beta_2_pair, sigma_beta_2_pair,
                                  mu_beta_3_pair, sigma_beta_3_pair,
                                  pair_agreement,
                                  loc = loc, scale = scale, N_channels,  N_conditions = 2)
        # y should be a visit x condition x channels array
      }
      dims = dim(y)
      results = array(c(results, y), c(dims, p))
      agreement = array(c(agreement, pair_agreement), c(p, 2))
      
    }
    
    out = list("agreement" = agreement, "results" = results)
  return(out)
}

simulate_data_2 = function(mu_beta_1_real_vec,
                         sigma_beta_1_real_vec,
                         mu_beta_1_surrogate_vec,
                         sigma_beta_1_surrogate_vec,
                         
                         mu_beta_2_real_vec, 
                         sigma_beta_2_real_vec,
                         mu_beta_2_surrogate_vec, 
                         sigma_beta_2_surrogate_vec, 
                         
                         mu_beta_3_real_vec, 
                         sigma_beta_3_real_vec,
                         mu_beta_3_surrogate_vec, 
                         sigma_beta_3_surrogate_vec, 
                         
                         loc =1, scale = 3, N_sims = 2, N_pairs = 44, 
                         N_pairs_real = 22, N_channels = 33){
  
  
  all_sims = c()
  all_agreements = c()
  
  parameters = c()
  for (sim in 1:N_sims){
    mu_beta_1_real = sample(mu_beta_1_real_vec, 1)
    sigma_beta_1_real = sample(sigma_beta_1_real_vec, 1)
    mu_beta_1_surrogate = sample(mu_beta_1_surrogate_vec, 1)
    sigma_beta_1_surrogate = sample(sigma_beta_1_surrogate_vec, 1)
    
    mu_beta_2_real = sample(mu_beta_2_real_vec, 1) 
    sigma_beta_2_real= sample(sigma_beta_2_real_vec, 1)
    mu_beta_2_surrogate = sample(mu_beta_2_surrogate_vec, 1) 
    sigma_beta_2_surrogate = sample(sigma_beta_2_surrogate_vec, 1)
    
    mu_beta_3_real = sample(mu_beta_3_real_vec, 1) 
    sigma_beta_3_real= sample(sigma_beta_3_real_vec, 1)
    mu_beta_3_surrogate = sample(mu_beta_3_surrogate_vec, 1) 
    sigma_beta_3_surrogate = sample(sigma_beta_3_surrogate_vec, 1)
    
    simulated_pair_data = simulate_pairs_data_2(mu_beta_1_real, sigma_beta_1_real, 
                                  mu_beta_1_surrogate, sigma_beta_1_surrogate,
                                  mu_beta_2_real, sigma_beta_2_real,
                                  mu_beta_2_surrogate, sigma_beta_2_surrogate,
                                  mu_beta_3_real, sigma_beta_3_real,
                                  mu_beta_3_surrogate, sigma_beta_3_surrogate,
                                  N_pairs,
                                  N_pairs_real, loc, scale, N_channels)
    
  results = simulated_pair_data$results
  agreement = simulated_pair_data$agreement
  
  all_agreements = array(c(all_agreements, agreement), c(dim(agreement), sim))
    
  all_sims = array(c(all_sims, results), c(dim(results), sim))
  
  
  parameters = array(c(parameters, c(mu_beta_1_real, sigma_beta_1_real, 
                                     mu_beta_1_surrogate, sigma_beta_1_surrogate,
                                     mu_beta_2_real, sigma_beta_2_real, 
                                     mu_beta_2_surrogate, sigma_beta_2_surrogate,
                                     mu_beta_3_real, sigma_beta_3_real, 
                                     mu_beta_3_surrogate, sigma_beta_3_surrogate)),c(12, sim))
  
  
  }
  
  # Turn parameters into dataframe
  
  parameters_df = parameters %>% t() %>% as.data.frame()
  colnames(parameters_df) = c("mu_beta_1_real", "sigma_beta_1_real", 
                              "mu_beta_1_surrogate", "sigma_beta_1_surrogate",
                              "mu_beta_2_real", "sigma_beta_2_real", 
                              "mu_beta_2_surrogate", "sigma_beta_2_surrogate",
                              "mu_beta_3_real", "sigma_beta_3_real", 
                              "mu_beta_3_surrogate", "sigma_beta_3_surrogate"
  )
  
  out = list("parameters" = parameters_df, "sims_array" = all_sims, "agreement_array" = all_agreements)
  return(out)
}




parameter_recovery_agreement = function(simulated_data_array, parameters, stan_model, chains, parallel_chains, iter_warmup, iter_sampling, filenames, directory = "ParameterRecovery", refresh = 5, seed = 1337){
  
  write_csv(as.data.frame(parameters), paste(directory, "/parameters.csv", sep = ""))
  
  dims = dim(simulated_data_array)
  
  N_visits = dims[1]
  N_conditions = dims[2]
  N_channels = dims[3]
  N_pairs = dims[4]
  
  N_sims = dims[5]
  print(dims)
  print(N_sims)
  print(length(filenames))
  N_pairs_real = N_pairs/2
  N_pairs_surrogate = N_pairs_real
  
  options(mc.cores=4)
  
  for (i in 1:N_sims){
    coherence = simulated_data_array[,,,,i]
    
    
    
    data <- list(N_channels = N_channels,
                 N_visits=N_visits,
                 N_pairs=N_pairs,
                 N_pairs_real=N_pairs_real,
                 N_pairs_surrogate=N_pairs_surrogate,
                 N_conditions=N_conditions,
                 coherence=coherence, visit = c(0,1))
    
    # Run model
    samples <- mod_coherence$sample(
      data = data,
      seed = 1000,
      chains = 3,
      parallel_chains = 3,
      threads_per_chain = 1,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = refresh,
      max_treedepth = 20,
      adapt_delta = 0.99,)
    
    samples$save_object(file =paste(directory, "/parameter_recovery_", filenames[i],".RDS", sep = "") )
    
    
  }
  
}

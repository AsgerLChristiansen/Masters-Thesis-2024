pacman::p_load(pacman, tidyverse, stringr, LaplacesDemon)

#### PREPROCESSING FUNCTIONS

# Helper function to turn a row in a dataframe into a vector instead of
# a list.
row_to_vector = function(dataframe_row){
  
  new_vector = unlist(c(dataframe_row))
  new_vector = new_vector[!is.na(new_vector)]
  return(new_vector)
}



simulate_condition_and_visits = function(mu_beta_1_channel, 
                                         sigma_beta_1_channel,
                                         mu_beta_2_channel,
                                         sigma_beta_2_channel,
                                         N_conditions = 2, loc =0.4, scale = 4){
  
  
  results = c()
  for (c in 1:N_conditions){
    # Generate  condition parameters
    beta_1 =  rnorm(n = 1,mean = mu_beta_1_channel, sd =  sigma_beta_1_channel)
    beta_2 =  rnorm(n = 1,mean = mu_beta_2_channel, sd =  sigma_beta_2_channel)
    sigma = rgamma(n = 1,loc, scale)
    visit_results = c()
    for (v in c(0,1)){
      
      mu =  beta_1 - v*beta_1 + v*beta_2
          
      y = invlogit(rnorm(n = 1, mean = mu, sd = sigma))
          
      visit_results = c(visit_results, y)
          
        }
      
    results = array(c(results, visit_results), dim = c(2, c))
        }
        return(results)
      }
      
  

dummy_code = function(vector){
  
  
  new_levels = c(1:length(unique(vector)))
  
  levels(vector) <- new_levels
  
  return(vector)
}


simulate_channel_data = function(mu_beta_1_pair, sigma_beta_1_pair, mu_beta_2_pair,
                                 sigma_beta_2_pair, N_channels = 33, N_conditions = 2, N_visits = 2,
                                 loc =0.4, scale = 4){
  
  results = c()
  for (c in 1:N_channels){
      
      mu_beta_1_channel = rnorm(n = 1,mean = mu_beta_1_pair, sd =  sigma_beta_1_pair)
      sigma_beta_1_channel = rgamma(n = 1, loc, scale)
      mu_beta_2_channel =  rnorm(n = 1, mean = mu_beta_2_pair, sd =  sigma_beta_2_pair)
      sigma_beta_2_channel = rgamma(n = 1,loc, scale)
      
      
      y = simulate_condition_and_visits(mu_beta_1_channel, sigma_beta_1_channel, mu_beta_2_channel, sigma_beta_2_channel,
                                        loc = loc, scale = scale,N_conditions = 2)
      
      
      dims = dim(y)
      results = array(c(results, y), dim = c(dims, c))
      
    }
    
    return(results)
}

simulate_pairs_data = function(mu_beta_1_real, sigma_beta_1_real, mu_beta_1_surrogate, sigma_beta_1_surrogate,
                         mu_beta_2_real, sigma_beta_2_real, mu_beta_2_surrogate, sigma_beta_2_surrogate,
                         N_pairs,
                         N_pairs_real,
                         loc =0.4, scale = 4, N_channels = 33){

    results = c()
    for (p in 1:N_pairs){
      
      if (p <= N_pairs_real){
        
        mu_beta_1_pair = rnorm(n = 1,mean = mu_beta_1_real, sd =  sigma_beta_1_real)
        sigma_beta_1_pair = rgamma(n = 1, loc, scale)
        mu_beta_2_pair =  rnorm(n = 1, mean = mu_beta_2_real, sd =  sigma_beta_2_real)
        sigma_beta_2_pair = rgamma(n = 1,loc, scale)
        
        
        y = simulate_channel_data(mu_beta_1_pair, sigma_beta_1_pair,
                                  mu_beta_2_pair, sigma_beta_2_pair, loc = loc, scale = scale, N_channels = N_channels,  N_conditions = 2)
        
        
        
      }
      
      else{
        mu_beta_1_pair = rnorm(n = 1,mean = mu_beta_1_surrogate, sd = sigma_beta_1_surrogate)
        sigma_beta_1_pair = rgamma(n = 1, loc, scale)
        mu_beta_2_pair = rnorm(n = 1, mean = mu_beta_2_surrogate, sd = sigma_beta_2_surrogate)
        
        y = simulate_channel_data(mu_beta_1_pair, sigma_beta_1_pair, mu_beta_2_pair, sigma_beta_2_pair)
        # y should be a visit by conditions by channels array
      }
      dims = dim(y)
      results = array(c(results, y), c(dims, p))
      
    }
    
  
  return(results)
}

simulate_data = function(mu_beta_1_real_vec,
                         sigma_beta_1_real_vec,
                         mu_beta_2_real_vec, 
                         sigma_beta_2_real_vec,
                         mu_beta_1_surrogate_vec,
                         sigma_beta_1_surrogate_vec,
                         mu_beta_2_surrogate_vec, 
                         sigma_beta_2_surrogate_vec, 
                         loc =0.4, scale = 4, N_sims = 2, N_pairs = 44, 
                         N_pairs_real = 22, N_channels = 33){
  
  
  all_sims = c()
  
  parameters = c()
  for (sim in 1:N_sims){
    mu_beta_1_real = sample(mu_beta_1_real_vec, 1)
    mu_beta_2_real = sample(mu_beta_2_real_vec, 1) 
    sigma_beta_1_real = sample(sigma_beta_1_real_vec, 1)
    sigma_beta_2_real= sample(sigma_beta_2_real_vec, 1)
    mu_beta_1_surrogate = sample(mu_beta_1_surrogate_vec, 1)
    mu_beta_2_surrogate = sample(mu_beta_2_surrogate_vec, 1) 
    sigma_beta_1_surrogate = sample(sigma_beta_1_surrogate_vec, 1)
    sigma_beta_2_surrogate = sample(sigma_beta_2_surrogate_vec, 1)
    
    results = simulate_pairs_data(mu_beta_1_real, sigma_beta_1_real, 
                                  mu_beta_1_surrogate, sigma_beta_1_surrogate,
                                  mu_beta_2_real, sigma_beta_2_real,
                                  mu_beta_2_surrogate, sigma_beta_2_surrogate,
                                  N_pairs,
                                  N_pairs_real, loc, scale, N_channels = N_channels)
    
  
  all_sims = array(c(all_sims, results), c(dim(results), sim))
  
  
  parameters = array(c(parameters, c(mu_beta_1_real, sigma_beta_1_real, 
                                     mu_beta_1_surrogate, sigma_beta_1_surrogate,
                                     mu_beta_2_real, sigma_beta_2_real, 
                                     mu_beta_2_surrogate, sigma_beta_2_surrogate)),c(8, sim))
  
  
  
  }
  parameters = t(parameters)%>% as.data.frame()
  colnames(parameters) = c("mu_beta_1_real", "sigma_beta_1_real", 
                              "mu_beta_1_surrogate", "sigma_beta_1_surrogate",
                              "mu_beta_2_real", "sigma_beta_2_real", 
                              "mu_beta_2_surrogate", "sigma_beta_2_surrogate"  )
  out = list("parameters" = parameters, "sims_array" = all_sims)
  return(out)
}




get_est = function(summary, var_name){
  est = unlist(c(summary[summary$variable==var_name,"mean"]))
  return(est)
}


############## ANALYSIS FUNCTIONS #####################


evaluate_dimension = function(dataset_column){
  
  return(dataset_column %>% unique() %>% length() %>% as.numeric())
}



prepare_for_analysis = function(dataset){
  N_visits = evaluate_dimension(RelevantHbo$visit) # {0, 1}, representing first and fourth conversation.
  N_conditions = evaluate_dimension(RelevantHbo$condition) # {1, 2}, representing rest and talking
  N_channels = evaluate_dimension(RelevantHbo$channel) # {1:33}.
  N_pairs =  evaluate_dimension(RelevantHbo$pair) # 1:44
  
  print(N_visits*N_conditions*N_channels*N_pairs)
  
  dataset = dataset %>%  arrange(pair, channel, condition, visit)
  
  data_array = array(data = NA, dim = c(N_visits, N_conditions, N_channels, N_pairs))
  
  
  row = 0
  for (p in 1:N_pairs){
    for (ch in 1:N_channels){
      for (c in 1:N_conditions){ 
        for (v in 1:N_visits){
          row = row +1
          data_array[v,c,ch, p] = dataset$coherence[row]
          
        }}}}
  
  
  return(data_array)
  
}

fit_to_data = function(data_array, stan_model, chains, parallel_chains, iter_warmup, iter_sampling, filename = "test", directory = "Results", refresh = 15){
  
  dims = dim(data_array)
  
  N_visits = dims[1]
  N_conditions = dims[2]
  N_channels = dims[3]
  N_pairs = dims[4]
  
  N_pairs_real = N_pairs/2
  N_pairs_surrogate = N_pairs_real
  
  options(mc.cores=4)
  
  
  
  data <- list(N_channels = N_channels,
               N_visits=N_visits,
               N_pairs=N_pairs,
               N_pairs_real=N_pairs_real,
               N_pairs_surrogate=N_pairs_surrogate,
               N_conditions=N_conditions,
               coherence=data_array, visit = c(0,1))
  
  # Run model
  samples <- mod_coherence$sample(
    data = data,
    seed = 1000,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    max_treedepth = 20,
    adapt_delta = 0.99,)
  
  samples$save_object(file =paste(directory, "/Results_", filename,".RDS", sep = "") )
  
}



convert_to_coherence = function(dataframe, column_indeces){
  
  for (index in column_indeces){
    
    dataframe[,index] = dataframe[,index] %>%  invlogit()
    
    
  }
  return(dataframe)
}



sanity_check = function(dataframe, condition){
  channel = c()
  beta_1_ER = c()
  beta_2_ER = c()
  beta_change_ER = c()
  beta_1_cred = c()
  beta_2_cred = c()
  beta_change_cred = c()
  
  map_beta_1_real = c()
  map_beta_1_surrogate = c()
  high_ci_beta_1_real = c()
  high_ci_beta_1_surrogate = c()
  low_ci_beta_1_real =c()
  low_ci_beta_1_surrogate = c()
  
  map_beta_2_real = c()
  map_beta_2_surrogate = c()
  high_ci_beta_2_real = c()
  high_ci_beta_2_surrogate = c()
  low_ci_beta_2_real = c()
  low_ci_beta_2_surrogate =c()
  
  
  
  map_beta_change_real = c()
  map_beta_change_surrogate = c()
  high_ci_beta_change_real = c()
  high_ci_beta_change_surrogate = c()
  low_ci_beta_change_real = c()
  low_ci_beta_change_surrogate =c()
  
  for (i in unique(dataframe$N_channels)){
    
    channel = append(channel, i)
    
    real_beta_1 =  dataframe[dataframe$surrogate == "real" & dataframe$N_channels == i & dataframe$N_conditions == condition, 5]
    surrogate_beta_1 =  dataframe[dataframe$surrogate == "surrogate" & dataframe$N_channels == i & dataframe$N_conditions == condition, 5]
    
    real_beta_2 =  dataframe[dataframe$surrogate == "real" & dataframe$N_channels == i & dataframe$N_conditions == condition, 6]
    surrogate_beta_2 =  dataframe[dataframe$surrogate == "surrogate" & dataframe$N_channels == i & dataframe$N_conditions == condition, 6]
    
    
    real_beta_change =  dataframe[dataframe$surrogate == "real" & dataframe$N_channels == i & dataframe$N_conditions == condition, 7]
    surrogate_beta_change =  dataframe[dataframe$surrogate == "surrogate" & dataframe$N_channels == i & dataframe$N_conditions == condition, 7]
    
    beta_1_cred = beta_1_cred %>% append(
      sum((  real_beta_1 - surrogate_beta_1) > 0)/3000
    )
    beta_1_ER = beta_1_ER %>% append( 
      sum((  real_beta_1 - surrogate_beta_1) > 0)/ # greater than 0 implies that real pairs have *more* initial coherence in channel 1 than surrogate pairs.
        sum(( real_beta_1 - surrogate_beta_1)<= 0))
    
    
    beta_2_cred = beta_2_cred %>% append(
      sum((  real_beta_2 - surrogate_beta_2) > 0)/3000
    )
    beta_2_ER = beta_2_ER %>% append(
      sum((real_beta_2 - surrogate_beta_2) > 0)/ # greater than 0 implies that real pairs have *more* initial coherence in channel 1 than surrogate pairs.
        sum((real_beta_2 - surrogate_beta_2) <= 0))
    
    
    
    beta_change_cred = beta_change_cred %>% append(
      sum((  real_beta_change - surrogate_beta_change) > 0)/3000
    )
    beta_change_ER = beta_change_ER %>% append(
      sum((real_beta_change - surrogate_beta_change) > 0)/ # greater than 0 implies that real pairs have *more* initial coherence in channel 1 than surrogate pairs.
        sum((real_beta_change - surrogate_beta_change) <= 0))
    
    
    map_beta_1_real = append(map_beta_1_real, (map_estimate(real_beta_1))$MAP_Estimate)
    map_beta_1_surrogate = append(map_beta_1_surrogate, (map_estimate(surrogate_beta_1))$MAP_Estimate)
    high_ci_beta_1_real = append(high_ci_beta_1_real, (ci(real_beta_1))$CI_high)
    high_ci_beta_1_surrogate = append(high_ci_beta_1_surrogate, (ci(surrogate_beta_1))$CI_high)
    low_ci_beta_1_real = append(low_ci_beta_1_real, (ci(real_beta_1))$CI_low)
    low_ci_beta_1_surrogate = append(low_ci_beta_1_surrogate, (ci(surrogate_beta_1))$CI_low)
    
    
    
    map_beta_2_surrogate = append(map_beta_2_surrogate, (map_estimate(surrogate_beta_2))$MAP_Estimate)
    map_beta_2_real = append(map_beta_2_real, (map_estimate(real_beta_2))$MAP_Estimate)
    high_ci_beta_2_real = append(high_ci_beta_2_real, (ci(real_beta_2))$CI_high)
    high_ci_beta_2_surrogate = append(high_ci_beta_2_surrogate, (ci(surrogate_beta_2))$CI_high)
    low_ci_beta_2_real = append(low_ci_beta_2_real, (ci(real_beta_2))$CI_low)
    low_ci_beta_2_surrogate = append(low_ci_beta_2_surrogate, (ci(surrogate_beta_2))$CI_low)
    
    
    
    
    map_beta_change_surrogate = append(map_beta_change_surrogate, (map_estimate(surrogate_beta_change))$MAP_Estimate)
    map_beta_change_real = append(map_beta_change_real, (map_estimate(real_beta_change))$MAP_Estimate)
    high_ci_beta_change_real = append(high_ci_beta_change_real, (ci(real_beta_change))$CI_high)
    high_ci_beta_change_surrogate = append(high_ci_beta_change_surrogate, (ci(surrogate_beta_change))$CI_high)
    low_ci_beta_change_real = append(low_ci_beta_change_real, (ci(real_beta_change))$CI_low)
    low_ci_beta_change_surrogate = append(low_ci_beta_change_surrogate, (ci(surrogate_beta_change))$CI_low)
    
    
  }
  
  df_beta_1 = data.frame("parameter" = "beta_1", channel,
                         "map_real" = map_beta_1_real %>% round(3),
                         "CI_low_real" = low_ci_beta_1_real %>% round(3),
                         "CI_high_real" = high_ci_beta_1_real %>% round(3),
                         "map_surrogate" = map_beta_1_surrogate %>% round(3), 
                         "CI_low_surrogate" = low_ci_beta_1_surrogate %>% round(3),
                         "CI_high_surrogate" = high_ci_beta_1_surrogate %>% round(3),
                         "ER" = round(beta_1_ER,3),
                         "credibility" = round(beta_1_cred, 3))
  
  
  df_beta_2 = data.frame("parameter" = "beta_2", channel,
                         "map_real" = map_beta_2_real %>% round(3),
                         "CI_low_real" = low_ci_beta_2_real %>% round(3),
                         "CI_high_real" = high_ci_beta_2_real %>% round(3),
                         "map_surrogate" = map_beta_2_surrogate %>% round(3), 
                         "CI_low_surrogate" = low_ci_beta_2_surrogate %>% round(3),
                         "CI_high_surrogate" = high_ci_beta_2_surrogate %>% round(3),
                         "ER" = round(beta_2_ER,3),
                         "credibility" = round(beta_2_cred, 3))
  
  
  df_beta_change = data.frame("parameter" = "beta_change", channel,
                              "map_real" = map_beta_change_real %>% round(3),
                              "CI_low_real" = low_ci_beta_change_real %>% round(3),
                              "CI_high_real" = high_ci_beta_change_real %>% round(3),
                              "map_surrogate" = map_beta_change_surrogate %>% round(3), 
                              "CI_low_surrogate" = low_ci_beta_change_surrogate %>% round(3),
                              "CI_high_surrogate" = high_ci_beta_change_surrogate %>% round(3),
                              "ER" = round(beta_change_ER,3),
                              "credibility" = round(beta_change_cred, 3))
  
  
  
  
  
  df = rbind(df_beta_1, df_beta_2, df_beta_change)
  return(df)
}



compare_conditions = function(dataframe, surrogate){
  channel = c()
  beta_2_ER = c()
  beta_1_ER = c()
  beta_change_ER = c()
  beta_1_cred = c()
  beta_2_cred = c()
  beta_change_cred = c()
  
  map_beta_1_rest = c()
  map_beta_1_task = c()
  high_ci_beta_1_rest = c()
  high_ci_beta_1_task = c()
  low_ci_beta_1_rest =c()
  low_ci_beta_1_task = c()
  
  map_beta_2_rest = c()
  map_beta_2_task = c()
  high_ci_beta_2_rest = c()
  high_ci_beta_2_task = c()
  low_ci_beta_2_rest = c()
  low_ci_beta_2_task =c()
  
  
  
  map_beta_change_rest = c()
  map_beta_change_task = c()
  high_ci_beta_change_rest = c()
  high_ci_beta_change_task = c()
  low_ci_beta_change_rest = c()
  low_ci_beta_change_task =c()
  
  for (i in unique(dataframe$N_channels)){
    
    channel = append(channel, i)
    
    rest_beta_1 =  dataframe[dataframe$N_conditions == 1 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 5]
    task_beta_1 =  dataframe[dataframe$N_conditions == 2 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 5]
    
    rest_beta_2 =  dataframe[dataframe$N_conditions == 1 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 6]
    task_beta_2 =  dataframe[dataframe$N_conditions == 2 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 6]
    
    
    rest_beta_change =  dataframe[dataframe$N_conditions == 1 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 7]
    task_beta_change =  dataframe[dataframe$N_conditions == 2 & dataframe$N_channels == i & dataframe$surrogate == surrogate, 7]
    
    beta_1_cred = beta_1_cred %>% append(
      sum(( task_beta_1 - rest_beta_1) > 0)/ 3000
    )
    beta_1_ER = beta_1_ER %>% append( 
      sum(( task_beta_1 - rest_beta_1) > 0)/ # greater than 0 implies that rest pairs have *more* initial coherence in channel 1 than task pairs.
        sum(( task_beta_1 - rest_beta_1)<= 0))
    
    beta_2_cred = beta_2_cred %>% append(
      sum(( task_beta_2 - rest_beta_2) > 0)/ 3000
    )
    beta_2_ER = beta_2_ER %>% append(
      sum((task_beta_2 - rest_beta_2) > 0)/ # greater than 0 implies that rest pairs have *more* initial coherence in channel 1 than task pairs.
        sum((task_beta_2 - rest_beta_2) <= 0))
    
    
    
    beta_change_cred = beta_change_cred %>% append(
      sum(( task_beta_change - rest_beta_change) > 0)/ 3000
    )
    beta_change_ER = beta_change_ER %>% append(
      sum((task_beta_change - rest_beta_change) > 0)/ # greater than 0 implies that rest pairs have *more* initial coherence in channel 1 than task pairs.
        sum((task_beta_change - rest_beta_change) <= 0))
    
    
    map_beta_1_rest = append(map_beta_1_rest, (map_estimate(rest_beta_1))$MAP_Estimate)
    map_beta_1_task = append(map_beta_1_task, (map_estimate(task_beta_1))$MAP_Estimate)
    high_ci_beta_1_rest = append(high_ci_beta_1_rest, (ci(rest_beta_1))$CI_high)
    high_ci_beta_1_task = append(high_ci_beta_1_task, (ci(task_beta_1))$CI_high)
    low_ci_beta_1_rest = append(low_ci_beta_1_rest, (ci(rest_beta_1))$CI_low)
    low_ci_beta_1_task = append(low_ci_beta_1_task, (ci(task_beta_1))$CI_low)
    
    
    
    map_beta_2_task = append(map_beta_2_task, (map_estimate(task_beta_2))$MAP_Estimate)
    map_beta_2_rest = append(map_beta_2_rest, (map_estimate(rest_beta_2))$MAP_Estimate)
    high_ci_beta_2_rest = append(high_ci_beta_2_rest, (ci(rest_beta_2))$CI_high)
    high_ci_beta_2_task = append(high_ci_beta_2_task, (ci(task_beta_2))$CI_high)
    low_ci_beta_2_rest = append(low_ci_beta_2_rest, (ci(rest_beta_2))$CI_low)
    low_ci_beta_2_task = append(low_ci_beta_2_task, (ci(task_beta_2))$CI_low)
    
    
    
    
    map_beta_change_task = append(map_beta_change_task, (map_estimate(task_beta_change))$MAP_Estimate)
    map_beta_change_rest = append(map_beta_change_rest, (map_estimate(rest_beta_change))$MAP_Estimate)
    high_ci_beta_change_rest = append(high_ci_beta_change_rest, (ci(rest_beta_change))$CI_high)
    high_ci_beta_change_task = append(high_ci_beta_change_task, (ci(task_beta_change))$CI_high)
    low_ci_beta_change_rest = append(low_ci_beta_change_rest, (ci(rest_beta_change))$CI_low)
    low_ci_beta_change_task = append(low_ci_beta_change_task, (ci(task_beta_change))$CI_low)
    
    
  }
  
  df_beta_1 = data.frame("parameter" = "beta_1", channel,
                         "map_rest" = map_beta_1_rest %>% round(3),
                         "CI_low_rest" = low_ci_beta_1_rest %>% round(3),
                         "CI_high_rest" = high_ci_beta_1_rest %>% round(3),
                         "map_task" = map_beta_1_task %>% round(3), 
                         "CI_low_task" = low_ci_beta_1_task %>% round(3),
                         "CI_high_task" = high_ci_beta_1_task %>% round(3),
                         "ER" = round(beta_1_ER,3),
                         "credibility" = round(beta_1_cred, 3))
  
  
  df_beta_2 = data.frame("parameter" = "beta_2", channel,
                         "map_rest" = map_beta_2_rest %>% round(3),
                         "CI_low_rest" = low_ci_beta_2_rest %>% round(3),
                         "CI_high_rest" = high_ci_beta_2_rest %>% round(3),
                         "map_task" = map_beta_2_task %>% round(3), 
                         "CI_low_task" = low_ci_beta_2_task %>% round(3),
                         "CI_high_task" = high_ci_beta_2_task %>% round(3),
                         "ER" = round(beta_2_ER,3),
                         "credibility" = round(beta_2_cred, 3))
  
  
  df_beta_change = data.frame("parameter" = "beta_change", channel,
                              "map_rest" = map_beta_change_rest %>% round(3),
                              "CI_low_rest" = low_ci_beta_change_rest %>% round(3),
                              "CI_high_rest" = high_ci_beta_change_rest %>% round(3),
                              "map_task" = map_beta_change_task %>% round(3), 
                              "CI_low_task" = low_ci_beta_change_task %>% round(3),
                              "CI_high_task" = high_ci_beta_change_task %>% round(3),
                              "ER" = round(beta_change_ER,3),
                              "credibility" = round(beta_change_cred, 3))
  
  
  
  
  
  df = rbind(df_beta_1, df_beta_2, df_beta_change)
  return(df)
}


################### PARAMETER RECOVERY ###############


parameter_recovery = function(simulated_data_array, parameters, stan_model, chains, parallel_chains, iter_warmup, iter_sampling, filenames, directory = "ParameterRecovery", refresh = 15, seed = 1337){
  
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

save_summary = function(directory, filename, number){
  
  data_path = paste(directory, "/", filename, "_", number, ".RDS", sep = "")
  out_path = paste(directory, "/", "sim", number, "_summary", ".RDS", sep = "")
  sim_summary = readRDS( data_path )
  sim_summary = sim_summary$summary(.cores = 4)
  write_csv(sim_summary, out_path)
  
}
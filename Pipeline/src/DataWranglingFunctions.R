require(pacman)
pacman::p_load(biwavelet, tidyverse, pacman)
suppressPackageStartupMessages(library(tidyverse))


index_similar_value = function(vector, value){
  
  new_vector = abs(vector - value)
  
  indeces = which(new_vector == min(new_vector))
  
  if (length(indeces) > 1){
    
    return(indeces[length(indeces)])
  }
  
  
  
  return(indeces)
}


coi_to_indeces = function(coi, period){
  
  coi$period_idx = NA
  coi$times_idx = NA
  
  for (n in 1:nrow(coi)){
    coi$period_idx[n] = index_similar_value(period, coi$X1[n])
    coi$times_idx[n] = n
    
  }
  return(coi)
}




remove_coi_from_data = function(coherence_data, formatted_coi){
  
  for (n in 1:nrow(formatted_coi)){
    column_index = formatted_coi$times_idx[n]
    row_index = formatted_coi$period_idx[n]
    
    coherence_data[row_index:nrow(coherence_data), column_index] = NA
    
    
  }
  
  return(coherence_data)
}



slice_dataframe = function(dataframe, period, times, trigger, duration, max_period, min_period){
  
  row_min = index_similar_value(period, min_period)
  row_max = index_similar_value(period, max_period)
  column_min = index_similar_value(times, trigger)
  column_max = index_similar_value(times, trigger+duration)
  
  
  return(dataframe[row_min:row_max, column_min:column_max])
}




collapse_and_average_dataframe = function(dataframe){
  
  vec = c()
  
  for (n in 1:nrow(dataframe)){
    vec = c(vec, dataframe[n,])
  }
  
  vec <- as.numeric(vec[!is.na(vec)])
  return(mean(vec))
}

get_surrogate_pair_triggers = function(root_dir, visit, pair, SurrogatePairList){
  
  participant1 = SurrogatePairList[SurrogatePairList$SurrogatePairNumber == pair, 2][1] %>% as.numeric()
  
  participant2 = SurrogatePairList[SurrogatePairList$SurrogatePairNumber == pair, 3][1] %>% as.numeric()
  
  triggers_participant_1 = paste(
    root_dir,
    "/",
    participant1,
    "/",
    visit, "/_triggers.tsv", sep = "") %>% read_tsv(col_names = FALSE, show_col_types = FALSE)
  
  triggers_participant_2= paste(
    root_dir,
    "/",
    participant2,
    "/",
    visit, "/_triggers.tsv", sep = "") %>% read_tsv(col_names = FALSE, show_col_types = FALSE)
  
  surrogate_triggers = c()
  
  for (n in 1:length(triggers_participant_1$X1)){
    t1 = triggers_participant_1$X1[n]
    t2 = triggers_participant_2$X1[n]
    surrogate_triggers = c(surrogate_triggers, mean(c(t1, t2)))
    
  }
  return(surrogate_triggers)
}



channel_task_and_rest_averages= function(root_dir, pair, visit, channel, Surrogate = FALSE, SurrogatePairList = NA){
  
  # Construct filepath
  filepath = paste(
    root_dir,
    "/HaemoglobinTimeSeries",
    "/",
    pair,
    "/",
    visit, sep = "")
  if (Surrogate == TRUE){
    triggers = get_surrogate_pair_triggers(root_dir, visit, pair, SurrogatePairList)
  }
  else {
    triggers = (filepath %>% paste("/_triggers.tsv", sep = "") %>% read_tsv(col_names = FALSE, show_col_types = FALSE))$X1
  }
  
  filepath = filepath %>% paste("/SynchronyResults",
                                "/file_", pair, "_", channel,
                                sep = ""  )
  
  rsq = paste(filepath, "_Rsq.csv", sep = "") %>% read_csv(col_names = FALSE, show_col_types = FALSE)
  period = paste(filepath, "_period.csv", sep = "") %>% read_csv(col_names = FALSE, show_col_types = FALSE)
  times = paste(filepath, "_time.csv", sep = "") %>% read_csv(col_names = FALSE, show_col_types = FALSE)
  coi = paste(filepath, "_coi.csv", sep = "") %>% read_csv(col_names = FALSE, show_col_types = FALSE)
  
  # format coi
  
  formatted_coi = coi_to_indeces(coi, period)
  rsq_minus_coi = remove_coi_from_data(rsq, formatted_coi)
  
  rest1_avg = slice_dataframe(rsq_minus_coi,
                              period$X1, times$X1, triggers[1], duration = 120, max_period = 20, min_period = 10) %>% collapse_and_average_dataframe()
  
  task1_avg= slice_dataframe(rsq_minus_coi,
                             period$X1, times$X1, triggers[2], duration = 60*5, max_period = 20, min_period = 10) %>% collapse_and_average_dataframe()
  
  rest2_avg = slice_dataframe(rsq_minus_coi,
                              period$X1, times$X1, triggers[3], duration = 60, max_period = 20, min_period = 10) %>% collapse_and_average_dataframe()
  
  
  task2_avg = slice_dataframe(rsq_minus_coi,
                              period$X1, times$X1, triggers[4], duration = 60*5, max_period = 20, min_period = 10) %>% collapse_and_average_dataframe()
  
  rest3_avg = slice_dataframe(rsq_minus_coi,
                              period$X1, times$X1, triggers[5], duration = 120, max_period = 20, min_period = 10) %>% collapse_and_average_dataframe()
  
  
  rest_avg = mean(rest1_avg, rest2_avg, rest3_avg)
  task_avg = mean(task1_avg, task2_avg)
  
  results = c(rest_avg, task_avg) # REST IS FIRST!
  return(results)
}


get_all_channel_results = function(channel_names, pair, visit, Surrogate = FALSE, SurrogatePairList = NA){
  root_dir = getwd()
  
  for (channel in channel_names){
    
    
    tryCatch({
      
      results = channel_task_and_rest_averages(root_dir, pair, visit, channel, Surrogate, SurrogatePairList)
      
      
      if (exists("total_results") == TRUE){
        more_results = data.frame("channel" = channel,"rest_avg" = results[1], "task_avg" = results[2])
        
        total_results = rbind(total_results, more_results)
        
      }
      else{
        total_results = data.frame("channel" = channel,"rest_avg" = results[1], "task_avg" = results[2])
        
      }
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    })
  }
  if(exists("total_results") == TRUE){
    return(total_results)
  }
}


analyse_both_visits = function(channel_names, pair, Surrogate = FALSE, SurrogatePairList = NA){
  
  
  visit1_df = get_all_channel_results(channel_names, pair, visit = "Visit1", Surrogate, SurrogatePairList)
  visit1_df$visit = 1
  
  
  visit4_df = get_all_channel_results(channel_names, pair, visit = "Visit4", Surrogate, SurrogatePairList)
  visit4_df$visit = 4
  
  
  return(rbind(visit1_df, visit4_df))
  
  
}


loop_over_all_pairs = function(pair_list, channel_names, filename, Surrogate = FALSE, SurrogatePairList = NA){
  
  for (pair in pair_list){
    
    pair_df = analyse_both_visits(CHANNEL_NAMES, pair, Surrogate = FALSE, SurrogatePairList = NA)
    pair_df$pair = pair
    
    print(pair)
    
    if (exists("all_pairs") == TRUE){
      all_pairs = rbind(all_pairs, pair_df)
      
    }
    else{
      all_pairs = pair_df
      
    }
    write_csv(all_pairs, filename)
    
  }
  write_csv(all_pairs, filename)
}




dummy_code = function(vector){
  
  
  new_levels = c(1:length(unique(vector)))
  
  levels(vector) <- new_levels
  
  return(vector)
}
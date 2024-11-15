require(pacman)
pacman::p_load(biwavelet, tidyverse, pacman)
suppressPackageStartupMessages(library(tidyverse))

# Load surrogate pairs
surrogate_pairs = read_csv("temp/SurrogatePairs.csv")


get_truncation_length = function(ts1, ts2){
  # Takes two timeseries (in dataframe format) and ensures they're of
  # equal length. Subsequently, takes the timepoint vector of the first time series and uses
  # it for the second one.
  
  length1 = length(ts1[,1])
  length2 = length(ts2[,1])
  
  desired_length = min(length1, length2)
  
  return(desired_length)
  
  
}


save_csvs = function(analysis, output_path, filename){
  # rsq
  data.frame(analysis["rsq"]) %>% 
  write_csv(paste(output_path, "/", filename, "_rsq.csv", sep = ""), col_names = FALSE)
  #coi
  data.frame(analysis["coi"]) %>% 
  write_csv(paste(output_path, "/", filename, "_coi.csv", sep = ""), col_names = FALSE)
  #phase
  data.frame(analysis["phase"]) %>% 
  write_csv(paste(output_path, "/", filename, "_phase.csv", sep = ""), col_names = FALSE)
  #period
  data.frame(analysis["period"]) %>% 
  write_csv(paste(output_path, "/", filename, "_period.csv", sep = ""), col_names = FALSE)
  #time
  data.frame(analysis["t"]) %>% 
  write_csv(paste(output_path, "/", filename, "_time.csv", sep = ""), col_names = FALSE)
  # signif = data.frame(analysis["signif])  
  
}

save_plot = function(analysis, output_path, filename){
  png(paste(output_path, "/", "plot", "_", filename, ".png", sep = ""))
  par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
  plot(analysis, plot.cb=TRUE)
  dev.off()
  
}

single_channel_coherence = function(input_path1, input_path2){
  
  d1 = data.frame(read_tsv(input_path1, col_types = cols()))
  d2 = data.frame(read_tsv(input_path2, col_types = cols()))
  
  truncation_length = get_truncation_length(d1, d2)
  
  d1 = d1[1:truncation_length,]
  d2 = d2[1:truncation_length,]
  
  # Use common time column
  
  d2[,1] = d1[,1]
  d1 = cbind(round(d1[,1],4), d1[,2])
  d2 = cbind(round(d2[,1],4), d2[,2])
  tryCatch(
    #try to do this
    {
      
      analysis = wtc(
        d1,
        d2,
        nrands = 0,
        quiet = TRUE
      )
      
      # rsq = data.frame(analysis["rsq"])
      
      return(analysis)
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
    })
  
  # rsq = data.frame(analysis["rsq"])
  
  return(analysis)
  
}



create_directory = function(new_path){
# check if sub directory exists 
if (file.exists(new_path) == FALSE){
  # create a new sub directory inside
  # the main path
  dir.create(file.path(getwd(), new_path))
  return(TRUE)
}
  return(FALSE)}



# Create fake directory list here.


#############################################################


analyse_visit_surrogate = function(surrogate_pair, input_path1, input_path2, results_path,  channel_names){
  
  
  for (channel_name in channel_names){
    results_filename_root = paste("file_" , as.character(surrogate_pair) , "_" , as.character(channel_name), sep = "")
    
    channel_path_1 = paste(input_path1, '/1_' , channel_name , '.tsv', sep = "")
    channel_path_2 = paste(input_path2, '/2_' , channel_name , '.tsv', sep = "")
    tryCatch(
      #try to do this
      {
        
        
        analysis = single_channel_coherence(channel_path_1, channel_path_2)
        
        save_csvs(analysis, results_path, results_filename_root)
        save_plot(analysis, results_path, results_filename_root)
        
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
      })
  }
  return()
}


analyse_pair_surrogate = function(surrogate_pair, participant1_nr, participant2_nr, channel_names){
  
  
  input_path1 = paste(participant1_nr, "/Visit1")
  
  # Visit 1
  tryCatch(
    #try to do this
    {
      
      results_path = paste("HaemoglobinTimeSeries/", surrogate_pair, "/Visit1", "/SynchronyResults", sep = "")
      create_directory(results_path)
      input_path1 = paste("HaemoglobinTimeSeries/", participant1_nr, "/Visit1", sep = "")
      input_path2 = paste("HaemoglobinTimeSeries/", participant2_nr, "/Visit1", sep = "")
      analyse_visit_surrogate(surrogate_pair, input_path1, input_path2, results_path,  channel_names)
      
      
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
    })
  
  # Visit 4
  results_path = paste("HaemoglobinTimeSeries/", surrogate_pair, "/Visit4", "/SynchronyResults", sep = "")
  
  tryCatch(
    #try to do this
    {
  create_directory(results_path)
  input_path1 = paste("HaemoglobinTimeSeries/", participant1_nr, "/Visit4", sep = "")
  input_path2 = paste("HaemoglobinTimeSeries/", participant2_nr, "/Visit4", sep = "")
  analyse_visit_surrogate(surrogate_pair, input_path1, input_path2, results_path,  channel_names)
  
  
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
  })
  
  return()
}


CHANNEL_NAMES =  c(
  'S1_D1 hbo',
  'S1_D2 hbo',
  'S1_D3 hbo',
  'S2_D1 hbo',
  'S2_D4 hbo',
  'S3_D1 hbo',
  'S3_D3 hbo',
  'S3_D4 hbo',
  'S4_D4 hbo',
  'S4_D5 hbo',
  'S4_D10 hbo',
  'S5_D2 hbo',
  'S5_D3 hbo',
  'S5_D7 hbo',
  'S6_D2 hbo',
  'S6_D7 hbo',
  'S7_D3 hbo',
  'S7_D5 hbo',
  'S7_D6 hbo',
  'S7_D8 hbo',
  'S8_D6 hbo',
  'S8_D7 hbo',
  'S8_D16 hbo',
  'S9_D4 hbo',
  'S9_D9 hbo',
  'S10_D4 hbo',
  'S10_D9 hbo',
  'S10_D10 hbo',
  'S10_D11 hbo',
  'S11_D5 hbo',
  'S11_D8 hbo',
  'S11_D10 hbo',
  'S12_D9 hbo',
  'S12_D11 hbo',
  'S12_D12 hbo',
  'S13_D10 hbo',
  'S13_D11 hbo',
  'S13_D13 hbo',
  'S14_D11 hbo',
  'S14_D12 hbo',
  'S14_D13 hbo',
  'S14_D14 hbo',
  'S15_D12 hbo',
  'S15_D14 hbo',
  'S15_D15 hbo',
  'S16_D13 hbo',
  'S16_D14 hbo',
  'S17_D7 hbo',
  'S17_D17 hbo',
  'S18_D6 hbo',
  'S18_D8 hbo',
  'S18_D16 hbo',
  'S19_D7 hbo',
  'S19_D16 hbo',
  'S19_D17 hbo',
  'S19_D18 hbo',
  'S20_D16 hbo',
  'S20_D18 hbo',
  'S20_D19 hbo',
  'S21_D17 hbo',
  'S21_D18 hbo',
  'S21_D20 hbo',
  'S22_D18 hbo',
  'S22_D19 hbo',
  'S22_D20 hbo',
  'S22_D21 hbo',
  'S23_D19 hbo',
  'S23_D21 hbo',
  'S24_D20 hbo',
  'S24_D21 hbo',
  'S24_D22 hbo',
  'S1_D1 hbr',
  'S1_D2 hbr',
  'S1_D3 hbr',
  'S2_D1 hbr',
  'S2_D4 hbr',
  'S3_D1 hbr',
  'S3_D3 hbr',
  'S3_D4 hbr',
  'S4_D4 hbr',
  'S4_D5 hbr',
  'S4_D10 hbr',
  'S5_D2 hbr',
  'S5_D3 hbr',
  'S5_D7 hbr',
  'S6_D2 hbr',
  'S6_D7 hbr',
  'S7_D3 hbr',
  'S7_D5 hbr',
  'S7_D6 hbr',
  'S7_D8 hbr',
  'S8_D6 hbr',
  'S8_D7 hbr',
  'S8_D16 hbr',
  'S9_D4 hbr',
  'S9_D9 hbr',
  'S10_D4 hbr',
  'S10_D9 hbr',
  'S10_D10 hbr',
  'S10_D11 hbr',
  'S11_D5 hbr',
  'S11_D8 hbr',
  'S11_D10 hbr',
  'S12_D9 hbr',
  'S12_D11 hbr',
  'S12_D12 hbr',
  'S13_D10 hbr',
  'S13_D11 hbr',
  'S13_D13 hbr',
  'S14_D11 hbr',
  'S14_D12 hbr',
  'S14_D13 hbr',
  'S14_D14 hbr',
  'S15_D12 hbr',
  'S15_D14 hbr',
  'S15_D15 hbr',
  'S16_D13 hbr',
  'S16_D14 hbr',
  'S17_D7 hbr',
  'S17_D17 hbr',
  'S18_D6 hbr',
  'S18_D8 hbr',
  'S18_D16 hbr',
  'S19_D7 hbr',
  'S19_D16 hbr',
  'S19_D17 hbr',
  'S19_D18 hbr',
  'S20_D16 hbr',
  'S20_D18 hbr',
  'S20_D19 hbr',
  'S21_D17 hbr',
  'S21_D18 hbr',
  'S21_D20 hbr',
  'S22_D18 hbr',
  'S22_D19 hbr',
  'S22_D20 hbr',
  'S22_D21 hbr',
  'S23_D19 hbr',
  'S23_D21 hbr',
  'S24_D20 hbr',
  'S24_D21 hbr',
  'S24_D22 hbr')


for (i in 1:nrow(surrogate_pairs)){
  analyse_pair_surrogate(surrogate_pair = surrogate_pairs[i,1],
                         participant1_nr = surrogate_pairs[i,2],
                         participant2_nr = surrogate_pairs[i,3],
                         channel_names = CHANNEL_NAMES)
}



# Pairs 3 and 8 have ARIMA errors. The trycatch function in the code should ensure these errors
# don't disturb the process. Note that this means missing data later!


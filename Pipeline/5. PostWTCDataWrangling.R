require(pacman)
pacman::p_load(biwavelet, tidyverse, pacman)
suppressPackageStartupMessages(library(tidyverse))
source("src/DataWranglingFunctions.R")


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
  "S10_D11",
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
  "S22_D21 hbr",
  'S23_D19 hbr',
  "S23_D21 hbr",
  'S24_D20 hbr',
  'S24_D21 hbr',
  'S24_D22 hbr')

# REAL PAIRS

pair_list = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23)
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
  write_csv(all_pairs, "temp/RealPairsTemp.csv")
  
}
write_csv(all_pairs, "temp/RealPairsTemp.csv")


# SURROGATE PAIRS

SurrogatePairList = read_csv("SurrogatePairs.csv", show_col_types = FALSE)

pair_list = c(
              "Surrogate1",
              "Surrogate2",
              "Surrogate3",
              "Surrogate4",
              "Surrogate5",
              "Surrogate6",
              "Surrogate7",
              "Surrogate8",
              "Surrogate9",
              "Surrogate10",
              "Surrogate11",
              "Surrogate12",
              "Surrogate13",
              "Surrogate14",
              "Surrogate15",
              "Surrogate16",
              "Surrogate17",
              "Surrogate18",
              "Surrogate19",
              "Surrogate20",
              "Surrogate21",
              "Surrogate22")



for (pair in pair_list){
  
  pair_df = analyse_both_visits(CHANNEL_NAMES, pair, Surrogate = TRUE, SurrogatePairList = SurrogatePairList)
  pair_df$pair = pair
  
  print(pair)
  
  if (exists("all_pairs") == TRUE){
    all_pairs = rbind(all_pairs, pair_df)
    
  }
  else{
    all_pairs = pair_df
    
  }
  write_csv(all_pairs, "SurrogatePairsTemp.csv")
  
}
write_csv(all_pairs, "SurrogatePairsTemp.csv")



real_pairs_processed = read_csv("RealPairsTemp.csv")

surrogate_pairs_processed = read_csv("SurrogatePairsTemp.csv")

###################################################
######### ADDING COLUMNS AND FORMATTING ###########
###################################################


real_pairs_processed$surrogate = "real"
surrogate_pairs_processed$surrogate = "surrogate"
processed_pairs = rbind(real_pairs_processed, surrogate_pairs_processed)
write_csv(processed_pairs, "temp/processed_pairs.csv")


processed_pairs = read_csv("temp/processed_pairs.csv")

all_channels = unique(processed_pairs$channel)

central_channels = c(
  "S1_D1",
  "S1_D2",
  "S1_D3",
  "S3_D3",
  "S5_D3",
  "S7_D3",
  "S7_D5",
  "S7_D6",
  "S7_D8",
  "S11_D8",
  "S18_D8"
)
right_channels = c(
  "S5_D2",
  "S6_D2",
  "S5_D7",
  "S6_D7",
  "S8_D7",
  "S8_D16",
  "S17_D7",
  "S17_D17",
  "S19_D7",
  "S19_D16",
  "S19_D17",
  "S19_D18",
  "S20_D16",
  "S20_D18",
  "S20_D19",
  "S21_D17",
  "S21_D18",
  "S21_D20",
  "S22_D18",
  "S22_D19",
  "S22_D20",
  "S22_D21",
  "S23_D19",
  "S23_D21",
  "S24_D20",
  "S24_D21",
  "S24_D22"
)

right_channels_hbo = paste(right_channels, " hbo", sep = "")
all(right_channels_hbo %in% all_channels)
right_channels_hbr = paste(right_channels, " hbr", sep = "")
all(right_channels_hbr %in% all_channels)
right_channels = c(right_channels_hbo, right_channels_hbr)

central_channels_hbo = paste(central_channels, " hbo", sep = "")
all(central_channels_hbo %in% all_channels)
central_channels_hbr = paste(central_channels, " hbr", sep = "")
all(central_channels_hbr %in% all_channels)
central_channels = c(central_channels_hbo, central_channels_hbr)


Fp_channels = c(
  "S1_D1",
  "S1_D2",
  "S1_D3",
  "S2_D1",
  "S3_D1"
  ,"S5_D2",
  "S6_D2"
)
Fp_channels_hbo = paste(Fp_channels, " hbo", sep = "")
all(Fp_channels_hbo %in% all_channels)
Fp_channels_hbr = paste(Fp_channels, " hbr", sep = "")
all(Fp_channels_hbr %in% all_channels)
Fp_channels = c(Fp_channels_hbo, Fp_channels_hbr)

Af_channels = c(
  "S2_D4",
  "S3_D4",
  "S3_D3",
  "S5_D3",
  "S5_D7",
  "S6_D7",
  "S7_D3"
)
Af_channels_hbo = paste(Af_channels, " hbo", sep = "")
all(Af_channels_hbo %in% all_channels)
Af_channels_hbr = paste(Af_channels, " hbr", sep = "")
all(Af_channels_hbr %in% all_channels)
Af_channels = c(Af_channels_hbo, Af_channels_hbr)

F_channels = c(
  "S9_D4",
  "S9_D9"
  ,"S4_D4",
  "S4_D5",
  "S4_D10",
  "S7_D5",
  "S7_D6",
  "S7_D8",
  "S8_D6",
  "S8_D7",
  "S8_D16",
  "S10_D4",
  "S11_D5",
  "S17_D7",
  "S17_D17",
  "S18_D6",
  "S19_D7"
)
F_channels_hbo = paste(F_channels, " hbo", sep = "")
all(F_channels_hbo %in% all_channels)
F_channels_hbr = paste(F_channels, " hbr", sep = "")
all(F_channels_hbr %in% all_channels)
F_channels = c(F_channels_hbo, F_channels_hbr)

FC_T_channels = c(
  "S10_D9",
  "S10_D10",
  "S10_D11",
  "S11_D10",
  "S11_D8",
  "S18_D8",
  "S18_D16",
  "S19_D16",
  "S19_D17",
  "S19_D18",
  "S12_D9",
  "S13_D10",
  "S20_D16",
  "S21_D17"
)
FC_T_channels_hbo = paste(FC_T_channels, " hbo", sep = "")
all(FC_T_channels_hbo %in% all_channels)
FC_T_channels_hbr = paste(FC_T_channels, " hbr", sep = "")
all(FC_T_channels_hbr %in% all_channels)
FC_T_channels = c(FC_T_channels_hbo, FC_T_channels_hbr)

C_T_channels = c(
  "S12_D11",
  "S13_D11",
  "S20_D18",
  "S21_D18",
  "S22_D18",
  "S20_D19",
  "S21_D20",
  "S12_D12",
  "S13_D13",
  "S14_D11"
)
C_T_channels_hbo = paste(C_T_channels, " hbo", sep = "")
all(C_T_channels_hbo %in% all_channels)
C_T_channels_hbr = paste(C_T_channels, " hbr", sep = "")
all(C_T_channels_hbr %in% all_channels)
C_T_channels = c(C_T_channels_hbo, C_T_channels_hbr)

C_TP_channels = c(
  "S14_D12",
  "S14_D13",
  "S14_D14",
  "S15_D12",
  "S16_D13",
  "S22_D19",
  "S22_D20",
  "S22_D21",
  "S23_D19",
  "S24_D20"
)
C_TP_channels_hbo = paste(C_TP_channels, " hbo", sep = "")
all(C_TP_channels_hbo %in% all_channels)
C_TP_channels_hbr = paste(C_TP_channels, " hbr", sep = "")
all(C_TP_channels_hbr %in% all_channels)
C_TP_channels = c(C_TP_channels_hbo, C_TP_channels_hbr)


P_channels = c(
  "S15_D14",
  "S16_D14",
  "S15_D15",
  "S23_D21",
  "S24_D21",
  "S24_D22"
)
P_channels_hbo = paste(P_channels, " hbo", sep = "")
all(P_channels_hbo %in% all_channels)
P_channels_hbr = paste(P_channels, " hbr", sep = "")
all(P_channels_hbr %in% all_channels)
P_channels = c(P_channels_hbo, P_channels_hbr)


###################################################
processed_pairs$hemisphere = NA
processed_pairs$cap_region = NA

uncategorized_cap_region = 0
for (row in 1:nrow(processed_pairs)){
  channel = processed_pairs$channel[row]
  
  if (channel %in% central_channels){
    processed_pairs$hemisphere[row] = "Central"
  }
  else if (channel %in% right_channels){
    processed_pairs$hemisphere[row] = "Right"}
  else{
    processed_pairs$hemisphere[row] = "Left"
    
  }
  if (channel %in% Fp_channels){
    
    processed_pairs$cap_region[row] = "Fp"
  }
  else if (channel %in% Af_channels){
    processed_pairs$cap_region[row] = "Af"}
  else if (channel %in% F_channels){
    processed_pairs$cap_region[row] = "F"}
  else if (channel %in% FC_T_channels){
    processed_pairs$cap_region[row] = "FC_T"}
  else if (channel %in% C_T_channels){
    processed_pairs$cap_region[row] = "C_T"}
  else if (channel %in% C_TP_channels){
    processed_pairs$cap_region[row] = "C_TP"}
  else if (channel %in% P_channels){
    processed_pairs$cap_region[row] = "P"}
  else{
    uncategorized_cap_region = uncategorized_cap_region +1 }
  

  
}

write_csv(processed_pairs, "temp/AverageCoherenceData.csv")


# Pivot_longer

AverageCoherenceData = read_csv("temp/AverageCoherenceData.csv")


AverageCoherenceDataLong =  pivot_longer(AverageCoherenceData,
                                         
                                         cols = c(rest_avg, task_avg), 
                                         names_to = "condition", 
                                         values_to = "coherence"
)



write_csv(AverageCoherenceDataLong, "temp/CoherenceLong.csv")
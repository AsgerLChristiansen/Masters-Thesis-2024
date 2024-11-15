require(pacman)
pacman::p_load(tidyverse, pacman)
suppressPackageStartupMessages(library(tidyverse))


set.seed(1337)

# Exclude pair 15 due to bad data
participant1 = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23)


participant2 = sample(participant1)

# Sanity check
all(participant1 != participant2)

# turn into dataframe

surrogate_pairs = data.frame("SurrogatePairNumber" = paste('Surrogate', c(1:22), sep = ""), "participant1" = participant1, "participant2" = participant2)
write_csv(surrogate_pairs,"temp/SurrogatePairs.csv",  col_names = TRUE)

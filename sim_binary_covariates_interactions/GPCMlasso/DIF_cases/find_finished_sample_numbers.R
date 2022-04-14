## After simulations copy each result file from "result_files" 
## into "result_files - copied from all tries"

## Then find codes of samples finished - and thereby of samples not finished - 
## with this code.

## Finally conduct simulations for the samples not finished (maybe with fewer
## cores for better RAM availability) and again copy results from 
## "result_files" to "result_files - copied from all tries"


# how to:
# read file names
# separat seed in the beginning
# put seed into numeric vector
# change seeds to number of sample (-1998)
# invert numbers to missing numbers in 1:10000



filenames <- list.files("result_files_all/")
head(filenames)

seeds_done <- unlist(lapply(filenames, function(x) as.numeric((strsplit(x,"_"))[[1]][1])))
seeds_done_sorted <- sort(seeds_done)

numbers_of_samples_done <- seeds_done_sorted - 1998

numbers_of_missing_samples <- (1:10000)[-numbers_of_samples_done]

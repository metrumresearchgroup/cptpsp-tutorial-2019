## This script reproduces the results for the mrgsolve introduction section in the manuscript. 
## Note: Before proceeding, make "script" your working directory.

# clear workspace
rm(list=ls())
gc()

# load libraries
.libPaths("lib")
library(tidyverse)
library(mrgsolve)

################################################################################################################
#################################################### Chunk 1 ###################################################
################################################################################################################

# compile pk1 model from model library modlib()
mod <- mread("pk1", "../model")

# show model info
mod

# look at model parameters
param(mod)

# create event object for a daily extravascular dose of 100 given for 10 days 
evnt <- ev(amt = 100, ii = 24, addl = 9)

# simulate till time == 480 with an output every 0.1 time units
out <- mod %>% ev(evnt) %>% mrgsim(end = 480, delta = 0.1)

# glimpse simulation output
head(out)

# plot
plot(out)

################################################################################################################
################################################################################################################


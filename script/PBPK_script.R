## This script reproduces the results from manuscript. 
## Note: Before proceeding, make "script" your working directory.

# clear workspace
rm(list=ls())
gc()

# load libraries
.libPaths("lib")
library(tidyverse)
library(cowplot)
library(mrgsolve)

##################################################################################################
##################################################################################################

## load observed data
obsA <- read.csv("../data/Adult_IV.csv")  #adult data
obsP <- read.csv("../data/Pediatric_IV.csv")  #adult data

## generate models
modA <- mread("../model/voriPBPK")






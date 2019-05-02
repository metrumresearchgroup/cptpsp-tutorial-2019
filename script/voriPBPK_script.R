## This script reproduces the results for voriconazole PBPK section in the manuscript. 
## Note: Before proceeding, make "script" your working directory.

# clear workspace
rm(list=ls())
gc()

# load libraries
.libPaths("lib")
library(tidyverse)
library(cowplot)
library(mrgsolve)

################################################################################################################
#################################################### Chunk 1 ###################################################
################################################################################################################

## compile and generate adult and pediatric PBPK models
modA <- mread("../model/voriPBPK")

# pediatric (5 yo) male physiology; https://www.ncbi.nlm.nih.gov/pubmed/14506981
pedPhys <- list(WEIGHT = 19,
                Vad = 5.5,
                Vbo = 2.43,
                Vbr = 1.31,
                VguWall = 0.22,
                VguLumen = 0.117,
                Vhe = 0.085,
                Vki = 0.11,
                Vli = 0.467,
                Vlu = 0.125,
                Vmu = 5.6,
                Vsp = 0.05,
                Vbl = 1.5,
                Qad = 0.05*3.4*60,
                Qbo = 0.05*3.4*60,
                Qbr = 0.12*3.4*60,
                Qgu = 0.15*3.4*60, 
                Qhe = 0.04*3.4*60,
                Qki = 0.19*3.4*60,
                Qmu = 0.17*3.4*60,
                Qsp = 0.03*3.4*60,
                Qha = 0.065*3.4*60, 
                Qlu = 3.4*60,
                MPPGL = 26,
                VmaxH = 120.5,
                KmH = 11)

modP <- param(modA, pedPhys)  #generate pediatric model by updating adult model with pediatric physiology

################################################################################################################
################################################################################################################

################################################################################################################
#################################################### Chunk 2 ###################################################
################################################################################################################

## get partition coefficients using Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
### voriconazole physicochemical properties
tissueComp <- read.csv("../data/tissue_comp_PT.csv")  # tissue composition
source("calcKp_PT.R")
Kp <- calcKp_PT(logP=2.56, pKa=1.76, fup=0.42, BP=1, type=3, dat=tissueComp)

#update model parameters partition coefficients
modA <- param(modA, Kp)
modP <- param(modP, Kp)

################################################################################################################
################################################################################################################

################################################################################################################
#################################################### Chunk 3 ###################################################
################################################################################################################

## Run simulations and compare the predicted plasma concentration to the observed data

### Adult: simulate a steady state plasma concentration for a voriconazole dose of 4 mg/kg IV dose 
### infused over an hour given twice a day for a week to a 73 kg adult male.
bw <- 73
simA <- 
  modA %>% 
  param(WEIGHT = bw) %>%
  ev(cmt="VEN", amt=4*bw, rate=4*bw, ii=12, addl=13, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)


### Pediatric: simulate a steady state plasma concentration for a voriconazole dose of 4 mg/kg IV dose 
### infused with a rate of 3 mg/kg given twice a day for a week to a 19 kg male child.
bw <- 19
simP <- 
  modP %>% 
  param(WEIGHT = bw) %>%
  ev(cmt="VEN", amt=4*bw, rate=3*bw, ii=12, addl=13, ss=1) %>% 
  mrgsim(delta = 0.1, end = 12) %>% 
  filter(row_number() != 1)

################################################################################################################
################################################################################################################

################################################################################################################
#################################################### Chunk 4 ###################################################
################################################################################################################

## Validate model predictions against observed data

### Adult
#### load observed data; digitized using webplotdigitizer https://automeris.io/WebPlotDigitizer/
obsA <- read.csv("../data/Adult_IV.csv")  #adult IV infusion data; dose and rate = 4 mg/kg

### plot
gp1 <- 
  ggplot() + 
  geom_point(data = obsA, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obsA, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = simA, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Adult 4 mg/kg IV", x="", y="Plasma concentration (mg/L)") +
  theme(plot.title=element_text(size=20),
        axis.title=element_text(size=15, face="bold"),
        axis.text=element_text(size=10),
        legend.text=element_text(size=15),
        legend.position="")


### Child
#### load observed data; digitized using webplotdigitizer https://automeris.io/WebPlotDigitizer/
obsP <- read.csv("../data/Pediatric_IV.csv")  #pediatric IV infusion data; dose = 4 mg/kg and rate = 3 mg/kg

### plot
gp2 <- 
  ggplot() + 
  geom_point(data = obsP, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obsP, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = simP, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg IV", x="", y="") +
  theme(plot.title=element_text(size=20),
        axis.title=element_text(size=15, face="bold"),
        axis.text=element_text(size=10),
        legend.text=element_text(size=12))

# same plot without legend
gp3 <- 
  ggplot() + 
  geom_point(data = obsP, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obsP, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = simP, aes(x=time, y=CP, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(title="Pediatric 4 mg/kg IV", x="", y="") +
  theme(plot.title=element_text(size=20),
        axis.title=element_text(size=15, face="bold"),
        axis.text=element_text(size=10),
        legend.text=element_text(size=12),
        legend.position="")


### generate manuscript plot
legend <- get_legend(gp2)

grid1 <- plot_grid(gp1, gp3,
                   labels=c("a","b"),
                   label_size=20)

grid2 <- plot_grid(grid1, legend, rel_widths = c(1,.1))

ggdraw(add_sub(grid2, "time (h)", 
               vpadding=grid::unit(0,"lines"), y=6, x=0.5, vjust=4.5,
               fontface = "bold",
               size=15))

################################################################################################################
################################################################################################################


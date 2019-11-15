## This R script compiles and runs simulations for voriPBPK model; tasks here follow `hands-on_voriPBPK`

## load libraries
rm(list=ls())
library(dplyr)
library(ggplot2)
library(mrgsolve)


## adjust general theme for plotting
th <- theme(plot.title=element_text(size=30),
            axis.title=element_text(size=25),
            axis.text=element_text(size=20),
            legend.text=element_text(size=20))

#############################################################################################
########################################  Chunk 1  ##########################################
#############################################################################################

# compile extended model
modA_ext <- mread("../model/voriPBPK_ext")

#############################################################################################
#############################################################################################

#############################################################################################
#########################################  Chunk 2  #########################################
#############################################################################################

## run the IV simulations and compare to observed data

# Adult
obs <- read.csv("../data/Adult_IV.csv")  #load observed data

### Run this:
wt <- 73  #adult body weight
dose <- 4*wt  
rate <- 4*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA_ext %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, rate=rate, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th 
gp


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
                Qgu = 0.16*3.4*60, 
                Qhe = 0.04*3.4*60,
                Qki = 0.19*3.4*60,
                Qmu = 0.17*3.4*60,
                Qsp = 0.03*3.4*60,
                Qha = 0.065*3.4*60, 
                Qlu = 3.4*60,
                MPPGL = 26,
                VmaxH = 120.5,
                KmH = 11,
                MPPGI = 26/25,
                VmaxG = 120.5,
                KmG = 11,
                L = 170)

modP_ext <- param(modA_ext, pedPhys)

obs <- read.csv("../data/Pediatric_IV.csv")  #load observed data

### Run this:
wt <- 19  #adult body weight
dose <- 4*wt  
rate <- 3*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modP_ext %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, rate=rate, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp


#############################################################################################
#########################################  Chunk 2  #########################################
#############################################################################################
## Oral simulations
## Finetune the influential absorption parameters and run oral simulations
## fperm was optimized while MPPGI was scaled down from MPPGL assuming a 25-fold lower enzyme 
## expression in gut https://www.jstage.jst.go.jp/article/yakushi/123/5/123_5_369/_article

# adult
modA_ext <- param(modA_ext, fperm=0.47, MPPGI=30.3/25)

## load observed data
obs <- read.csv("../data/Adult_PO.csv")  #load observed data

wt <- 73
dose <- 200  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modA_ext %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=13, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp


## pediatric
modP_ext <- param(modA_ext, pedPhys)
obs <- read.csv("../data/Pediatric_PO.csv")  #load observed data

wt <- 19
dose <- 4*wt  
cmt <- "GUTLUMEN"  #intravenous infusion

# simulate
sim <- as.data.frame(modP_ext %>%
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

# plot
gp <- ggplot() + 
  geom_point(data = obs, aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, aes(x=time, y=Cvenous, col="sim"), lwd=1) + 
  scale_colour_manual(name='', 
                      values=c('sim'='black', 'observed'='black'), 
                      breaks=c("observed","sim"),
                      labels=c("observed","predicted")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
  labs(x="time (h)", y="Plasma concentration (mg/L)") +
  th
gp


#############################################################################################
#############################################################################################


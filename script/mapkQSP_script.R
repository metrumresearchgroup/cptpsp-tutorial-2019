## This script reproduces the results for MAPK QSP section in the manuscript. 
## Note: Before proceeding, make "script" your working directory.

# clear workspace
rm(list=ls())
gc()

# load libraries
.libPaths("lib")
library(tidyverse)
library(mrgsolve)
library(parallel)
#source("functions.R")
# mclapply <- lapply

# useful plotting functions
noline <- element_blank()
theme_plain <- function(...) {
  theme_bw() + theme(panel.grid.major=noline,panel.grid.minor=noline, 
                     plot.margin=margin(0.5,0.5,1,0.5,unit="cm"),...)
}
rotx <- function(angle=30) theme(axis.text.x = element_text(angle = angle, hjust = 1))
roty <- function(angle=30) theme(axis.text.y = element_text(angle = angle, hjust = 1))

################################################################################################################
#################################################### Chunk 1 ###################################################
################################################################################################################

# Read in the virtual population
vp <- readRDS("../data/s10vpop_pk.RDS") %>% 
  mutate(VPOP2 = seq(n()))


# Load the model and pick one parameter set from vpop

mod <- mread("../model/mapkQSP") %>% 
  update(end = 56)
mod <- param(mod, filter(vp, VPOP2==41))

################################################################################################################
################################################################################################################

################################################################################################################
#################################################### Chunk 2 ###################################################
################################################################################################################

# Predicting clinical outcomes for combination therapies

# Re-create figure 6B in the publication https://www.ncbi.nlm.nih.gov/pubmed/28649441

## Generate dosing regimens

### No treatment
data0 <- ev(amt=0, cmt=8)

### BRAF inhibitor - vemurafanib (VEMU) - Compartment 8
dataV <- ev(amt=960,  cmt=8, ii=0.5, addl=120)

### ERK inhibitor - GCD-994 (GDC) - Compartment 12
dataG <- ev(amt = 400, cmt = 12, ii = 1, addl = 20)
dataG <- seq(dataG, wait = 7, dataG) 

### MEK inhibitor - cobimetinib (COBI) - Compartment 10
dataCO <- mutate(dataG,amt=60,cmt=10)

### EGFR inihbitor - cetuximab (CETUX) - Compartment 7
dataCE <- ev(cmt=7,ii=7,addl=7,amt=450)


# We create two functions: one to combine dosing regimens and the other to simulate from a dosing regimen
comb <- function(...) {
  x <- lapply(list(...), as.data.frame)
  bind_rows(x) %>% arrange(time)
}

sim <- function(Data,Vp,Mod) {
  Mod %>%
    ev(as.ev(Data)) %>%
    mrgsim(idata=Vp, end=-1, add = 56) %>%
    filter(time==56) 
}


## Simulate all combination therapies

### Generate a data frame of runs to do

sims <- 
  tribble(
    ~label, ~object, 
    "No Treatment",        data0,
    "CETUX",               dataCE, 
    "VEMU",                dataV,
    "COBI",                dataCO, 
    "GDC",                 dataG,
    "CETUX+VEMU",          comb(dataCE, dataV), 
    "CETUX+COBI",          comb(dataCE, dataCO), 
    "CETUX+GDC",           comb(dataCE, dataG),
    "VEMU+COBI",           comb(dataV, dataG), 
    "VEMU+GDC",            comb(dataV, dataG),
    "COBI+GDC",            comb(dataCO, dataG),
    "CETUX+VEMU+COBI",     comb(dataCE, dataV,dataCO), 
    "CETUX+VEMU+GDC",      comb(dataCE, dataV,dataG), 
    "CETUX+COBI+GDC",      comb(dataCE, dataCO,dataG), 
    "VEMU+COBI+GDC",       comb(dataV, dataCO,dataG),
    "CETUX+VEMU+COBI+GDC", comb(dataCE, dataV, dataCO, dataG)
  ) %>% mutate(object = map(object,as.data.frame))



### Run the simulation
sims <- mutate(sims, out = parallel::mclapply(object,sim, Vp = vp, Mod = mod))


### Summarize and plot
sms <- select(sims, label, out) %>% unnest()
sms <- mutate(sms, 
              labelf = fct_inorder(label), 
              gdc = factor(grepl("GDC", label)))
head(sms)

p1 <- 
  ggplot(data=sms) + 
  geom_point(aes(x=labelf, y=TUMOR),position=position_jitter(width=0.15),col="grey") +
  geom_hline(yintercept=0.7,col="black", lty=1,lwd=0.7)  +
  scale_y_continuous(limits=c(0,2.5), name="Tumor size",breaks=c(0,0.5,1,1.5,2,2.5,3)) +
  scale_x_discrete(name="") + 
  geom_boxplot(aes(x=labelf,y=TUMOR,col=gdc),fill="darkslateblue",alpha=0.2) +
  scale_color_manual(values = c("darkslateblue", "firebrick"), guide = FALSE) + 
  theme_plain() + 
  rotx(30)
p1

################################################################################################################
################################################################################################################

################################################################################################################
#################################################### Chunk 3 ###################################################
################################################################################################################

# ORR analysis

## ORR in full population: GDC +/- COBI
sms %>%
  filter(label %in% c("GDC", "COBI+GDC")) %>%
  group_by(label) %>%
  summarise(orr = mean(TUMOR < 0.7)) %>% 
  knitr::kable(digits = 3)

## ORR in select patients: GDC +/- COBI
vp_select <- filter(vp, wOR > median(wOR))
re_run <- 
  sims %>%
  select(label,object) %>%
  filter(label %in% c("GDC", "COBI+GDC")) %>% 
  mutate(out = parallel::mclapply(object,sim,Vp = vp_select,Mod = mod)) %>%
  select(label, out) %>% 
  unnest()

re_run %>%
  group_by(label) %>%
  summarise(orr = mean(TUMOR < 0.7)) %>% 
  knitr::kable(digits = 3)

################################################################################################################
################################################################################################################

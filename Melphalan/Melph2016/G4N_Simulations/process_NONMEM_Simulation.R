#Script for simulating a population from a PK/PD model for melphalan
#------------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(R2HTML)
library(ggplot2)
library(doBy)
library(stringr)
library(Hmisc)
library(grid)
library(plyr)
library(reshape2)

#Use custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 16))  
theme_bw2 <- theme_update(plot.margin = unit(c(1.1,1.1,3,1.1), "lines"),
axis.title.x=element_text(size = 16, vjust = 0),
axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
strip.text.x=element_text(size = 14),
strip.text.y=element_text(size = 14, angle = 90))
	
#------------------------------------------------------------------------------------
#Set working directory
setwd("D:/Wojciechowski/Melphalan/G4N_Simulations/")

#Source functions_utility	
source("functions_utility.R")

# #------------------------------------------------------------------------------------
# #Run name
# runname1 <- "PKPD_Final_SIM_test_data_female"
  
# #Process the fit file
# processSIMdata(paste(runname1,".ctl",sep=""))
   
# #Read the simulated data
# SIM.data1 <- read.csv(paste(runname1,".nm7/",runname1,".fit.csv",sep=""), stringsAsFactors=F)
# SIM.data1 <- subset(SIM.data1, CMT == 4)  
	
# #Run name
# runname2 <- "PKPD_Final_SIM_test_data_male"
  
# #Process the fit file
# processSIMdata(paste(runname2,".ctl",sep=""))
   
# #Read the simulated data
# SIM.data2 <- read.csv(paste(runname2,".nm7/",runname2,".fit.csv",sep=""), stringsAsFactors=F)
# SIM.data2 <- subset(SIM.data2, CMT == 4) 

#------------------------------------------------------------------------------------
#Run name
runname1 <- "PKPD_Final_SIM_points_duration"
  
#Process the fit file
processSIMdata(paste(runname1,".ctl",sep=""))
   
#Read the simulated data
SIM.data1 <- read.csv(paste(runname1,".nm7/",runname1,".fit.csv",sep=""), stringsAsFactors=F)
SIM.data1 <- subset(SIM.data1, CMT == 4)  
 
#------------------------------------------------------------------------------------
#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)
	
#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

plotobj2 <- NULL
plotobj2 <- ggplot()
#plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = PRED), data = SIM.data1, colour = "red")
#plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = PRED), data = SIM.data2, colour = "blue")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = IPRED), data = SIM.data1, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.2)
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = IPRED), data = SIM.data2, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "blue", alpha = 0.2)
plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = PRED), data = SIM.data1[SIM.data1$ID == 1,], colour = "blue")
plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = IPRED), data = SIM.data1[SIM.data1$ID == 1,])
plotobj2 <- plotobj2 + geom_point(aes(x = TIME, y = DV), data = SIM.data1[SIM.data1$ID == 1,], colour = "red")
plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.5), linetype = "dashed")
plotobj2 <- plotobj2 + scale_y_log10("ANC\n")
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)", breaks = seq(from = 0, to = 800, by = 100))
plotobj2

head(SIM.data1)


G4N_EXIT <- max(SIM.data1$G4N_EXIT)
G4N_EXIT

G4NDUR <- max(SIM.data1$G4ND)
G4NDUR

G4N_ENTRY <- G4N_EXIT-G4NDUR
G4N_ENTRY



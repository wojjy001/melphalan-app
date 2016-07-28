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
setwd("/Volumes/Prosecutor/Melphalan/Simulation/")

#Source functions_utility	
source("functions_utility.R")

#------------------------------------------------------------------------------------
#Read the original data
ORG.data <- read.csv("PKPD_OSU11055_Neutropenia_lnplus4.csv", stringsAsFactors=F, na.strings=".")
ORG.data <- subset(ORG.data, X.Note != "#")
ORG.data <- ORG.data[,-1]	#Delete the "Notes" column
  
#Run name
runname <- "119pt_PKPD_Neutropenia_INPUTdelay_sim100"
  
#Process the fit file
#processSIMdata(paste(runname,".ctl",sep=""))
   
#Read the simulated data
SIM.data <- read.csv(paste(runname,".nm7/",runname,".fit.csv",sep=""), stringsAsFactors=F)
   
#Change working directory
setwd(paste(master.dir,"/",runname,".nm7",sep="")) 

#------------------------------------------------------------------------------------
#Subset PK and PD data
ORG.data$TIME <- as.numeric(ORG.data$TIME)
ORG.data$DV <- as.numeric(ORG.data$DV)
ORG.PK.data <- subset(ORG.data, DVID != 2)
ORG.PD.data <- subset(ORG.data, DVID == 2)
#Bin time - ORG.PK.data
ORG.PK.data$TIMEBIN <- cut2(ORG.PK.data$TIME, g=10, levels.mean=T)
ORG.PK.data$TIMEBIN <- as.numeric(paste(ORG.PK.data$TIMEBIN))
#Bin time - ORG.PD.data
ORG.PD.data$TIMEBIN <- cut2(ORG.PD.data$TIME, g=10, levels.mean=T)
ORG.PD.data$TIMEBIN <- as.numeric(paste(ORG.PD.data$TIMEBIN))
  
#------------------------------------------------------------------------------------
#Bin time - SIM.data
#Subset PK and PD data
SIM.PK.data <- subset(SIM.data, DVID != 2)
SIM.PD.data <- subset(SIM.data, DVID == 2)
#Bin time - ORG.PK.data
SIM.PK.data$TIMEBIN <- cut2(SIM.PK.data$TIME, g=10, levels.mean=T)
SIM.PK.data$TIMEBIN <- as.numeric(paste(SIM.PK.data$TIMEBIN))
#Bin time - ORG.PD.data
SIM.PD.data$TIMEBIN <- cut2(SIM.PD.data$TIME, g=10, levels.mean=T)
SIM.PD.data$TIMEBIN <- as.numeric(paste(SIM.PD.data$TIMEBIN))

#------------------------------------------------------------------------------------
#Plot - PK data
#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)
	
#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

sim.data.bystudy.median <- ddply(SIM.PK.data, .(SIM,TIMEBIN), function(df) median(df$DV)) 
sim.data.bystudy.median <- rename(sim.data.bystudy.median, c("V1"="medianS")) 

sim.data.bystudy.loCI <- ddply(SIM.PK.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
sim.data.bystudy.loCI <- rename(sim.data.bystudy.loCI, c("5%"="loCI90S")) 

sim.data.bystudy.hiCI <- ddply(SIM.PK.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
sim.data.bystudy.hiCI <- rename(sim.data.bystudy.hiCI, c("95%"="hiCI90S"))

sim.data.bystudy <- data.frame(sim.data.bystudy.median, "loCI90S"=sim.data.bystudy.loCI$loCI90S, "hiCI90S"=sim.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the SIM.data
plotobj1 <- NULL
plotobj1 <- ggplot(ORG.PK.data)
plotobj1 <- plotobj1 + geom_point(aes(x = TIME, y = DV), colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", size = 1)

#Lower 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

#Upper 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = median, geom = "line", colour = "red", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1) 

plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n")
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", breaks = c(0,2,4,6,8))
print(plotobj1)
	
#------------------------------------------------------------------------------------
#Plot - PD data
sim.data.bystudy.median <- ddply(SIM.PD.data, .(SIM,TIMEBIN), function(df) median(df$DV)) 
sim.data.bystudy.median <- rename(sim.data.bystudy.median, c("V1"="medianS")) 

sim.data.bystudy.loCI <- ddply(SIM.PD.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
sim.data.bystudy.loCI <- rename(sim.data.bystudy.loCI, c("5%"="loCI90S")) 

sim.data.bystudy.hiCI <- ddply(SIM.PD.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
sim.data.bystudy.hiCI <- rename(sim.data.bystudy.hiCI, c("95%"="hiCI90S"))

sim.data.bystudy <- data.frame(sim.data.bystudy.median, "loCI90S"=sim.data.bystudy.loCI$loCI90S, "hiCI90S"=sim.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the SIM.data
plotobj2 <- NULL
plotobj2 <- ggplot(ORG.PD.data)
plotobj2 <- plotobj2 + geom_point(aes(x = TIME, y = DV), colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", size = 1)

#Lower 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

#Upper 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = median, geom = "line", colour = "red", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1) 

plotobj2 <- plotobj2 + scale_y_continuous()
plotobj2 <- plotobj2 + ylab(expression(paste("Neutrophils (", 10^9, "/L)")))
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
print(plotobj2)
	  
  

#Script for simulating a population from a PK/PD model for melphalan
#------------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(deSolve) #Differential equation solver
library(ggplot2) #Plotting
library(plyr) #Split and rearrange data, ddply function
library(grid) #Plotting
library(MASS) #mvrnorm function
library(MBESS) #cor2cov function
library(compiler) #Compile repeatedly-called functions
library(Hmisc)  #cut2 function

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
#Read in the observed data
data <- read.csv("PKPD_OSU11055_Neutropenia_lnplus4.csv", stringsAsFactors=F, na.strings=".")
data <- subset(data, X.Note != "#")
data <- data[,-1]	#Delete the "Notes" column
data$TIME <- as.numeric(data$TIME)
data$ANC <- as.numeric(data$ANC)
data$DV <- as.numeric(data$DV)
	
#Repeat data frame for x number of simulations
x <- 200
SIM <- seq(from = 1, to = x, by = 1)
SIM <- rep(SIM, times=length(data$ID))
SIM <- sort(SIM)
par.data <- lapply(data, rep.int, times=x)
par.data$SIM <- SIM
par.data <- as.data.frame(par.data)

#------------------------------------------------------------------------------------
#Population Parameters
#THETAs
	#PK Parameters
	POPCL <- 29.3	#THETA1
	POPV1 <- 19.6	#THETA2
	POPQ <- 28.3	#THETA3
	POPV2 <- 20.8	#THETA4
	THETA_SHARE <- 1.31	#THETA5
	#PD Parameters
	POPBASE <- 5.2	#THETA6
	POPMTT <- 104	#THETA7
	POWER1 <- 0.209	#THETA8
	POPSLOPE <- 10.3	#THETA9
	KE <- 7	#THETA10
	POPIPT <- 12.9	#THETA11
	POPIP0 <- 19.9	#THETA12
	W <- 0.624	#THETA13
	#Covariate Effects
	COVCRCL <- 0.148	#THETA14
	COVSLC7A5 <- -0.116	#THETA15
		par.data$SLC7A5[is.na(par.data$SLC7A5 == T)] <- 0
	COVHCT <- 0.205	#THETA16
	
#OMEGAs (SDs)
	#PK Parameters
	PPVCL_V2 <- sqrt(0.0777)	#ETA shared between 2 parameters (CL and V2)
	PPVV1 <- sqrt(0.0681)
	PPVQ <- sqrt(0.287)
	#PD Parameters
	PPVBASE <- sqrt(0.111)
	PPVMTT <- sqrt(0.00987)
	PPVSLOPE <- sqrt(0.0778)
	PPVIPT <- 0
	PPVIP0 <- sqrt(2.39)
	
#SIGMAs (SDs)
	#PK Parameters
	ERR1 <- sqrt(0.0588)
	#PD Parameters
	ERR2 <- 1

#------------------------------------------------------------------------------------
#Calculate ETA values for each individual
ETA.data <- subset(par.data, DVID == 0)	#Faster than oneperID function
ETA.data <- ETA.data[-c(2:19)]	#Remove dosing/concentration information
	
#Use nrorm to calculate ETA values
ETA.data$ETA1 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVCL_V2)
ETA.data$ETA2 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVV1)
ETA.data$ETA3 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVQ)
ETA.data$ETA4 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVBASE)
ETA.data$ETA5 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVMTT)
ETA.data$ETA6 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVSLOPE)
ETA.data$ETA7 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVIPT)
ETA.data$ETA8 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVIP0)
	
#Merge ETA.data with par.data
new.par.data <- merge(par.data,ETA.data,by=c("ID","SIM"),all=T)
	
#Calculate individual parameter values
	#Population
		#PK Parameters
		new.par.data$TVCL <- POPCL*((new.par.data$CrCL/91.94)^COVCRCL)*((new.par.data$FFM/59.9)^0.75)*((new.par.data$HCT/32.5)^COVHCT)
		new.par.data$TVV1 <- POPV1*(new.par.data$FFM/59.9)
		new.par.data$TVQ <- POPQ*((new.par.data$FFM/59.9)^0.75)
		new.par.data$TVV2 <- POPV2*(1+new.par.data$SLC7A5*COVSLC7A5)*(new.par.data$FFM/59.9)
		#PD Parameters
		new.par.data$TVBASE <- POPBASE
		new.par.data$TVMTT <- POPMTT
		new.par.data$TVSLOPE <- POPSLOPE
		new.par.data$TVIPT <- POPIPT
		new.par.data$TVIP0 <-  POPIP0
	#Individual
		#PK Parameters
		new.par.data$CL <- new.par.data$TVCL*exp(new.par.data$ETA1)
		new.par.data$V1 <- new.par.data$TVV1*exp(new.par.data$ETA2)
		new.par.data$Q <- new.par.data$TVQ*exp(new.par.data$ETA3)
		new.par.data$V2 <- new.par.data$TVV2*exp(new.par.data$ETA1*THETA_SHARE)
		#PD Parameters
		new.par.data$BASE <- new.par.data$TVBASE*exp(new.par.data$ETA4)
		new.par.data$MTT <- new.par.data$TVMTT*exp(new.par.data$ETA5)
		new.par.data$SLOPE <- new.par.data$TVSLOPE*exp(new.par.data$ETA6)
		new.par.data$IPT <- new.par.data$TVIPT*exp(new.par.data$ETA7)
		new.par.data$IP0 <- new.par.data$TVIP0*exp(new.par.data$ETA8)
		
#------------------------------------------------------------------------------------
#Make a reduced data frame to be used in differential equation solver
PK.data <- data.frame(ID = new.par.data$ID,
											SIM = new.par.data$SIM,
											DVID = new.par.data$DVID,
											AMT = new.par.data$AMT,
											RATE = new.par.data$RATE,
											CL = new.par.data$CL,
											V1 = new.par.data$V1,
											Q = new.par.data$Q,
											V2 = new.par.data$V2)
PK.data <- subset(PK.data, DVID == 0)
PK.data <- PK.data[,-3]	#Remove the DVID column

PK.data <- PK.data[with(PK.data, order(PK.data$SIM,PK.data$ID)), ]

#Find the unique times from the PK data
TIME <- unique(new.par.data$TIME[new.par.data$DVID != 2])
TIME <- sort(TIME)

#Add in the time-points that signify when an individuals infusion has finished
#Required for the differential equation solver
INFDUR <- PK.data$AMT/PK.data$RATE
TIMEs <- sort(unique(c(TIME,INFDUR)))

#Test with just a few individuals
#PK.data <- PK.data[PK.data$SIM < 11, ]

#------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, THETA) {
	#Infusion function
	RATE <- THETA[6]
	INFEND <- THETA[5]/RATE	#Duration = AMT/RATE
	END <- max(T)+1
	#Vector marking infusion's time events
	TIMEinf <- c(0,INFEND,END)
	RATEinf <- c(RATE,0,0)  #Vector marking infusion's rates
				
	#Define an interpolation function that returns rate when given time - "const"
	Cstep.doseinf <- approxfun(TIMEinf, RATEinf, method = "const")
	RATEC <- Cstep.doseinf(T)
	
	#Differential equations
	K10 <- THETA[1]/THETA[2]  #Elimination rate constant (CL/V1)
	K12 <- THETA[3]/THETA[2]
	K21 <- THETA[3]/THETA[4]
	dAdT <- vector(length = 2)
	dAdT[1] = RATEC -K12*A[1] +K21*A[2] -K10*A[1]	#Central compartment
	dAdT[2] = K12*A[1] -K21*A[2]	#Peripheral compartment		
	list(dAdT)
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#Function for simulating concentrations for the ith patient
simulate.conc <- function(PK.data) {		
	#List of parameters from input for the differential equation solver			
	THETAlist <- c(CL = PK.data$CL,
									V1 = PK.data$V1,
									Q = PK.data$Q,
									V2 = PK.data$V2,
									AMT = PK.data$AMT,
									RATE = PK.data$RATE)	
	#Set initial compartment conditions
	A_0 <- c(A1 = 0, A2 = 0)	
	#Run differential equation solver for simulated variability data	
	var.data <- lsoda(y = A_0, TIMEs, func = DES.cmpf, parms = THETAlist)
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)

#Apply simulate.conc.cmpf to each individual in par.data
sim.data <- ddply(PK.data, .(ID,SIM,V1), simulate.conc.cmpf)

#------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$IPRE <- sim.data$A1/sim.data$V1
	
#------------------------------------------------------------------------------------
#Calculate individual score with RUV
EPS <- rnorm(n = length(sim.data$IPRE), mean = 0, sd = ERR1)
sim.data$DV <- sim.data$IPRE*exp(EPS)

#------------------------------------------------------------------------------------
#Plot
#Bin time - sim.data (same bins as NONMEM simulation for comparison!)
names(sim.data)[4] <- "TIME"
sim.data$TIMEBIN <- NA
sim.data$TIMEBIN[sim.data$time == 0] <- 0
sim.data$TIMEBIN[sim.data$time >= 0.42 & sim.data$time < 0.63] <- 0.53395
sim.data$TIMEBIN[sim.data$time >= 0.63 & sim.data$time < 0.78] <- 0.69291
sim.data$TIMEBIN[sim.data$time >= 0.78 & sim.data$time < 1.05] <- 0.93180
sim.data$TIMEBIN[sim.data$time >= 1.05 & sim.data$time < 1.28] <- 1.15598
sim.data$TIMEBIN[sim.data$time >= 1.28 & sim.data$time < 1.48] <- 1.35821
sim.data$TIMEBIN[sim.data$time >= 1.48 & sim.data$time < 1.65] <- 1.55054
sim.data$TIMEBIN[sim.data$time >= 1.65 & sim.data$time < 3.58] <- 2.82141
sim.data$TIMEBIN[sim.data$time >= 3.58 & sim.data$time < 6.52] <- 4.53580
sim.data$TIMEBIN[sim.data$time >= 6.52 & sim.data$time <= 7.75] <- 6.64495

#Pull out just the time-points from the original dataset for each individual
ID.data <- new.par.data[new.par.data$DVID != 2, ]
ID.data <- data.frame(ID = ID.data$ID,
													SIM = ID.data$SIM,
													TIME = ID.data$TIME)
new.sim.data <- merge(ID.data,sim.data,by=c("ID","SIM","TIME"),all=F)

#Bin time - data (original data)
data$TIMEBIN <- NA
data$TIMEBIN[data$TIME == 0] <- 0
data$TIMEBIN[data$TIME >= 0.42 & data$TIME < 0.63] <- 0.53395
data$TIMEBIN[data$TIME >= 0.63 & data$TIME < 0.78] <- 0.69291
data$TIMEBIN[data$TIME >= 0.78 & data$TIME < 1.05] <- 0.93180
data$TIMEBIN[data$TIME >= 1.05 & data$TIME < 1.28] <- 1.15598
data$TIMEBIN[data$TIME >= 1.28 & data$TIME < 1.48] <- 1.35821
data$TIMEBIN[data$TIME >= 1.48 & data$TIME < 1.65] <- 1.55054
data$TIMEBIN[data$TIME >= 1.65 & data$TIME < 3.58] <- 2.82141
data$TIMEBIN[data$TIME >= 3.58 & data$TIME < 6.52] <- 4.53580
data$TIMEBIN[data$TIME >= 6.52 & data$TIME <= 7.75] <- 6.64495

#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)

#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

sim.data.bystudy.median <- ddply(new.sim.data, .(SIM,TIMEBIN), function(df) median(df$DV)) 
sim.data.bystudy.median <- rename(sim.data.bystudy.median, c("V1"="medianS")) 

sim.data.bystudy.loCI <- ddply(new.sim.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
sim.data.bystudy.loCI <- rename(sim.data.bystudy.loCI, c("5%"="loCI90S")) 

sim.data.bystudy.hiCI <- ddply(new.sim.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
sim.data.bystudy.hiCI <- rename(sim.data.bystudy.hiCI, c("95%"="hiCI90S"))

sim.data.bystudy <- data.frame(sim.data.bystudy.median, "loCI90S"=sim.data.bystudy.loCI$loCI90S, "hiCI90S"=sim.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the sim.data
plotobj1 <- NULL
plotobj1 <- ggplot(new.sim.data)
plotobj1 <- plotobj1 + geom_point(aes(x = TIME, y = DV), data = data[data$DVID == 1,], colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", size = 1)

#Lower 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

#Upper 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = sim.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = data[data$DVID == 1,], fun.y = median, geom = "line", colour = "red", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = data[data$DVID == 1,], fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = data[data$DVID == 1,], fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1) 
	
plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n")
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", lim = c(0,8), breaks = c(0,2,4,6,8))
print(plotobj1)
	


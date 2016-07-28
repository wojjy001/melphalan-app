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
	K <- 4/MTT
	POWER1 <- 0.209	#THETA8
	POPSLOPE <- 10.3	#THETA9
	KE <- log(2)/7	#THETA10
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
data <- data.frame(ID = new.par.data$ID,
											SIM = new.par.data$SIM,
											DVID = new.par.data$DVID,
											AMT = new.par.data$AMT,
											RATE = new.par.data$RATE,
											CL = new.par.data$CL,
											V1 = new.par.data$V1,
											Q = new.par.data$Q,
											V2 = new.par.data$V2,
											BASE = new.par.data$BASE,
											SLOPE = new.par.data$SLOPE,
											MTT = new.par.data$MTT,
											IPT = new.par.data$IPT,
											IP0 = new.par.data$IP0)
data <- subset(data, DVID == 0)
data <- data[,-3]	#Remove the DVID column

data <- data[with(data, order(data$SIM,data$ID)), ]

#Find the unique times from the PK data
TIME <- unique(new.par.data$TIME)
TIME <- sort(TIME)

#Add in the time-points that signify when an individuals infusion has finished
#Required for the differential equation solver
INFDUR <- data$AMT/data$RATE
TIMEs <- sort(unique(TIME,INFDUR))

#Test with just a few individuals
#data <- data[data$ID == 1, ]

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

	#Differential equations for PK
	K10 <- THETA[1]/THETA[2]  #Elimination rate constant (CL/V1)
	K12 <- THETA[3]/THETA[2]
	K21 <- THETA[3]/THETA[4]
	dAdT <- vector(length = 8)
	dAdT[1] = RATEC -K12*A[1] +K21*A[2] -K10*A[1]	#Central compartment
	dAdT[2] = K12*A[1] -K21*A[2]	#Peripheral compartment
	
	#PD parameters
	CP <- A[1]/THETA[2]	#Plasma concentration in the central compartment
	BASE <- THETA[7]
	SLOPE <- THETA[8]
	MTT <- THETA[9]
	IPT <- THETA[10]
	IP0 <- THETA[11]
	DRUG <- SLOPE*CP
	K <- 4/MTT
	POWER1 <- 0.209
	KE <- log(2)/7
	
	#KIN function
	TIMEkin <- c(0,72,END)
	RATEkin <- c(0,1/IPT,1/IPT)
	Kstep.dosekin <- approxfun(TIMEkin, RATEkin, method = "const")
	KIN <- Kstep.dosekin(T)
	
	#FN
	if (A[4] > 0.001) { FN <- (BASE/A[4])^POWER1 }
	if (A[4] < 0.001) { FN <- (BASE/0.001)^POWER1 }
	
	#Differential equations for PD
	dAdT[3] = K*A[3]*(1-DRUG)*FN -K*A[3]
	dAdT[4] = K*A[7] -KE*A[4] +KIN*A[8]
	dAdT[5] = K*A[3] -K*A[5]
	dAdT[6] = K*A[5] -K*A[6]
	dAdT[7] = K*A[6] -K*A[7]
	dAdT[8] = -KIN*A[8]
	list(dAdT)
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#Function for simulating concentrations for the ith patient
simulate.conc <- function(data) {		
	#List of parameters from input for the differential equation solver			
	THETAlist <- c(CL = data$CL,
									V1 = data$V1,
									Q = data$Q,
									V2 = data$V2,
									AMT = data$AMT,
									RATE = data$RATE,
									BASE = data$BASE,
									SLOPE = data$SLOPE,
									MTT = data$MTT,
									IPT = data$IPT,
									IP0 = data$IP0)	
	#Set initial compartment conditions
	A_0 <- c(A1 = 0, A2 = 0, C3 = (log(2)/7*data$BASE)/(4/data$MTT), C4 = data$BASE, C5 = (log(2)/7*data$BASE)/(4/data$MTT), C6 = (log(2)/7*data$BASE)/(4/data$MTT), C7 = (log(2)/7*data$BASE)/(4/data$MTT), C8 = data$IP0)
	#Run differential equation solver for simulated variability data	
	var.data <- lsoda(y = A_0, TIMEs, func = DES.cmpf, parms = THETAlist)
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)

#Apply simulate.conc.cmpf to each individual in par.data
sim.data <- ddply(data, .(ID,SIM,V1), simulate.conc.cmpf)

#------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$IPRE <- sim.data$A1/sim.data$V1
sim.data$IPRED <- log(sim.data$C4)+4
	
#------------------------------------------------------------------------------------
#Calculate individual predictions with RUV
EPS1 <- rnorm(n = length(sim.data$IPRE), mean = 0, sd = ERR1)
sim.data$DV <- sim.data$IPRE*exp(EPS1)

EPS2 <- rnorm(n = length(sim.data$IPRED), mean = 0, sd = ERR2)
sim.data$DVD <- sim.data$IPRED + W*EPS2

#------------------------------------------------------------------------------------
#Plot - PK
#Bin time - PK.data (same bins as NONMEM simulation for comparison!)
PK.data <- sim.data
names(PK.data)[4] <- "TIME"
PK.data$TIMEBIN <- NA
PK.data$TIMEBIN[PK.data$TIME == 0] <- 0
PK.data$TIMEBIN[PK.data$TIME >= 0.42 & PK.data$TIME < 0.63] <- 0.53395
PK.data$TIMEBIN[PK.data$TIME >= 0.63 & PK.data$TIME < 0.78] <- 0.69291
PK.data$TIMEBIN[PK.data$TIME >= 0.78 & PK.data$TIME < 1.05] <- 0.93180
PK.data$TIMEBIN[PK.data$TIME >= 1.05 & PK.data$TIME < 1.28] <- 1.15598
PK.data$TIMEBIN[PK.data$TIME >= 1.28 & PK.data$TIME < 1.48] <- 1.35821
PK.data$TIMEBIN[PK.data$TIME >= 1.48 & PK.data$TIME < 1.65] <- 1.55054
PK.data$TIMEBIN[PK.data$TIME >= 1.65 & PK.data$TIME < 3.58] <- 2.82141
PK.data$TIMEBIN[PK.data$TIME >= 3.58 & PK.data$TIME < 6.52] <- 4.53580
PK.data$TIMEBIN[PK.data$TIME >= 6.52 & PK.data$TIME <= 7.75] <- 6.64495

#Pull out just the time-points from the original dataset for each individual
ID.data <- new.par.data[new.par.data$DVID != 2, ]
ID.data <- data.frame(ID = ID.data$ID,
													SIM = ID.data$SIM,
													TIME = ID.data$TIME)
new.PK.data <- merge(ID.data,PK.data,by=c("ID","SIM","TIME"),all=F)

#Bin time - data (original data)
new.par.data$TIMEBIN <- NA
new.par.data$TIMEBIN[new.par.data$TIME == 0] <- 0
new.par.data$TIMEBIN[new.par.data$TIME >= 0.42 & new.par.data$TIME < 0.63] <- 0.53395
new.par.data$TIMEBIN[new.par.data$TIME >= 0.63 & new.par.data$TIME < 0.78] <- 0.69291
new.par.data$TIMEBIN[new.par.data$TIME >= 0.78 & new.par.data$TIME < 1.05] <- 0.93180
new.par.data$TIMEBIN[new.par.data$TIME >= 1.05 & new.par.data$TIME < 1.28] <- 1.15598
new.par.data$TIMEBIN[new.par.data$TIME >= 1.28 & new.par.data$TIME < 1.48] <- 1.35821
new.par.data$TIMEBIN[new.par.data$TIME >= 1.48 & new.par.data$TIME < 1.65] <- 1.55054
new.par.data$TIMEBIN[new.par.data$TIME >= 1.65 & new.par.data$TIME < 3.58] <- 2.82141
new.par.data$TIMEBIN[new.par.data$TIME >= 3.58 & new.par.data$TIME < 6.52] <- 4.53580
new.par.data$TIMEBIN[new.par.data$TIME >= 6.52 & new.par.data$TIME <= 7.75] <- 6.64495

#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)

#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

PK.data.bystudy.median <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) median(df$DV)) 
PK.data.bystudy.median <- rename(PK.data.bystudy.median, c("V1"="medianS")) 

PK.data.bystudy.loCI <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
PK.data.bystudy.loCI <- rename(PK.data.bystudy.loCI, c("5%"="loCI90S")) 

PK.data.bystudy.hiCI <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
PK.data.bystudy.hiCI <- rename(PK.data.bystudy.hiCI, c("95%"="hiCI90S"))

PK.data.bystudy <- data.frame(PK.data.bystudy.median, "loCI90S"=PK.data.bystudy.loCI$loCI90S, "hiCI90S"=PK.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the PK.data (PK data)
plotobj1 <- NULL
plotobj1 <- ggplot(new.PK.data)
plotobj1 <- plotobj1 + geom_point(aes(x = TIME, y = DV), data = new.par.data[new.par.data$DVID == 1,], colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = PK.data.bystudy, fun.y = median, geom = "line", colour = "black", size = 1)

#Lower 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = PK.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

#Upper 90% CI simulated with confidence band
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = PK.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 1,], fun.y = median, geom = "line", colour = "red", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 1,], fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 1,], fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1) 
	
plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n")
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", lim = c(0,8), breaks = c(0,2,4,6,8))
print(plotobj1)
	
#------------------------------------------------------------------------------------
#Plot - PD
#Bin time - PD.data (same bins as NONMEM simulation for comparison!)
PD.data <- sim.data
names(PD.data)[4] <- "TIME"
PD.data$TIMEBIN <- NA
PD.data$TIMEBIN[PD.data$TIME < 72] <- 24.096
PD.data$TIMEBIN[PD.data$TIME >= 72 & PD.data$TIME < 120] <- 83.949
PD.data$TIMEBIN[PD.data$TIME == 120] <- 120
PD.data$TIMEBIN[PD.data$TIME >= 144 & PD.data$TIME < 192] <- 156
PD.data$TIMEBIN[PD.data$TIME >= 192 & PD.data$TIME < 204] <- 204
PD.data$TIMEBIN[PD.data$TIME >= 240] <- 240
PD.data$TIMEBIN[PD.data$TIME >= 264 & PD.data$TIME < 312] <- 276
PD.data$TIMEBIN[PD.data$TIME >= 312 & PD.data$TIME < 360] <- 323.63
PD.data$TIMEBIN[PD.data$TIME >= 360 & PD.data$TIME < 408] <- 370.447
PD.data$TIMEBIN[PD.data$TIME >= 408 & PD.data$TIME <= 552] <- 456.736

#Pull out just the time-points from the original dataset for each individual
ID.data <- new.par.data[new.par.data$DVID == 2, ]
ID.data <- data.frame(ID = ID.data$ID,
													SIM = ID.data$SIM,
													TIME = ID.data$TIME)
new.PD.data <- merge(ID.data,PD.data,by=c("ID","SIM","TIME"),all=F)

#Bin time - data (original data)
new.par.data$TIMEBIN <- NA
new.par.data$TIMEBIN[new.par.data$TIME < 72] <- 24.096
new.par.data$TIMEBIN[new.par.data$TIME >= 72 & new.par.data$TIME < 120] <- 83.949
new.par.data$TIMEBIN[new.par.data$TIME == 120] <- 120
new.par.data$TIMEBIN[new.par.data$TIME >= 144 & new.par.data$TIME < 192] <- 156
new.par.data$TIMEBIN[new.par.data$TIME >= 192 & new.par.data$TIME < 204] <- 204
new.par.data$TIMEBIN[new.par.data$TIME >= 240] <- 240
new.par.data$TIMEBIN[new.par.data$TIME >= 264 & new.par.data$TIME < 312] <- 276
new.par.data$TIMEBIN[new.par.data$TIME >= 312 & new.par.data$TIME < 360] <- 323.63
new.par.data$TIMEBIN[new.par.data$TIME >= 360 & new.par.data$TIME < 408] <- 370.447
new.par.data$TIMEBIN[new.par.data$TIME >= 408 & new.par.data$TIME <= 552] <- 456.736

#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)

#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

PD.data.bystudy.median <- ddply(new.PD.data, .(SIM,TIMEBIN), function(df) median(df$DVD)) 
PD.data.bystudy.median <- rename(PD.data.bystudy.median, c("V1"="medianS")) 

PD.data.bystudy.loCI <- ddply(new.PD.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DVD))
PD.data.bystudy.loCI <- rename(PD.data.bystudy.loCI, c("5%"="loCI90S")) 

PD.data.bystudy.hiCI <- ddply(new.PD.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DVD))
PD.data.bystudy.hiCI <- rename(PD.data.bystudy.hiCI, c("95%"="hiCI90S"))

PD.data.bystudy <- data.frame(PD.data.bystudy.median, "loCI90S"=PD.data.bystudy.loCI$loCI90S, "hiCI90S"=PD.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the PK.data (PK data)
plotobj2 <- NULL
plotobj2 <- ggplot(new.PD.data)
plotobj2 <- plotobj2 + geom_point(aes(x = TIME, y = DV), data = new.par.data[new.par.data$DVID == 2,], colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = medianS), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = medianS), data = PD.data.bystudy, fun.y = median, geom = "line", colour = "black", size = 1)

#Lower 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = PD.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

#Upper 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
#plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = PD.data.bystudy, fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 2,], fun.y = median, geom = "line", colour = "red", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 2,], fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIMEBIN, y = DV), data = new.par.data[new.par.data$DVID == 2,], fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1) 
	
plotobj2 <- plotobj2 + scale_y_continuous("Cell Count\n")
#plotobj2 <- plotobj2 + ylab(expression(paste("Neutrophils (", 10^9, "/L)")))
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
print(plotobj2)


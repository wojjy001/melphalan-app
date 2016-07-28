#Simulating Population
#R Script for simulating warfarin concentrations and INR for a population
#Model reference: Hamberg et al (2010) "A Pharmacometric Model Describing the Relationship Between Warfarin Dose and INR Response with Respect to Variations in CYP2C9, VKORC1 and Age"
#------------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(deSolve)	#Differential equation solver
library(ggplot2)	#Plotting
library(plyr)	#Split and rearrange data, ddply function
library(grid)	#Plotting
library(MASS)		#mvrnorm function
library(MBESS)		#cor2cov function
library(compiler)	#Compile repeatedly-called functions

#------------------------------------------------------------------------------------
#Set a directory for where plots can be saved (best where this R script is saved)
setwd("/Volumes/Prosecutor/PhD/BayesEstimation/Time-Weighting/NovProject1/")

#Define a custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 16))  
theme_bw2 <- theme_update(plot.margin = unit(c(1.1,1.1,3,1.1), "lines"),
axis.title.x=element_text(size = 16, vjust = 0),
axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
strip.text.x=element_text(size = 14),
strip.text.y=element_text(size = 14, angle = 90))

#Function for calculating the median, and 5th and 95th percentiles for plotting simulation results
sumfuncx <- function(x) {
	stat1 <- median(x)
	stat2 <- quantile(x, probs=0.05, names=F) 
	stat3 <- quantile(x, probs=0.95, names=F)
	stat4 <- length(x)
	result <- c("median"=stat1, "low"=stat2, "hi"=stat3, "n"=stat4)
	result
}

#------------------------------------------------------------------------------------
#Assign patient population characteristics
#Parameters in this section are regarded as "input" and will become widgets when converting to a Shiny application
n <- 2	#Number of individuals to be simulated

#Set up a TIME sequence
SINT <- 24	#Simulation interval (hours)
TIME <- seq(from=0, to=30*24, by=SINT)	#30x24-hour days every SINT hour(s)

CYP2C9 <- 11	#CYP2C9 Genotype - 11, 12, 13, 22, 23, or 33
AGE <- 70+TIME/24/365	#Age (years)
VKORC1 <- "GG"	#VKORC1 Genotype - GG, GA, AA

DOSE <- 7.6
AMT <- DOSE/2	#Amount (mg) of S-warfarin

#------------------------------------------------------------------------------------
#Define the values for the model's population parameters
#THETAs
	#PK Parameters
	if (CYP2C9 == 11) POPCL <- 0.174 + 0.174
	if (CYP2C9 == 12) POPCL <- 0.174 + 0.0879	
	if (CYP2C9 == 13) POPCL <- 0.174 + 0.0422	
	if (CYP2C9 == 22) POPCL <- 0.0879 + 0.0879 
	if (CYP2C9 == 23) POPCL <- 0.0879 + 0.0422 
	if (CYP2C9 == 33) POPCL <- 0.0422 + 0.0422
	POPV1 <- 14.3	
	POPKA <- 2
	#PD Parameters
	POPEMAX <- 1
	POPHILL <- 1.15	
	if (VKORC1 == "GG") POPEC50 <- 2.05 + 2.05
	if (VKORC1 == "GA") POPEC50 <- 2.05 + 0.96
	if (VKORC1 == "AA") POPEC50 <- 0.96 + 0.96
	POPMTT1 <- 28.6
	POPMTT2 <- 118.3
	#Covariate Effects
	COVAGE <- -0.00571	#Effect of age (change/year) on CL
	
#OMEGAs (SDs)
	#PK Parameters
	PPVCL <- 29.8/100
	PPVV1 <- 23.2/100
		#Covariance
		RCLV1 <- 0.84
	#PD Parameters
	PPVEC50 <- 33.2/100
	
#SIGMAs (SDs)
	#PK Parameters
	ERR1 <- 0.099/100
	#PD Parameters
	ERR2 <- 20/100

#------------------------------------------------------------------------------------
#Use mvrnorm and rnorm(generates random numbers from the specified normal distribution density) to simulate ETA values for the number of individuals in the population
	#Calculate ETA values for each subject
	CORR <- matrix(c(1,RCLV1,RCLV1,1),2,2)
	#Specify the between subject variability for CL, V1
	SDVAL <- c(PPVCL,PPVV1)
	#Use this function to turn CORR and SDVAL into a covariance matrix
	OMEGA <- cor2cov(CORR,SDVAL)
	#Now use multivariate rnorm to turn the covariance matrix into ETA values
	ETAmat <- mvrnorm(n=n, mu=c(0,0), OMEGA)   
	if (n > 1) ETA1 <- ETAmat[,1]
	if (n > 1) ETA2 <- ETAmat[,2]
	if (n > 1) ETA3 <- rnorm(n=n, mean=0, sd=PPVEC50)
		#When number of individuals (n) == 1, only use the population values
		if (n == 1) ETA1 <- 0
		if (n == 1) ETA2 <- 0
		if (n == 1) ETA3 <- 0
	
#Calculate individual parameter values from population typical values
	#Population Typical Values
		#PK Parameters
		TVCL <- POPCL*(1+COVAGE*(AGE-71))
	#Individual Parameter Value
		#PK Parameters
		V1 <- POPV1*exp(ETA2)	#THETA1
		#PD Parameters
		EC50 <- POPEC50*exp(ETA3)	#THETA7
		
#------------------------------------------------------------------------------------
#Make a data frame of individual parameter information to be entered into the differential equation solver - "input.data"
#Each individual needs one row in the data frame
#Firstly, make a list of ID numbers so that one individual can differentiated from the others
ID <- seq(from=1, to=n, by=1)
#Make a data frame of individual-specific parameters linked with ID's
input.data <- data.frame(ID,V1,EC50,ETA1)

#Define dose events and times when CL changes
if(SINT < 24) FACT <- 1	#Factor for affecting replicate times
if(SINT == 24) FACT <- 0
event.data <- data.frame(var=c(rep(1,times=length(TIME)/24*SINT+FACT),	#Dose
																rep(9,times=length(TIME)) 	#Changing CL
																),
													time=c(seq(from=min(TIME),to=max(TIME),by=24),
																TIME
																),
													value=c(rep(AMT,times=length(TIME)/24*SINT+FACT),
																TVCL
																),
													method=c(rep("add",times=length(TIME)/24*SINT+FACT),
																rep("rep",times=length(TIME))
																))

#------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, THETA) {
	#9 differential equations
	#2 x PK compartments, 6 x PD compartments, 1 x time-varying parameter(s)	
	dAdt <- vector(length=9)

	##############
	##_PK BLOCK_##
	##############
	#Define which values in the THETA vector (defined as THETAlist below) are the PK parameters
	CL <- A[9]*exp(THETA[3])
	V1 <- THETA[1]
	KA <- POPKA
	KDE <- CL/V1

	#Differential equations for PK
	dAdt[1] = -KA*A[1]	#Absorption compartment
	dAdt[2] = KA*A[1] -KDE*A[2]	#Central compartment

	##############
	##_PD BLOCK_##
	##############
	#Define which values in the THETA vector (defined as THETAlist below) are the PD parameters
	EMAX <- POPEMAX
	HILL <- POPHILL
	EC50 <- THETA[2]
	MTT1 <- 3/POPMTT1
	MTT2 <- 3/POPMTT2
	
	#Calculate the dose driving rate (DR) and drug effect (EFF)
	DR <- KDE*A[2]
	EDK50 <- CL*EC50	#Dose rate that leads to 50% inhibition
	EFF <- (EMAX*(DR^HILL))/((EDK50^HILL)+(DR^HILL))

	#Differential equations for PD
	dAdt[3] = MTT1*(1-EFF) -MTT1*A[3]	#TRANSIT11
	dAdt[4] = MTT1*A[3] -MTT1*A[4]	#TRANSIT12
	dAdt[5] = MTT1*A[4] -MTT1*A[5]	#TRANSIT13
	dAdt[6] = MTT2*(1-EFF) -MTT2*A[6]	#TRANSIT21
	dAdt[7] = MTT2*A[6] -MTT2*A[7]	#TRANSIT22
	dAdt[8] = MTT2*A[7] -MTT2*A[8]	#TRANSIT23
	
	################
	##_PARAMETERS_##
	################
	#Set rate to zero so parameters don't change unless there is an event
	dAdt[9] = 0		#CL
		
	list(dAdt)
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#------------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.INR <- function(input.data) {		
	#List of parameters from input for the differential equation solver			
	THETAlist <- c(V1 = input.data$V1,
									EC50 = input.data$EC50,
									ETA1 = input.data$ETA1)
									
	#Set initial compartment conditions
	A1_0 <- 0	#Absorption compartment (PK)
	A2_0 <- 0	#Central compartment (PK)
	C3_0 <- 1	#TRANSIT11 (PD)
	C4_0 <- 1 #TRANSIT12 (PD)
	C5_0 <- 1	#TRANSIT13 (PD)
	C6_0 <- 1	#TRANSIT21 (PD)
	C7_0 <- 1 #TRANSIT22 (PD)
	C8_0 <- 1	#TRANSIT23 (PD)
	P9_0 <- event.data$value[event.data$var == 9 & event.data$time == 0]	#CL (Parameter)
	
	#Vector it!
	A_0 <- c(A1=A1_0, A2=A2_0, C3=C3_0, C4=C4_0, C5=C5_0, C6=C6_0, C7=C7_0, C8=C8_0, P9=P9_0)
	
	#Run differential equation solver for simulated variability data	
	var.data <- ode(y=A_0, TIME, func=DES.cmpf, parms=THETAlist, events=list(data=event.data), method="lsodar")
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.INR.cmpf <- cmpfun(simulate.INR)

#Apply simulate.conc.cmpf to each individual in par.data
#Maintain their individual value for V1 for later calculations
sim.data <- ddply(input.data, .(ID,V1), simulate.INR.cmpf)

#------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$CONC <- sim.data$A2/sim.data$V1
sim.data$INR <- 1+20*(1-(sim.data$C5+sim.data$C8)/2)
	
#------------------------------------------------------------------------------------
#At each time-point ("time") in sim.data, calculate the median, and 5th and 95th percentiles for predictions
statsINR <- ddply(sim.data, .(time), function(sim.data) sumfuncx(sim.data$INR))
names(statsINR)[c(2,3,4)] <- c("median","low","hi")

#Plot INR simulation results		
plotobj1 <- NULL
plotobj1 <- ggplot(statsINR)
plotobj1 <- plotobj1 + geom_line(aes(x=time/24, y=median), colour="red", size=1)
plotobj1 <- plotobj1 + geom_ribbon(aes(x=time/24, ymin=low, ymax=hi), fill="red", alpha=0.3)
plotobj1 <- plotobj1 + scale_y_continuous("INR\n")
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (days)")
print(plotobj1)

#------------------------------------------------------------------------------------
#At each time-point ("time") in sim.data, calculate the median, and 5th and 95th percentiles for predictions
statsCONC <- ddply(sim.data, .(time), function(sim.data) sumfuncx(sim.data$CONC))
names(statsCONC)[c(2,3,4)] <- c("median","low","hi")

#Plot CONC simulation results		
plotobj2 <- NULL
plotobj2 <- ggplot(statsCONC)
plotobj2 <- plotobj2 + geom_line(aes(x=time/24, y=median), colour="blue", size=1)
plotobj2 <- plotobj2 + geom_ribbon(aes(x=time/24, ymin=low, ymax=hi), fill="blue", alpha=0.3)
plotobj2 <- plotobj2 + scale_y_continuous("S-Warfarin Concentration (mg/L)\n")
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (days)")
print(plotobj2)

##------------------------------------------------------------------------------------
##Write individual simulated data to .csv file
#ind.data <- sim.data[sim.data$ID == 2, ]
##Add RUV to CONC and INR predictions
#ind.data$CONC_DV <- ind.data$CONC*exp(rnorm(length(ind.data$CONC),mean=0,sd=ERR1))
#ind.data$INR_DV <- ind.data$INR*exp(rnorm(length(ind.data$INR),mean=0,sd=ERR2))
#ind.data <- data.frame(ID = ind.data$ID,
#				TIME = ind.data$time,
#				CYP2C9,
#				VKORC1,
#				AGE,
#				V1 = V1[2],
#				EC50 = EC50[2],
#				EXP_ETA1 = exp(ETA1[2]),
#				CONC = ind.data$CONC_DV,
#				INR = ind.data$INR_DV)
#write.csv(ind.data, file="single_individual_data.csv", na=".", quote=F, row.names=F)

#Script that simulates melphalan concentrations and neutrophils for a single individual "n" times with specified characteristics 
#Based on the NONMEM control stream: 119pt_PKPD_Neutropenia_INPUTdelay_sim100.ctl
#Due to the number of differential equations, simulation time is quite slow
#------------------------------------------------------------------------------------
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(deSolve)	#Differential equation solver
library(ggplot2)	#Plotting
library(plyr)	#Split and rearrange data, ddply function
library(grid)	#Plotting
library(compiler)	#Compile repeatedly-called functions

#------------------------------------------------------------------------------------
#Set a directory for where plots can be saved (best where this R script is saved)
setwd("/Volumes/Prosecutor/Melphalan/Simulation/")

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
n <- 100	#Number of individuals to be simulated

FFM <- 59.9	#Fat free mass
CRCL <- 91.94	#Creatinine clearance
SLC7A5 <- 0	#0 or 1
HCT <- 32.5	#HCT

AMT <- 400	#Amount to be infused
INFD <- 0.5	#Infusion duration
RATE <- AMT/INFD	#Rate of infusion

#------------------------------------------------------------------------------------
#Define the values for the model's population parameters
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
	KE <- log(2)/7	#THETA10
	POPIPT <- 12.9	#THETA11
	POPIP0 <- 19.9	#THETA12
	W <- 0.624	#THETA13
	#Covariate Effects
	COVCRCL <- 0.148	#THETA14
	COVSLC7A5 <- -0.116	#THETA15
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
#Use nrorm (generates random numbers from the specified normal distribution density) to simulate ETA values for the number of individuals in the population
ETA1 <- rnorm(n = n, mean = 0, sd = PPVCL_V2)
ETA2 <- rnorm(n = n, mean = 0, sd = PPVV1)
ETA3 <- rnorm(n = n, mean = 0, sd = PPVQ)
ETA4 <- rnorm(n = n, mean = 0, sd = PPVBASE)
ETA5 <- rnorm(n = n, mean = 0, sd = PPVMTT)
ETA6 <- rnorm(n = n, mean = 0, sd = PPVSLOPE)
ETA7 <- rnorm(n = n, mean = 0, sd = PPVIPT)
ETA8 <- rnorm(n = n, mean = 0, sd = PPVIP0)
	
#Calculate individual parameter values from population typical values
	#Population Typical Values
		#PK Parameters
		TVCL <- POPCL*((CRCL/91.94)^COVCRCL)*((FFM/59.9)^0.75)*((HCT/32.5)^COVHCT)
		TVV1 <- POPV1*(FFM/59.9)
		TVQ <- POPQ*((FFM/59.9)^0.75)
		TVV2 <- POPV2*(1+SLC7A5*COVSLC7A5)*(FFM/59.9)
		#PD Parameters
		TVBASE <- POPBASE
		TVMTT <- POPMTT
		TVSLOPE <- POPSLOPE
		TVIPT <- POPIPT
		TVIP0 <-  POPIP0
	#Individual Parameter Value
		#PK Parameters
		CL <- TVCL*exp(ETA1)
		V1 <- TVV1*exp(ETA2)
		Q <- TVQ*exp(ETA3)
		V2 <- TVV2*exp(ETA1*THETA_SHARE)
		#PD Parameters
		BASE <- TVBASE*exp(ETA4)
		MTT <- TVMTT*exp(ETA5)
		K <- 4/MTT
		SLOPE <- TVSLOPE*exp(ETA6)
		IPT <- TVIPT*exp(ETA7)
		IP0 <- TVIP0*exp(ETA8)
		
#------------------------------------------------------------------------------------
#Make a data frame of individual parameter information to be entered into the differential equation solver - "input.data"
#Each individual needs one row in the data frame
#Firstly, make a list of ID numbers so that one individual can differentiated from the others
ID <- seq(from = 1, to = n, by = 1)
input.data <- data.frame(ID,CL,V1,Q,V2,BASE,MTT,K,SLOPE,IPT,IP0)

#Set up a TIME sequence specifying PK sampling times
PK.TIME <- seq(from = 0, to = 8, by = 0.5)

#Set up a TIME sequence specifying PD sampling times
PD.TIME <- seq(from = 0, to = 408, by = 24)

#Combine the unique PK and PD times (only one TIME sequence can be input to the differential equation solver)
#Also, add in the time-point that signifies when the infusion has finished (INFD, if not already done so)
TIME <- sort(unique(c(PK.TIME,PD.TIME,INFD)))

#------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, THETA) {
	#8 differential equations
	#2 x PK compartments, 6 x PD compartments	
	dAdT <- vector(length = 8)

	##############
	##_PK BLOCK_##
	##############
	#Define which values in the THETA vector (defined as THETAlist below) are the PK parameters
	CL <- THETA[1]
	V1 <- THETA[2]
	Q <- THETA[3]
	V2 <- THETA[4]
	#Specify the rate constants
	K10 <- CL/V1
	K12 <- Q/V1
	K21 <- Q/V2

	#Define the infusion (when it starts, when it finishes, and the rate) - this uses the "approxfun" function to make a "forcing function" for infusion rate in the differential equations
	#The function needs to continue long after the infusion ends - specify an "end" time for the function
	END <- max(T)+1	#Maximum value in the TIME sequence + 1
	#Specify a vector that marks the infusion's time events
	TIMEinf <- c(0,INFD,END)	#Something happens at TIME = 0 and TIME = INFD
	#Specify a vector marking the infusion's rates	
	RATEinf <- c(RATE,0,0)	#At TIME = 0, RATE = RATE and TIME = INFD, RATE = 0 (i.e., infusion has ended, and continues to be zero until the end of the function) 
	#Define an interpolation function that returns the rate when given a time - "const"
	Cstep.doseinf <- approxfun(TIMEinf, RATEinf, method = "const")
	RATEC <- Cstep.doseinf(T)

	#Differential equations for PK
	dAdT[1] = RATEC -K12*A[1] +K21*A[2] -K10*A[1]	#Central compartment, PKCENTR
	dAdT[2] = K12*A[1] -K21*A[2]	#Peripheral compartment, PKPERI

	##############
	##_PD BLOCK_##
	##############
	#Define which values in the THETA vector (defined as THETAlist below) are the PD parameters
	BASE <- THETA[5]
	MTT <- THETA[6]
	K <- THETA[7]
	POWER1 <- 0.209
	SLOPE <- THETA[8]
	KE <- log(2)/7
	IPT <- THETA[9]
	IP0 <- THETA[10]
	
	#Calculate the plasma concentration in the central compartment and drug effect
	CP <- A[1]/V1
	DRUG <- SLOPE*CP
	
	#FN
	if (A[4] > 0.001) {FN <- (BASE/A[4])^POWER1}
	if (A[4] < 0.001) {FN <- (BASE/0.001)^POWER1}

	#KIN function - same concept as for coding the infusion
	TIMEkin <- c(0,72,END)
	RATEkin <- c(0,1/IPT,1/IPT)
	Kstep.dosekin <- approxfun(TIMEkin, RATEkin, method = "const")
	KIN <- Kstep.dosekin(T)

	#Differential equations for PD
	dAdT[3] = K*A[3]*(1-DRUG)*FN -K*A[3]	#STEM
	dAdT[4] = K*A[7] -KE*A[4] +KIN*A[8]	#ANC
	dAdT[5] = K*A[3] -K*A[5]	#TRANSIT1
	dAdT[6] = K*A[5] -K*A[6]	#TRANSIT2
	dAdT[7] = K*A[6] -K*A[7]	#TRANSIT3
	dAdT[8] = -KIN*A[8]	#INPUT2
		
	list(dAdT)
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#------------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.conc <- function(input.data) {		
	#List of parameters from input for the differential equation solver			
	THETAlist <- c(CL = input.data$CL,
									V1 = input.data$V1,
									Q = input.data$Q,
									V2 = input.data$V2,
									BASE = input.data$BASE,
									MTT = input.data$MTT,
									K = input.data$K,
									SLOPE = input.data$SLOPE,
									IPT = input.data$IPT,
									IP0 = input.data$IP0)	
	#Set initial compartment conditions
	A1_0 <- 0	#Central compartment, PKCENTR (PK)
	A2_0 <- 0	#Peripheral compartment, PKPERI (PK)
	C3_0 <- (log(2)/7*input.data$BASE)/(input.data$K)	#STEM (PD)
	C4_0 <- input.data$BASE	#ANC (PD)
	C5_0 <- (log(2)/7*input.data$BASE)/(input.data$K)	#TRANSIT1 (PD)
	C6_0 <- (log(2)/7*input.data$BASE)/(input.data$K)	#TRANSIT2 (PD)
	C7_0 <- (log(2)/7*input.data$BASE)/(input.data$K)	#TRANSIT3 (PD)
	C8_0 <- input.data$IP0	#INPUT2 (PD)
	#Vector it!
	A_0 <- c(A1 = A1_0, A2 = A2_0, C3 = C3_0, C4 = C4_0, C5 = C5_0, C6 = C6_0, C7 = C7_0, C8 = C8_0)
	#Run differential equation solver for simulated variability data	
	var.data <- lsoda(y = A_0, TIME, func = DES.cmpf, parms = THETAlist)
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)

#Apply simulate.conc.cmpf to each individual in par.data
#Maintain their individual value for V1 for later calculations
sim.data <- ddply(input.data, .(ID,V1), simulate.conc.cmpf)

#------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$PKIPRE <- sim.data$A1/sim.data$V1
sim.data$PDIPRE <- log(sim.data$C4)+4
	
#------------------------------------------------------------------------------------
#Calculate individual predictions with RUV
#Simulate residuals using a random number generator (rnorm) again
EPS1 <- rnorm(n = length(sim.data$PKIPRE), mean = 0, sd = ERR1)
sim.data$PKDV <- sim.data$PKIPRE*exp(EPS1)

EPS2 <- rnorm(n = length(sim.data$PDIPRE), mean = 0, sd = ERR2)
sim.data$PDDV <- sim.data$PDIPRE + W*EPS2

#------------------------------------------------------------------------------------
#Make a dataset only containing time-points relevant for the PK part (i.e., those in PK.TIME)
PK.data <- sim.data[sim.data$time %in% PK.TIME, ]
#At each time-point ("time") in PK.data, calculate the median, and 5th and 95th percentiles for predictions
statsPK <- ddply(PK.data, .(time), function(PK.data) sumfuncx(PK.data$PKDV))
names(statsPK)[c(2,3,4)] <- c("median","low","hi")

#Plot PK simulation results
plotobj1 <- NULL
plotobj1 <- ggplot(statsPK)
plotobj1 <- plotobj1 + geom_line(aes(x = time, y = median), colour = "red", size = 1)
plotobj1 <- plotobj1 + geom_ribbon(aes(x = time, ymin = low, ymax = hi), fill = "red", alpha = 0.3)
plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n")
plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", lim = c(0,8), breaks = c(0,2,4,6,8))
print(plotobj1)
ggsave("PKSIM.png")
	
#------------------------------------------------------------------------------------
#Make a dataset only containing time-points relevant for the PD part (i.e., those in PD.TIME)
PD.data <- sim.data[sim.data$time %in% PD.TIME, ]
#At each time-point ("time") in PD.data, calculate the median, and 5th and 95th percentiles for predictions
statsPD <- ddply(PD.data, .(time), function(PD.data) sumfuncx(PD.data$PDDV))
names(statsPD)[c(2,3,4)] <- c("median","low","hi")

#Plot PD simulation results		
plotobj2 <- NULL
plotobj2 <- ggplot(statsPD)
plotobj2 <- plotobj2 + geom_line(aes(x = time, y = median), colour = "blue", size = 1)
plotobj2 <- plotobj2 + geom_ribbon(aes(x = time, ymin = low, ymax = hi), fill = "blue", alpha = 0.3)
plotobj2 <- plotobj2 + ylab(expression(paste("Neutrophils (", 10^9, "/L)")))
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
print(plotobj2)
ggsave("PDSIM.png")

#------------------------------------------------------------------------------------
#Set up a data frame that can be input into NONMEM
#1. Repeat the IDs for the length of PK times and length of PD times
	PK.ID <- rep(ID, times=length(PK.TIME))
	PK.ID <- sort(PK.ID)
	PD.ID <- rep(ID, times=length(PD.TIME))
	PD.ID <- sort(PD.ID)
	
#2. Make separate PK and PD data frames of ID and TIME and DVID
	PK.ID.data <- data.frame(ID = PK.ID, TIME = PK.TIME, DVID = 1)
	PK.ID.data <- PK.ID.data[PK.ID.data$TIME != 0, ]
	PD.ID.data <- data.frame(ID = PD.ID, TIME = PD.TIME, DVID = 2)
		
#3. Make a data frame of IDs and TIME = 0 and AMT and RATE and DVID = 0	
	DOSE.data <- data.frame(ID, TIME = 0, DVID = 0)
	
#4. rbind PK.ID.data and PD.ID.data, and merge with DOSE.data
	new.data <- rbind(PK.ID.data,PD.ID.data,DOSE.data)
	new.data <- new.data[with(new.data, order(new.data$ID,new.data$TIME,new.data$DVID)),]

#5. Add AMT and RATE columns and corresponding values for when DVID == 0
	new.data$AMT <- "."
	new.data$AMT[new.data$DVID == 0] <- AMT
	new.data$RATE <- "."
	new.data$RATE[new.data$DVID == 0] <- RATE
	
#6. Set up DV, MDV and CMT columns
	new.data$DV <- "."
	new.data$MDV <- 0
	new.data$MDV[new.data$DVID == 0] <- 1
	new.data$CMT <- 1
	new.data$CMT[new.data$DVID == 2] <- 4

#7. Add in required covariate values
	new.data$FFM <- 59.9	#Fat free mass
	new.data$CRCL <- 91.94	#Creatinine clearance
	new.data$SLC7A5 <- 0	#0 or 1
	new.data$HCT <- 32.5	#HCT
	
#8. Write to a .csv file
	write.csv(new.data, file="dummy_data.csv", na=".", quote=F, row.names=F)





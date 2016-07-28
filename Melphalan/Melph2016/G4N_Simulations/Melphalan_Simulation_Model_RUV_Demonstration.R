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
#setwd("/Volumes/Prosecutor/Melphalan/Melph2016/")

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

#Function for calculating the mean, and 95% confidence intervals
sumfuncy <- function(x) {
	stat1 <- mean(x)
	stat2 <- sd(x)
	stat3 <- length(x)	
	stat4 <- stat1 - 1.96*(stat2/sqrt(stat3))
	stat5 <- stat1 + 1.96*(stat2/sqrt(stat3))
	result <- c("mean"=stat1, "sd"=stat2, "n"=stat3, "lower95CI"=stat4, "upper95CI"=stat5)
	result
}

#------------------------------------------------------------------------------------
#Assign patient population characteristics
#Parameters in this section are regarded as "input" and will become widgets when converting to a Shiny application
n <- 1000	#Number of individuals to be simulated

AMT <- 400	#Amount to be infused
INFD <- 0.5	#Infusion duration
RATE <- AMT/INFD	#Rate of infusion

#Set up a TIME sequence specifying PK sampling times
PK.TIME <- seq(from = 0, to = 8, by = 0.5)

#Set up a TIME sequence specifying PD sampling times
PD.TIME <- seq(from = 0, to = 408, by = 24)

#Combine the unique PK and PD times (only one TIME sequence can be input to the differential equation solver)
#Also, add in the time-point that signifies when the infusion has finished (INFD, if not already done so)
TIME <- sort(unique(c(PK.TIME,PD.TIME,INFD)))

FFM <- 59.9	#Fat free mass
CRCL <- 91.94	#Creatinine clearance
SLC7A5 <- 0	#0 or 1
HCT <- 32.5	#Haematocrit
ANCBASE <- 3.5	#Baseline absolute neutrophil count
SEX <- 1	#Gender, male = 0, female = 1
BUN <- 14	#Blood urea nitrogen
LNP53FOLD <- 2.62	#If values were -99, they were made to be 1
WBC <- 4.9	#White blood cells
GCSF <- 0	#GCSF administered, yes or no

#------------------------------------------------------------------------------------
#Define the values for the model's population parameters
#THETAs
	#PK Parameters
	POPCL <- 28.9	#THETA1
	POPV1 <- 19.5	#THETA2
	POPQ <- 28.5	#THETA3
	POPV2 <- 20.9	#THETA4
	THETA_SHARE <- 1.31		#THETA5
	COVCRCL <- 0.193 #THETA6
	COVSLC7A5 <- -0.130 #THETA7
	COVHCT <- 0.216 #THETA8
	#PD Parameters
	W <- 0.693	#THETA9
	POPBASE <- 5.65	#THETA10
	POPSLOPE <- 7.4	#THETA11
	POPMTT <- 95.8	#THETA12
	COVGCSF <- 0.182	#THETA13
	COVNEUP <- 0.0364	#THETA14
	KE <- log(2)/7	#THETA15
	POPIP01 <- 110 #THETA16
	POPIP02 <- 0.035	#THETA17
	POPIPT <- 14.7	#THETA18
	GCSF_MTT <- 0.153	#THETA19
	ANC_BASE <- 0.225	#THETA20
	GCSF_BASE <- 0.405 #THETA21
	BUN_MTT <- 0.0629	#THETA22
	SEX_SLOPE <- 0.208	#THETA23
	HCT_SLOPE <- 0.653	#THETA24
	LNP53_IP0 <- -0.777	#THETA25
	HCT_BASE <- 0.584	#THETA26
	FAC9 <- 0.411 #THETA27
	
#OMEGAs (SDs)
	#PK Parameters
	PPVCL_V2 <- sqrt(0.0792)	#ETA1 shared between 2 parameters (CL and V2)
	PPVV1 <- sqrt(0.0740)	#ETA2
	PPVQ <- sqrt(0.297)	#ETA3
	#PD Parameters
	PPVBASE <- sqrt(0.101)	#ETA4
	PPVSLOPE <- sqrt(0.0608)	#ETA5
	PPVMTT <- sqrt(0.00674)	#ETA6
	PPVIP0 <- sqrt(0.159)	#ETA7
	PPVIPT <- 0	#ETA8
	
#SIGMAs (SDs)
	#PK Parameters
	ERR1 <- sqrt(0.0584)	#EPS1
	#PD Parameters
	ERR2 <- 1	#EPS2

#------------------------------------------------------------------------------------
#Use nrorm (generates random numbers from the specified normal distribution density) to simulate ETA values for the number of individuals in the population
# ETA1 <- rnorm(n = n, mean = 0, sd = PPVCL_V2)
# ETA2 <- rnorm(n = n, mean = 0, sd = PPVV1)
# ETA3 <- rnorm(n = n, mean = 0, sd = PPVQ)
# ETA4 <- rnorm(n = n, mean = 0, sd = PPVBASE)
# ETA5 <- rnorm(n = n, mean = 0, sd = PPVSLOPE)
# ETA6 <- rnorm(n = n, mean = 0, sd = PPVMTT)
# ETA7 <- rnorm(n = n, mean = 0, sd = PPVIP0)
# ETA8 <- rnorm(n = n, mean = 0, sd = PPVIPT)

ETA1 <- 0
ETA2 <- 0
ETA3 <- 0
ETA4 <- 0
ETA5 <- 0
ETA6 <- 0
ETA7 <- 0
ETA8 <- 0
		
#------------------------------------------------------------------------------------
#Make a data frame of individual parameter information to be entered into the differential equation solver - "input.data"
#Each individual needs one row in the data frame
#Firstly, make a list of ID numbers so that one individual can differentiated from the others
ID <- seq(from = 1, to = n, by = 1)
input.data <- data.frame(ID,ETA1,ETA2,ETA3,ETA4,ETA5,ETA6,ETA7,ETA8,
												CRCL,FFM,HCT,SLC7A5,BUN,ANCBASE,SEX,WBC,LNP53FOLD,GCSF)
												
# #Set up an "events" data frame for when GCSF (Neupogen) is administered
# event.data <- data.frame(var = rep(9, times = length(TIME)),	#Goes into compartment 9
													# time = TIME,
													# value = GCSF,
													# method = rep("rep", times = length(TIME)))

#------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, PAR) {
	#9 differential equations
	#2 x PK compartments, 6 x PD compartments, 1 x Statistical compartments
	dAdT <- vector(length = 9)
	
	################
	##_PARAMETERS_##
	################
	#Define which values in the vector (defined as PARAMETER.list below) are the parameters
	#Time-independent covariates
	ETA1 <- PAR[1]	#Shared ETA between CL and V2
	ETA2 <- PAR[2]	#PPV for V1
	ETA3 <- PAR[3]	#PPV for Q
	ETA4 <- PAR[4]	#PPV for BASE
	ETA5 <- PAR[5]	#PPV for SLOPE
	ETA6 <- PAR[6]	#PPV for MTT
	ETA7 <- PAR[7]	#PPV for IP0
	ETA8 <- PAR[8]	#PPV for IPT
	CRCL <- PAR[9]	#Creatinine clearance
	FFM <- PAR[10]	#Fat free mass
	HCT <- PAR[11]	#Haematocrit
	SLC7A5 <- PAR[12]	#Some genetic thing
	BUN <- PAR[13]	#Blood urea nitrogen
	ANCBASE <- PAR[14]	#Baseline absolute neutrophil count
	SEX <- PAR[15]	#Gender
	WBC <- PAR[16]	#White blood cells
	LNP53FOLD <- PAR[17]	#Some other genetic things
	GCSF <- PAR[18]	#Neupogen

	##############
	##_PK BLOCK_##
	##############
	#Individual PK parameters
	CL <- POPCL*((CRCL/91.94)^COVCRCL)*((FFM/59.9)^0.75)*((HCT/32.5)^COVHCT)*exp(ETA1)
	V1 <- POPV1*(FFM/59.9)*exp(ETA2)
	Q <- POPQ*((FFM/59.9)^0.75)*exp(ETA3)
	V2 <- POPV2*(1+SLC7A5*COVSLC7A5)*(FFM/59.9)*exp(ETA1*THETA_SHARE)
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
	#PD Parameters
	BASE <- POPBASE*exp(ETA4)*((ANCBASE/3.5)^ANC_BASE)*((HCT/32.5)^HCT_BASE)
	SLOPE <- POPSLOPE*exp(ETA5)*(1+GCSF*GCSF_BASE)*(1+SEX*SEX_SLOPE)*((HCT/32.5)^HCT_SLOPE)
	MTT <- POPMTT*exp(ETA6)*(1+GCSF*GCSF_MTT)*((BUN/14)^BUN_MTT)	
	K <- 4/MTT
	KE <- log(2)/7
	
	ON <- 1
	if (GCSF == 1) {ON <- 0}
	IP0 <- (POPIP01*ON+(1-ON)*POPIP02)*exp(ETA7)*((LNP53FOLD/2.62)^LNP53_IP0)*((WBC/4.9)^FAC9)	
	IPT <- POPIPT*exp(ETA8)
	
	#Calculate the plasma concentration in the central compartment and drug effect
	CP <- A[1]/V1
	DRUG <- SLOPE*CP
	
	#QN functions
	QN <- 0
	if (GCSF == 0 & T > 72) {QN <- 1}
	if (GCSF == 1 & T > 216) {QN <- 1}
	POWER1 <- COVGCSF + QN*COVNEUP
	
	#FN
	FN <- (BASE/A[4])^POWER1
	if (A[4] < 0.001) {FN <- (BASE/0.001)^POWER1}

	#KIN functions
	KIN <- 0
	if (GCSF == 0 & T > 72) {KIN <- 1/IPT}
	if (GCSF == 1 & T > 216) {KIN <- 1/IPT}

	#Differential equations for PD
	dAdT[3] = K*A[3]*(1-DRUG)*FN -K*A[3]	#STEM
	dAdT[4] = K*A[7] -KE*A[4] +KIN*A[8]	#ANC
	dAdT[5] = K*A[3] -K*A[5]	#TRANSIT1
	dAdT[6] = K*A[5] -K*A[6]	#TRANSIT2
	dAdT[7] = K*A[6] -K*A[7]	#TRANSIT3
	dAdT[8] = -KIN*A[8]	#INPUT
		
	################
	##_STAT BLOCK_##
	################	
	dAdT[9] = 0	#STAT
	if (A[4] < 0.5) {dAdT[9] = 1} 
		
	list(dAdT)	#Get it to give the "amount" in each compartment
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#------------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.conc <- function(input.data) {		
	#Print a message stating which ID and SIM it is up to
	print(paste("ID:",input.data$ID,sep=" "))

	#List of parameters from input for the differential equation solver			
	PARAMETER.list <- c(ETA1 = input.data$ETA1,
											ETA2 = input.data$ETA2,
											ETA3 = input.data$ETA3,
											ETA4 = input.data$ETA4,
											ETA5 = input.data$ETA5,
											ETA6 = input.data$ETA6,
											ETA7 = input.data$ETA7,
											ETA8 = input.data$ETA8,
											CRCL = input.data$CRCL,
											FFM = input.data$FFM,
											HCT = input.data$HCT,
											SLC7A5 = input.data$SLC7A5,
											BUN = input.data$BUN,
											ANCBASE = input.data$ANCBASE,
											SEX = input.data$SEX,
											WBC = input.data$WBC,
											LNP53FOLD = input.data$LNP53FOLD,
											GCSF = input.data$GCSF)	
	
	#Some individual parameters required for initial conditions
	BASEI <- POPBASE*exp(input.data$ETA4)*((input.data$ANCBASE/3.5)^ANC_BASE)*((input.data$HCT/32.5)^HCT_BASE)
	MTTI <- POPMTT*exp(input.data$ETA6)*(1+input.data$GCSF*GCSF_MTT)*((input.data$BUN/14)^BUN_MTT)
			
	ON <- 1
	if (input.data$GCSF == 1) {ON <- 0}
	IP0I <- (POPIP01*ON+(1-ON)*POPIP02)*exp(input.data$ETA7)*((input.data$LNP53FOLD/2.62)^LNP53_IP0)*((input.data$WBC/4.9)^FAC9)
			
	#Set initial compartment conditions
	A1_0 <- 0	#Central compartment, PKCENTR (PK)
	A2_0 <- 0	#Peripheral compartment, PKPERI (PK)
	C3_0 <- (KE*BASEI)/(4/MTTI)	#STEM (PD)
	C4_0 <- BASEI	#ANC (PD)
	C5_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT1 (PD)
	C6_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT2 (PD)
	C7_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT3 (PD)
	C8_0 <- IP0I	#INPUT (PD)
	S9_0 <- 0	#STAT compartment for calculating duration of severe neutropenia
	
	#Vector it!
	A_0 <- c(A1 = A1_0,A2 = A2_0,C3 = C3_0,C4 = C4_0,C5 = C5_0,C6 = C6_0,C7 = C7_0,C8 = C8_0,S9 = S9_0)
	
	#Run differential equation solver for simulated variability data	
	var.data <- lsoda(y=A_0, TIME, func=DES.cmpf, parms=PARAMETER.list)
	var.data <- as.data.frame(var.data)	
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)

#Apply simulate.conc.cmpf to each individual in par.data
#Maintain their individual value for FFM and ETA2 for later calculations of V1
sim.data <- ddply(input.data, .(ID,FFM,ETA2), simulate.conc.cmpf)

#------------------------------------------------------------------------------------
#Calculate individual predictions
sim.data$V1 <- POPV1*(sim.data$FFM/59.9)*exp(sim.data$ETA2)
sim.data$PKIPRE <- sim.data$A1/sim.data$V1
sim.data$PDIPRE <- sim.data$C4
	
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
#ggsave("R_PKSIM.png")
	
#------------------------------------------------------------------------------------
#Make a dataset only containing time-points relevant for the PD part (i.e., those in PD.TIME)
PD.data <- sim.data[sim.data$time %in% PD.TIME, ]
#At each time-point ("time") in PD.data, calculate the median, and 5th and 95th percentiles for predictions
PD.data$PDDV[PD.data$PDDV < 0] <- NA
statsPD <- ddply(PD.data, .(time), function(PD.data) sumfuncx(na.omit(PD.data$PDDV)))
names(statsPD)[c(2,3,4)] <- c("median","low","hi")

#At each time-point ("time") in PD.data, calculate the mean and 95%CI
statsPD2 <- ddply(PD.data, .(time), function(PD.data) sumfuncy(na.omit(PD.data$PDDV)))
names(statsPD2)[c(2,3,4,5,6)] <- c("mean","sd","n","lower95CI","upper95CI")

#Plot PD simulation results		
plotobj2 <- NULL
plotobj2 <- ggplot()
plotobj2 <- plotobj2 + geom_line(aes(x = time, y = median), data = statsPD, colour = "blue")
plotobj2 <- plotobj2 + geom_ribbon(aes(x = time, ymin = low, ymax = hi), data = statsPD, fill = "blue", alpha = 0.3)
plotobj2 <- plotobj2 + geom_line(aes(x = time, y = mean), data = statsPD2, colour = "red")
plotobj2 <- plotobj2 + geom_ribbon(aes(x = time, ymin = lower95CI, ymax = upper95CI), data = statsPD2, fill = "red", alpha = 0.3)
plotobj2 <- plotobj2 + geom_line(aes(x = time, y = PDIPRE), data = sim.data[sim.data$PDIPRE >= 0,])
plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.5),linetype = "dashed")
plotobj2 <- plotobj2 + scale_y_log10("ANC (K/uL)\n", breaks = c(0.1,1,10))
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
print(plotobj2)
#ggsave("R_PDSIM.png")

#Pull out the duration in grade 4 neutropenia
G4ND <- max(PD.data$S9)
G4ND

#Find the time when the patient enters grade 4 neutropenia
G4ND_ENTRY <- PD.data$time[PD.data$S9 == min(PD.data$S9[PD.data$S9 > 0])] - min(PD.data$S9[PD.data$S9 > 0])
G4ND_ENTRY

#Find the time when the patient exits grade 4 neutropenia
G4ND_EXIT <- G4ND_ENTRY + G4ND
G4ND_EXIT



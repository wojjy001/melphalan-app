#Define global functions for MelphPrototypeApp1
#Code for functions and variables which are not reactive (i.e., not dependent on user input)
#Elements that only need to run once on application start-up
#Supplies these functions and objects to server.R
#-----------------------------------------------------------------------------------
#Load package libraries
library(shiny)	#Infrastructure for the application
library(deSolve)	#Differential equation solver package
library(plyr)	#Split, apply and rearrange data package
library(ggplot2)	#Plotting package
library(grid)	#Extra plotting features
library(compiler)	#Compile repeatedly called functions for speed benefits

#Custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 16))
theme_bw2 <- theme_update(plot.margin = unit(c(1,1.3,1,1),"lines"),
							axis.title.x=element_text(size = 16,vjust = 0),
							axis.title.y=element_text(size = 16,vjust = 0,angle = 90),
							strip.text.x=element_text(size = 14),
							strip.text.y=element_text(size = 14,angle = 90),
							legend.direction = "vertical",
              legend.position = "bottom",
              legend.box = "horizontal")

n <- 50	#Number of individuals to simulate for 95% prediction intervals

#Function for calculating the median, and 2.5th and 97.5th percentiles for plotting simulation results
CI95lo <- function(x) quantile(x,probs = 0.025)
CI95hi <- function(x) quantile(x,probs = 0.975)

#Function for taking the last row of each individual's profile
lastperID <- function(x) {tail(x,1)}

#------------------------------------------------------------------------------------
#Define time sequence
#Set up a TIME sequence specifying PK sampling times
PK.TIME <- seq(from = 0,to = 10,by = 0.5)
#Set up a TIME sequence specifying PD sampling times
PD.TIME <- seq(from = 0,to = 720,by = 24)

#------------------------------------------------------------------------------------
#Define the values for the model's population parameters
#THETAs
	#PK Parameters
	POPCL <- 28.9	#THETA1
	POPV1 <- 19.5	#THETA2
	POPQ <- 28.5	#THETA3
	POPV2 <- 20.9	#THETA4
	THETA_SHARE <- 1.31	#THETA5
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
#Simulate string of ETA values for prediction intervals
#Not built into server.R yet
ETA1 <- rnorm(n,mean = 0,sd = PPVCL_V2)
ETA2 <- rnorm(n,mean = 0,sd = PPVV1)
ETA3 <- rnorm(n,mean = 0,sd = PPVQ)
ETA4 <- rnorm(n,mean = 0,sd = PPVBASE)
ETA5 <- rnorm(n,mean = 0,sd = PPVSLOPE)
ETA6 <- rnorm(n,mean = 0,sd = PPVMTT)
ETA7 <- rnorm(n,mean = 0,sd = PPVIP0)
ETA8 <- rnorm(n,mean = 0,sd = PPVIPT)

#-----------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T,A,PAR) {
	#T = TIME, A = AMOUNT, PAR = PARAMETER
	#9 differential equations
	#2 x PK compartments, 6 x PD compartments, 1 x STAT compartment
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
	BSA <- PAR[19]	#Body surface area (m^2)
	DOSE <- PAR[20]	#Dose of melphalan to be infused (mg/m^2)
	AMT <- BSA*DOSE	#Amount of melphalan to be infused (mg)
	INFD <- 0.5	#Infusion duration (fixed to 30 minutes)
	RATE <- AMT/INFD	#Infusion rate
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
simulate.conc <- function(data) {
	#List of parameters from input for the differential equation solver
	PARAMETER.list <- c(ETA1 = data$ETA1,
											ETA2 = data$ETA2,
											ETA3 = data$ETA3,
											ETA4 = data$ETA4,
											ETA5 = data$ETA5,
											ETA6 = data$ETA6,
											ETA7 = data$ETA7,
											ETA8 = data$ETA8,
											CRCL = data$CRCL,
											FFM = data$FFM,
											HCT = data$HCT,
											SLC7A5 = data$SLC7A5,
											BUN = data$BUN,
											ANCBASE = data$ANCBASE,
											SEX = data$SEX,
											WBC = data$WBC,
											LNP53FOLD = data$LNP53FOLD,
											GCSF = data$GCSF,
											BSA = data$BSA,
											DOSE = data$DOSE)

	#Some individual parameters required for initial conditions
	BASEI <- POPBASE*exp(data$ETA4)*((data$ANCBASE/3.5)^ANC_BASE)*((data$HCT/32.5)^HCT_BASE)
	MTTI <- POPMTT*exp(data$ETA6)*(1+data$GCSF*GCSF_MTT)*((data$BUN/14)^BUN_MTT)

	ON <- 1
	if (data$GCSF == 1) {ON <- 0}
	IP0I <- (POPIP01*ON+(1-ON)*POPIP02)*exp(data$ETA7)*((data$LNP53FOLD/2.62)^LNP53_IP0)*((data$WBC/4.9)^FAC9)

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

	#Combine the unique PK and PD times (only one TIME sequence can be input to the differential equation solver)
	#Also, add in the time-point that signifies when the infusion has finished (INFD, if not already done so)
	INFD <- 0.5	#Fixed infusion duration of 30 minutes
	TIME <- sort(unique(c(PK.TIME,PD.TIME,INFD)))

	#Run differential equation solver for simulated variability data
	sim.data <- lsoda(y = A_0,TIME,func = DES.cmpf,parms = PARAMETER.list)
	sim.data <- as.data.frame(sim.data)
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset
simulate.conc.cmpf <- cmpfun(simulate.conc)

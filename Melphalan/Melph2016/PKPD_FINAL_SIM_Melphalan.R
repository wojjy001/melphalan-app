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
setwd("/Volumes/Prosecutor/Melphalan/Melph2016/")

#------------------------------------------------------------------------------------
#Read in the observed data
data <- read.csv("PKPD_OSU11055_Neutropenia_Boxcox_GCSF_lnf_f4_1.csv", stringsAsFactors=F, na.strings=".")
data <- subset(data, X.Note != "#")
data <- data[,-1]	#Delete the "Notes" column
data$TIME <- as.numeric(data$TIME)
data$ANC <- as.numeric(data$ANC)
data$DV <- as.numeric(data$DV)
data$lnp53fold[data$lnp53fold == -99] <- 2.62 #Original NONMEM control stream said to make the FCOV0 variable = 1, making this covariate value = 2.62 will make FCOV0 = 1
data$SLC7A5[is.na(data$SLC7A5)==T] <- 0
	
#Repeat data frame for x number of simulations
x <- 100
SIM <- seq(from = 1, to = x, by = 1)
SIM <- rep(SIM, times=length(data$ID))
SIM <- sort(SIM)
par.data <- lapply(data, rep.int, times=x)
par.data$SIM <- SIM
par.data <- as.data.frame(par.data)

#par.data <- par.data[par.data$ID == 119,]

#Identify PK sample times
PK.times <- sort(unique(par.data$TIME[par.data$DVID != 2]))
#Identify PD sample times
PD.times <- sort(unique(par.data$TIME[par.data$DVID == 2]))

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
#Calculate ETA values for each individual
ETA.data <- subset(par.data, DVID == 0)	#Faster than oneperID function, first line of each individual
ETA.data <- ETA.data[-c(2:20)]	#Remove dosing/concentration information
	
#Use nrorm to calculate ETA values
ETA.data$ETA1 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVCL_V2)
ETA.data$ETA2 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVV1)
ETA.data$ETA3 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVQ)
ETA.data$ETA4 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVBASE)
ETA.data$ETA5 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVSLOPE)
ETA.data$ETA6 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVMTT)
ETA.data$ETA7 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVIP0)
ETA.data$ETA8 <- rnorm(n = length(ETA.data$ID), mean = 0, sd = PPVIPT)
	
#Merge ETA.data with par.data
new.par.data <- merge(par.data,ETA.data,by=c("ID","SIM"),all=T)
		
#------------------------------------------------------------------------------------
#Make a reduced data frame to be used in differential equation solver
input.data <- subset(new.par.data, select = -c(ANC,DV,MDV,CMT,DVID))	#Remove the TIME, ANC, DV, MDV, CMT, DVID columns
input.data <- input.data[with(input.data, order(input.data$SIM,input.data$ID,input.data$TIME)), ]

#------------------------------------------------------------------------------------
#Function containing differential equations for amount in each compartment
DES <- function(T, A, PAR) {
	#9 differential equations
	#2 x PK compartments, 6 x PD compartments
	dAdT <- vector(length = 8)
	
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
	AMT <- PAR[19]	#Amount to be infused
	RATE <- PAR[20]	#Infusion rate
	INFD <- AMT/RATE	#Infusion duration
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
		
	list(dAdT)	#Get it to give the "amount" in each compartment
}

#Compile DES function - it's called by lsoda for each individual in the dataset	
DES.cmpf <- cmpfun(DES)

#------------------------------------------------------------------------------------
#Function for simulating concentrations for the ith patient
simulate.conc <- function(input.data) {		
	#Print a message stating which ID and SIM it is up to
	print(paste("ID:",input.data$ID[1],"SIM:",input.data$SIM[1],sep=" "))
	
	#Some individual parameters required for initial conditions
	BASEI <- POPBASE*exp(input.data$ETA4[1])*((input.data$ANCBASE[1]/3.5)^ANC_BASE)*((input.data$HCT[1]/32.5)^HCT_BASE)
	MTTI <- POPMTT*exp(input.data$ETA6[1])*(1+input.data$GCSF[1]*GCSF_MTT)*((input.data$BUN[1]/14)^BUN_MTT)
	
	ON <- 1
	if (input.data$GCSF[1] == 1) {ON <- 0}
	IP0I <- (POPIP01*ON+(1-ON)*POPIP02)*exp(input.data$ETA7[1])*((input.data$lnp53fold[1]/2.62)^LNP53_IP0)*((input.data$WBC[1]/4.9)^FAC9)
			
	#Set initial compartment conditions
	A1_0 <- 0	#Central compartment, PKCENTR (PK)
	A2_0 <- 0	#Peripheral compartment, PKPERI (PK)
	C3_0 <- (KE*BASEI)/(4/MTTI)	#STEM (PD)
	C4_0 <- BASEI	#ANC (PD)
	C5_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT1 (PD)
	C6_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT2 (PD)
	C7_0 <- (KE*BASEI)/(4/MTTI)	#TRANSIT3 (PD)
	C8_0 <- IP0I	#INPUT (PD)
	
	#Vector it!
	A_0 <- c(A1=A1_0, A2=A2_0, C3=C3_0, C4=C4_0, C5=C5_0, C6=C6_0, C7=C7_0, C8=C8_0)

	#Specify individual time sequences
	TIMEi <- input.data$TIME
	#Add in the time-points that signify when an individuals infusion has finished
	#Required for the differential equation solver
	INFDUR <- input.data$AMT[1]/input.data$RATE[1]
	TIMEs <- sort(unique(c(TIMEi,INFDUR)))
	
	#List of parameters from input for the differential equation solver			
	PARAMETER.list <- c(ETA1 = input.data$ETA1[1],
											ETA2 = input.data$ETA2[1],
											ETA3 = input.data$ETA3[1],
											ETA4 = input.data$ETA4[1],
											ETA5 = input.data$ETA5[1],
											ETA6 = input.data$ETA6[1],
											ETA7 = input.data$ETA7[1],
											ETA8 = input.data$ETA8[1],
											CRCL = input.data$CrCL[1],
											FFM = input.data$FFM[1],
											HCT = input.data$HCT[1],
											SLC7A5 = input.data$SLC7A5[1],
											BUN = input.data$BUN[1],
											ANCBASE = input.data$ANCBASE[1],
											SEX = input.data$SEX[1],
											WBC = input.data$WBC[1],
											LNP53FOLD = input.data$lnp53fold[1],
											GCSF = input.data$GCSF[1],
											AMT = input.data$AMT[1],
											RATE = input.data$RATE[1])	
	
	#Run differential equation solver for simulated variability data	
	var.data <- lsoda(y=A_0, TIMEs, func=DES.cmpf, parms=PARAMETER.list)
	var.data <- as.data.frame(var.data)
}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
simulate.conc.cmpf <- cmpfun(simulate.conc)

#Apply simulate.conc.cmpf to each individual in par.data
#Maintain their individual value for FFM and ETA2 for later calculations of V1
sim.data <- ddply(input.data, .(ID,SIM,FFM,ETA2), simulate.conc.cmpf)

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

# #------------------------------------------------------------------------------------
# #Plot - PK
# #Bin time - PK.data (same bins as NONMEM simulation for comparison!)
# PK.data <- sim.data
# names(PK.data)[5] <- "TIME"
# PK.data <- PK.data[PK.data$TIME %in% PK.times,]
# PK.data$TIMEBIN <- cut2(PK.data$TIME, g=12, levels.mean=T)
# PK.data$TIMEBIN <- as.numeric(paste(PK.data$TIMEBIN))

# #Bin time - data (original data)
# ORG.PK.data <- new.par.data[new.par.data$TIME %in% PK.times,]
# ORG.PK.data$TIMEBIN <- cut2(ORG.PK.data$TIME, g=12, levels.mean=T)
# ORG.PK.data$TIMEBIN <- as.numeric(paste(ORG.PK.data$TIMEBIN))

# #Function for calculating 5th and 95th percentiles for plotting concentrations
# CI90lo <- function(x) quantile(x, probs = 0.05)
# CI90hi <- function(x) quantile(x, probs = 0.95)

# #Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
# CI95lo <- function(x) quantile(x, probs = 0.025)
# CI95hi <- function(x) quantile(x, probs = 0.975)

# # PK.data.bystudy.median <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) median(df$PKDV)) 
# # PK.data.bystudy.median <- rename(PK.data.bystudy.median, c("V1"="medianS")) 

# # PK.data.bystudy.loCI <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) CI90lo(df$PKDV))
# # PK.data.bystudy.loCI <- rename(PK.data.bystudy.loCI, c("5%"="loCI90S")) 

# # PK.data.bystudy.hiCI <- ddply(new.PK.data, .(SIM,TIMEBIN), function(df) CI90hi(df$PKDV))
# # PK.data.bystudy.hiCI <- rename(PK.data.bystudy.hiCI, c("95%"="hiCI90S"))

# # PK.data.bystudy <- data.frame(PK.data.bystudy.median, "loCI90S"=PK.data.bystudy.loCI$loCI90S, "hiCI90S"=PK.data.bystudy.hiCI$hiCI90S)
	
# #Generate a plot of the PK.data (PK data)
# plotobj1 <- NULL
# plotobj1 <- ggplot()
# plotobj1 <- plotobj1 + geom_point(aes(x = TIME, y = DV), data = ORG.PK.data, colour = "blue", shape = 1, size = 2)

# # #Median simulated with confidence band
# # plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = medianS), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")

# # #Lower 90% CI simulated with confidence band
# # plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = loCI90S), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")

# # #Upper 90% CI simulated with confidence band
# # plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = hiCI90S), data = PK.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")

# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = ORG.PK.data, fun.y = median, geom = "line", colour = "red", size = 1)
# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = ORG.PK.data, fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = DV), data = ORG.PK.data, fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1)

# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = PKDV), data = PK.data, fun.y = median, geom = "line", size = 1)
# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = PKDV), data = PK.data, fun.y = "CI90lo", geom = "line", linetype = "dashed", size = 1)
# plotobj1 <- plotobj1 + stat_summary(aes(x = TIMEBIN, y = PKDV), data = PK.data, fun.y = "CI90hi", geom = "line", linetype = "dashed", size = 1)
	
# plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n")
# plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", lim = c(0,8), breaks = c(0,2,4,6,8))
# print(plotobj1)
	
#------------------------------------------------------------------------------------
#Plot - PD
#Bin time - PD.data (same bins as NONMEM simulation for comparison!)
PD.data <- sim.data
names(PD.data)[5] <- "TIME"
PD.data <- PD.data[PD.data$TIME %in% PD.times,]

#Pull out original PD data
ORG.PD.data <- new.par.data[new.par.data$TIME %in% PD.times,]

#Function for calculating 5th and 95th percentiles for plotting concentrations
CI90lo <- function(x) quantile(x, probs = 0.05)
CI90hi <- function(x) quantile(x, probs = 0.95)

#Function for calculating 2.5 and 97.5 percentiles for plotting concentrations
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

PD.data.bystudy.median <- ddply(PD.data, .(SIM,TIME), function(df) median(df$PDDV)) 
PD.data.bystudy.median <- rename(PD.data.bystudy.median, c("V1"="medianS")) 

PD.data.bystudy.loCI <- ddply(PD.data, .(SIM,TIME), function(df) CI90lo(df$PDDV))
PD.data.bystudy.loCI <- rename(PD.data.bystudy.loCI, c("5%"="loCI90S")) 

PD.data.bystudy.hiCI <- ddply(PD.data, .(SIM,TIME), function(df) CI90hi(df$PDDV))
PD.data.bystudy.hiCI <- rename(PD.data.bystudy.hiCI, c("95%"="hiCI90S"))

PD.data.bystudy <- data.frame(PD.data.bystudy.median, "loCI90S"=PD.data.bystudy.loCI$loCI90S, "hiCI90S"=PD.data.bystudy.hiCI$hiCI90S)
	
#Generate a plot of the PK.data (PK data)
plotobj2 <- NULL
plotobj2 <- ggplot()
plotobj2 <- plotobj2 + geom_point(aes(x = TIME, y = ANC), data = ORG.PD.data, colour = "blue", shape = 1, size = 2)

#Median simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = medianS), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")

#Lower 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = loCI90S), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")

#Upper 90% CI simulated with confidence band
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = hiCI90S), data = PD.data.bystudy, geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")

plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = ANC), data = ORG.PD.data, fun.y = median, geom = "line", colour = "red", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = ANC), data = ORG.PD.data, fun.y = "CI90lo", geom = "line", colour = "red", linetype = "dashed", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = ANC), data = ORG.PD.data, fun.y = "CI90hi", geom = "line", colour = "red", linetype = "dashed", size = 1)

plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = PDDV), data = PD.data, fun.y = median, geom = "line", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = PDDV), data = PD.data, fun.y = "CI90lo", geom = "line", linetype = "dashed", size = 1)
plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = PDDV), data = PD.data, fun.y = "CI90hi", geom = "line", linetype = "dashed", size = 1) 
	
plotobj2 <- plotobj2 + scale_y_log10("ANC (K/uL)\n", breaks = c(0.1,1,10))
plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
print(plotobj2)

#write.csv(sim.data, file="PKPD_FINAL_SIM_Melphalan_100_R.csv", na=".", quote=F, row.names=F)


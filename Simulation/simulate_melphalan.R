# R script for simulating a population from Yu Kyoung's melphalan-ANC model
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	#Plotting
	library(grid)	#Plotting
	library(dplyr)	#New plyr - required for mrgsolve
	library(mrgsolve)	#Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# ------------------------------------------------------------------------------
# Define time sequence - using mrgsolve's tgrid function
	TIME.tgrid <- c(tgrid(0.5,10,0.5),tgrid(24,720,24))
# Set number of individuals that make up the 95% prediction intervals
	n <- 50
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)
# Set seed for reproducible numbers
	# set.seed(123456)

# ------------------------------------------------------------------------------
# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# This compiled model is used for simulating n individuals and their concentration-time profiles
code <- '
  $CMT
          // Initial conditions for compartments
          CENT,	// Central compartment for PK
          PERI,	// Peripheral compartment for PK
          STEM,	// Stem compartment for PD
          ANC,	// Absolute neutrophil count for PD
          TRANSIT1,	// Transit compartment 1 for PD
          TRANSIT2,	// Transit compartment 2 for PD
          TRANSIT3, // Transit compartment 3 for PD
          INPUT, // Input compartment for PD
          G4N1	// Time spent in Grade 4 Neutropenia
  
  $PARAM
          // Population pharmacokinetic parameters
          POPCL = 28.9,	// Clearance, L/h (THETA1)
          POPV1 = 19.5,	// Volume of central compartment, L (THETA2)
          POPQ = 28.5,	// Inter-compartmental clearance, L/h (THETA3)
          POPV2 = 20.9,	// Volume of peripheral compartment, L (THETA4)
          SHARE_ETA = 1.31,	// ETA shared between CL and V2 (THETA5)
          W = 0.693,	// Standard deviation of residuals for PK (THETA9)
  
          // Population pharmacodynamic parameters
          POPBASE = 5.66, // (THETA10)
          POPSLOPE = 7.41,	// (THETA11)
          POPMTT = 96,	// Mean transit time, days (THETA12)
          POPTHALF = 7,	// (THETA15)
          POPIP01 = 109,	// (THETA16)
          POPIP02 = 0.0564,	// (THETA17)
          POPIPT = 14.7,	// (THETA18)
  
          // Covariate effects
          CRCL_CL = 0.193,	// Effect of CrCL on CL (THETA6)
          SLC7A5_V2 = -0.130,	// Effect of SLC7A5 on V2 (THETA7)
          HCT_CL = 0.216,	// Effect of HCT on CL (THETA8)
          ENDO_POWER1 = 0.188,	// Endogenous GSCF effect (THETA13)
          NEUP_POWER1 = 0.0313,	// Neupogen GSCF effect (THETA14)
          GCSF_MTT = 0.157,	// Effect of GCSF on MTT (THETA19)
          ANCBASE_BASE = 0.226,	// Effect of ANC on BASE (THETA20)
          GCSF_SLOPE = 0.421,	// Effect of GSCF on SLOPE (THETA21)
          RACE_MTT = -0.0894,	// Effect of RACE on MTT (THETA22)
          CRCL_MTT = -0.056,	// Effect of CrCL on MTT (THETA23)
          SEX_SLOPE = 0.214,	// Effect of SEX on SLOPE (THETA24)
          HCT_SLOPE = 0.646,	// Effect of HCT on SLOPE (THETA25)
          HCT_BASE = 0.582,	// Effect of HCT on BASE (THETA26)
  
          // Default covariate values for simulation
          CRCL = 91.94,	// Creatinine clearance (mL/min)
          HCT = 32.5,	// Haematocrit (%)
          SEX = 0,	// Male (1) or female (0)
          FFM = 59.9,	// Fat free mass (kg)
          SLC7A5 = 0,	// AA or AG (0) versus GG (1)
          GCSF = 0,	// Neupogen administered on Day 1 (0) or Day 7 (1)
          RACE = 0,	// Caucasian or unknown (0) versus African-American (1)
          ANCBASE = 3.5	// Baseline ANC (K/µL)
  
  $OMEGA
          name = "BSV"
          block = FALSE
          labels = s(BSVCL,BSVV1,BSVQ,BSVBASE,BSVSLOPE,BSVMTT,BSVIP0)
          0.0792	// BSVCL
          0.0740	// BSVV1
          0.297	// BSVQ
          0.1	// BSVBASE
          0.0631	// BSVSLOPE
          0.00561	// BSVMTT
          0.23	// BSVIP0
  
  $SIGMA
          block = FALSE
          labels = s(ERRPK,ERRPD)
          0.0584	// Log-normal proportional error for PK observations
          1	// Additive error for PD observations
  
  $MAIN
          //////////////////////
          // PHARMACOKINETICS //
          //////////////////////
  
          // Infusion duration
          D_CENT = 0.5;	// 30 minutes
          
          // Covariate effects
          double SLC7A5COV_V2 = 1;
          if (SLC7A5 == 1) SLC7A5COV_V2 = 1 + SLC7A5_V2;
          
          // Individual parameter values
          if (ID == 1) {
          BSVCL = 0;
          BSVV1 = 0;
          BSVQ = 0;
          }
          double CL = POPCL*pow(FFM/59.9,0.75)*pow(CRCL/91.94,CRCL_CL)*pow(HCT/32.5,HCT_CL)*exp(BSVCL);
          double V1 = POPV1*(FFM/59.9)*exp(BSVV1);
          double Q = POPQ*pow(FFM/59.9,0.75)*exp(BSVQ);
          double V2 = POPV2*(FFM/59.9)*SLC7A5COV_V2*exp(BSVCL*SHARE_ETA);
  
          // Individual micro-rate constants
          double K10 = CL/V1;
          double K12 = Q/V1;
          double K21 = Q/V2;
          
          // Initial conditions
          CENT_0 = 0;
          PERI_0 = 0;
  
          //////////////////////
          // PHARMACODYNAMICS //
          //////////////////////
          
          // Covariate effects
          double GCSFCOV_SLOPE = 1;
          if (GCSF == 1) GCSFCOV_SLOPE = 1 + GCSF_SLOPE;
          double GCSFCOV_MTT = 1;
          if (GCSF == 1) GCSFCOV_MTT = 1 + GCSF_MTT;
          double SEXCOV_SLOPE = 1;
          if (SEX == 1) SEXCOV_SLOPE = 1 + SEX_SLOPE;
          double RACECOV_MTT = 1;
          if (RACE == 1) RACECOV_MTT = 1 + RACE_MTT;
          double ON_IP0 = 1;
          if (GCSF == 1) ON_IP0 = 0;
          
          // Population parameter values
          if (ID == 1) {
          BSVBASE = 0;
          BSVSLOPE = 0;
          BSVMTT = 0;
          BSVIP0 = 0;
          }
          double BASE = POPBASE*pow(ANCBASE/3.5,ANCBASE_BASE)*pow(HCT/32.5,HCT_BASE)*exp(BSVBASE);
          double SLOPE = POPSLOPE*GCSFCOV_SLOPE*SEXCOV_SLOPE*pow(HCT/32.5,HCT_SLOPE)*exp(BSVSLOPE);
          double MTT = POPMTT*GCSFCOV_MTT*RACECOV_MTT*pow(CRCL/91.94,CRCL_MTT)*exp(BSVMTT);
          double KE = log(2)/POPTHALF;
          double IP0 = ON_IP0*POPIP01+(1-ON_IP0)*POPIP02*exp(BSVIP0);
          double IPT = POPIPT;
          
          // Individual rate constants
          double K = 4/MTT;
          
          // Initial conditions
          STEM_0 = KE*BASE/K;
          ANC_0 = BASE;
          TRANSIT1_0 = KE*BASE/K;
          TRANSIT2_0 = KE*BASE/K;
          TRANSIT3_0 = KE*BASE/K;
          INPUT_0 = IP0;
          G4N1_0 = 0;
  
  $ODE
          //////////////////////
          // PHARMACOKINETICS //
          //////////////////////
          
          // Differential equations
          dxdt_CENT = K21*PERI -K10*CENT -K12*CENT;
          dxdt_PERI = K12*CENT -K21*PERI;
          
          //////////////////////
          // PHARMACODYNAMICS //
          //////////////////////
          
          // Other variables required for differential equations
          double QN = 0;
          if ((GCSF == 0) & (SOLVERTIME > 72)) QN = 1;
          if ((GCSF == 1) & (SOLVERTIME > 216)) QN = 1;
          double POWER1 = ENDO_POWER1 + QN*NEUP_POWER1;
          double FN = pow(BASE/0.001,POWER1);
          if (ANC > 0.001) FN = pow(BASE/ANC,POWER1);
          
          double KIN = 0;
          if ((GCSF == 0) & (SOLVERTIME > 72)) KIN = 1/IPT;
          if ((GCSF == 1) & (SOLVERTIME > 216)) KIN = 1/IPT;
          
          double CP = CENT/V1;
          double DRUG = SLOPE*CP;
  
          // Differential equations
          dxdt_STEM = K*STEM*(1-DRUG)*FN -K*STEM;
          dxdt_ANC = K*TRANSIT3 -KE*ANC +KIN*INPUT;
          dxdt_TRANSIT1 = K*STEM -K*TRANSIT1;
          dxdt_TRANSIT2 = K*TRANSIT1 -K*TRANSIT2;
          dxdt_TRANSIT3 = K*TRANSIT2 -K*TRANSIT3;
          dxdt_INPUT = -KIN*INPUT;
          dxdt_G4N1 = 0;
          if (ANC < 0.5) dxdt_G4N1 = 1;
  
  $TABLE 
          table(IPRE) = CENT/V1;
'
	# Compile the model code
		mod <- mcode("popMELPH",code)
		# There is opportunity to simply update model parameters after the model code has been compiled

# ------------------------------------------------------------------------------
# Simulate concentration-time profiles for the population
	input.conc.data <- expand.ev(
		ID = 1:n,	# n individuals
		amt = 100*1.8,	# amt in mg
		evid = 1,	# evid = 1; dosing event
		cmt = 1,	# cmt = 1; dose goes into compartment 1 = cent
		time = 0,	# time = 0; dose at time = 0
		rate = -2
	)
	conc.data <- mod %>% data_set(input.conc.data) %>% mrgsim(tgrid = TIME.tgrid)
	conc.data <- as.data.frame(conc.data)	#Convert to a data frame so that it is more useful for me!

# ------------------------------------------------------------------------------
# Plot results - melphalan concentrations
	plotobj1 <- NULL
	plotobj1 <- ggplot(conc.data)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),data = conc.data[conc.data$ID == 1,],geom = "line",fun.y = median,colour = "red")
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),data = conc.data[conc.data$ID != 1,],geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)",lim = c(0,10))
	plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n",breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
	print(plotobj1)

# Plot results - ANC
	plotobj2 <- NULL
	plotobj2 <- ggplot(conc.data)
	plotobj2 <- plotobj2 + stat_summary(aes(x = time,y = ANC),geom = "line",fun.y = median,colour = "red")
	plotobj2 <- plotobj2 + stat_summary(aes(x = time,y = ANC),geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.5),linetype = "dashed")
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
	plotobj2 <- plotobj2 + scale_y_log10("Melphalan Concentration (mg/L)\n",breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
	print(plotobj2)

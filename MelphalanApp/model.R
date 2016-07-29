# Define the model parameters and equations
# Using mrgsolve - differential equations
# This compiled model is used for simulating n individuals and their concentration-time profiles
code <- '
	$INIT	// Initial conditions for compartments
		DEPOT = 0, // Depot - dose is added here
		TRANS1 = 0, // Transit 1
		TRANS2 = 0, // Transit 2
		CENT = 0, // Central
		PERI = 0, // Peripheral
		AUC = 0 // Area under the curve compartment

	$PARAM	// Population parameters
		POPCL = 4.63, //THETA4
		POPV = 55.2, //THETA3
		POPCLP1 = 11.3, //THETA6
		POPVP1 = 49.8, //THETA5
		POPKTR = 2, //THETA7
		POPF = 1,

		// Covariate effects
		F1XC = 0.863, //THETA8
		F1CAP = 0.978, //THETA9
		COVFED = -0.209, //THETA10
		COVFED2 = -0.549, //THETA11
		ALAG1 = 0.115, //THETA12
		FTLAG2 = 0.203, //THETA13
		COVFEDF = 0.105, //THETA14
		COVSEX = 0.144, //THETA15

		// Averaged values for study effects on parameters
		COVSTDF = 1.157, //Study averaged F for simulation
		COVSTDKTR = 0.894, //Study averaged KTR for simulation
		COVSTDCL = 0.789, //Study averaged CL for simulation
		COVSTDV = 0.786, //Study averaged V for simulation

		// Default covariate values for simulation
		FED = 1, // Fed (1) or fasted (0)
		SEX = 1, // Male (1) or female (0)
		FFM = 55.49, // Fat free mass (kg)
		TRT = 1, // Formulation; Doryx MPC (1), Doryx tablet (2), Doryx capsule (3)
		PER = 1	 // Period; first occassion (1) or second occassion (2)

	$OMEGA	name = "BSV"
		block = TRUE
		labels = s(BSV_CL,BSV_KTR,BSV_VP1,BSV_V)
		0.0373 // BSV for CL
		0.0229 0.0796 // BSV for KTR
		0.0106 -0.01 0.0229 // BSV for VP1
		0.0522 0.0936 -0.00506 0.141 // BSV for V

	$OMEGA	name = "BOV"
		block = FALSE
		labels = s(BOV_CL1,BOV_CL2,BOV_KTR1,BOV_KTR2)
		0.0183 // BOV on first occassion for CL
		0.0183 // BOV on second occassion for CL
		0.0738 // BOV on first occassion for KTR
		0.0738 // BOV on second occassion for KTR

	$SIGMA	block = FALSE
		labels = s(ERR_PRO,ERR_ADD)
		0.0384 // Proportional error
		392.04 // Additive error

	$MAIN	// Covariate effects
		double FEDCOV2 = 1; // Fasted
		if (FED == 1) FEDCOV2 = 1+COVFEDF; // Fed
		double SEXCOV = 1; // Male
		if (SEX == 0) SEXCOV = 1+COVSEX; // Female
		double FEDCOV = 1; // Fasted
		if (FED == 1 & TRT == 1) FEDCOV = 1+COVFED2; // Fed and Doryx MPC
		if (FED == 1 & TRT == 2) FEDCOV = 1+COVFED; // Fed and Doryx tablet
		if (FED == 1 & TRT == 3) FEDCOV = 1+COVFED; // Fed and Doryx capsule

		// Between-occassion variability
			// Clearance
			double BOV_CL = BOV_CL1;
			if (PER == 2) BOV_CL = BOV_CL2;
			double ETA_CL = BSV_CL+BOV_CL;
			// Transit
			double BOV_KTR = BOV_KTR1;
			if (PER == 2) BOV_KTR = BOV_KTR2;
			double ETA_KTR = BSV_KTR+BOV_KTR;
			// Volume - peripheral
			double ETA_VP1 = BSV_VP1;
			// Volume - central
			double ETA_V = BSV_V;

		// Individual parameter values
		double CL = POPCL*pow(FFM/70,0.75)*exp(ETA_CL)*COVSTDF*COVSTDCL*FEDCOV2*SEXCOV;
		double V = POPV*(FFM/70)*exp(ETA_V)*COVSTDF*COVSTDV*FEDCOV2;
		double CLP1 = POPCLP1*pow(FFM/70,0.75)*COVSTDF*FEDCOV2;
		double VP1 = POPVP1*(FFM/70)*exp(ETA_VP1)*COVSTDF*FEDCOV2;
		double KTR = POPKTR*exp(ETA_KTR)*FEDCOV*COVSTDKTR;
		double F = POPF; // Doryx tablet
		if (TRT == 1) F = F1XC; // Doryx MPC
		if (TRT == 3) F = F1CAP; // Doryx capsule
		F_DEPOT = F;

		// Micro-rate constants
		double K12 = KTR; //DEPOT to TRANS1
		double K23 = KTR; //TRANS1 to TRANS2
		double K34 = KTR; //TRANS2 to CENT
		double K45 = CLP1/V; //CENT to PERI
		double K54 = CLP1/VP1; //PERI to CENT
		double K46 = CL/V; //CENT to elimination

		// Lag time
		double TLAG = 0; // Doryx tablet
		if (TRT == 1) TLAG = ALAG1; // Doryx MPC
		if (TRT == 3) TLAG = ALAG1; // Doryx capsule
		double FTLAG = 0; // Fasted
		if (FED == 1) FTLAG = FTLAG2; // Fed
		double ALAG_DEPOT = TLAG+FTLAG;

	$ODE	// Differential equations
		dxdt_DEPOT = -K12*DEPOT;
		dxdt_TRANS1 = K12*DEPOT -K23*TRANS1;
		dxdt_TRANS2 = K23*TRANS1 -K34*TRANS2;
		dxdt_CENT = K34*TRANS2 -K45*CENT +K54*PERI -K46*CENT;
		dxdt_PERI = K45*CENT -K54*PERI;
		dxdt_AUC = CENT/V;

		// Cmax and Tmax
		double CP = CENT/V;
		if (SOLVERTIME == 0) double Cmax = 0;
		if (SOLVERTIME == 0) double Tmax = 0;
		if (CP > Cmax) {
			Cmax = CP;
			Tmax = SOLVERTIME;
		}

	$TABLE	table(IPRE) = CENT/V;
		table(DV) = table(IPRE)*(1+ERR_PRO)+ERR_ADD;

	$CAPTURE CL V CLP1 VP1 KTR F ETA_CL ETA_V ETA_VP1 ETA_KTR FED SEX FFM TRT PER Cmax Tmax
'
# Compile the model code
mod <- mcode("popDOXY",code) # There is opportunity to simply update model parameters after the model code has been compiled

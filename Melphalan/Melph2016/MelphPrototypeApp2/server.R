#Define server for MelphPrototypeApp1
#Calls on reactive user-defined input from ui.R
#Calls on non-reactive user-defined functions and objects from global.R
#-----------------------------------------------------------------------------------
#Define user-input dependent functions for output
shinyServer(function(input, output, session) {
	###########
	##_INPUT_##
	###########
	#Rinput.data calls in user-defined widget input and makes a data frame
	Rinput.data <- reactive({
		#Call in user-defined widget values
		AGE <- input$AGE	#Numeric input for patient's age
		TBW <- input$TBW	#Numeric input for patient's total body weight
		HT <- input$HT	#Numeric input for patient's height
		SECR <- input$SECR	#Numeric input for patient's serum creatinine
		HCT <- input$HCT	#Numeric input for patient's haematocrit
		BUN <- input$BUN	#Numeric input for patient's blood urea nitrogen
		WBC <- input$WBC	#Numeric input for patient's white blood cell count
		ANCBASE <- input$ANC	#Numeric input for patient's absolute neutrophil count
		if (input$SEX == 1) {SEX <- 0}	#Select input for patient's gender, Female = 0
		if (input$SEX == 2) {SEX <- 1}	#Select input for patient's gender, Male = 1
		if (input$SLC7A5 == 1) {SLC7A5 <- 0}	#Select input for patient's SLC7A5
		if (input$SLC7A5 == 2) {SLC7A5 <- 1}	#Select input for patient's SLC7A5
		LNP53FOLD <- input$LNP53FOLD	#Numeric input for patient's LNP53FOLD
		#Calculate body mass index (BMI), body surface area (BSA), creatinine clearance (CRCL), and fat free mass (FFM) based on the above input
		BMI <- TBW/(HT/100)^2
		BSA <- sqrt((HT*TBW)/3600)
		if (SEX == 0) {CRCL <- (((140-AGE)*TBW)/(SECR*0.815))*0.85}	#Creatinine clearance (mL/min) for females
		if (SEX == 1)	{CRCL <- ((140-AGE)*TBW)/(SECR*0.815)}	#Creatinine clearance (mL/min) for males
		if (SEX == 0) {FFM <- 9270*TBW/(8780+(244*BMI))}	#Fat free mass (kg) for females
		if (SEX == 1)	{FFM <- 9270*TBW/(6680+(216*BMI))}	#Fat free mass (kg) for males
		#Make a data frame of widget values to be input into differential equation solver functions
		#Simulate 1 individual for PRED (population predicted) calculation, set ETAs to 0
		pred.data <- data.frame(ID = 1,
														ETA1 = 0,	#PPV for CL_V2
														ETA2 = 0,	#PPV for V1
														ETA3 = 0,	#PPV for Q
														ETA4 = 0,	#PPV for BASE
														ETA5 = 0,	#PPV for SLOPE
														ETA6 = 0,	#PPV for MTT
														ETA7 = 0,	#PPV for IP0
														ETA8 = 0,	#PPV for IPT
														BMI,
														BSA,
														CRCL,
														FFM,
														HCT,
														SLC7A5,
														BUN,
														ANCBASE,
														SEX,
														WBC,
														LNP53FOLD)
		#Make a data frame of widget values to be input into differential equation solver functions
		#Contains multiple individuals to calculate 95% prediction intervals
		var.data <- data.frame(ID = 2:(n+1), #Simulate n individuals
														ETA1,	#PPV for CL_V2
														ETA2,	#PPV for V1
														ETA3,	#PPV for Q
														ETA4,	#PPV for BASE
														ETA5,	#PPV for SLOPE
														ETA6,	#PPV for MTT
														ETA7,	#PPV for IP0
														ETA8,	#PPV for IPT
														BMI,
														BSA,
														CRCL,
														FFM,
														HCT,
														SLC7A5,
														BUN,
														ANCBASE,
														SEX,
														WBC,
														LNP53FOLD)
		input.data <- rbind(pred.data,var.data)	#Create "input.data," combination of PRED and variability individuals
		as.data.frame(input.data)
	})	#Brackets closing "reactive" expression for "Rinput.data"
	#Rsim.data1 simulates PK and PD profiles based on Rinput.data and DOSE1 and GCSF1
	Rsim.data1 <- reactive({
		withProgress(
			message = "Updating Regimen 1",
			value = 0,
			{
				Rsim.data1.function <- function(input.data) {
					if (input$GCSF1 == 1) {input.data$GCSF <- 0}	#Select input for when to administer G-CSF (Neupogen), Day 1 = 0
					if (input$GCSF1 == 2) {input.data$GCSF <- 1}	#Select input for when to administer G-CSF (Neupogen), Day 7 = 1
					input.data$DOSE <- input$DOSE1	#Slider input for Melphalan dose
					#Simulate
					if (input$PI951 == FALSE) {input.data <- input.data[input.data$ID == 1,]}
					sim.data1 <- ddply(input.data, .(ID,FFM,ETA2), simulate.conc.cmpf)
					sim.data1 <- as.data.frame(sim.data1)
				}
				sim.data1 <- Rsim.data1.function(Rinput.data())
				sim.data1 <- as.data.frame(sim.data1)
			}	#Brackets closing expression for "withProgress"
		)	#Brackets closing "withProgress"
	})	#Brackets closing "reactive" expression for "Rsim.data1"
	#Rsim.data2 simulates PK and PD profiles based on Rinput.data and DOSE2 and GCSF2
	Rsim.data2 <- reactive({
		withProgress(
			message = "Updating Regimen 2",
			value = 0,
			{
				Rsim.data2.function <- function(input.data) {
					if (input$GCSF2 == 1) {input.data$GCSF <- 0}	#Select input for when to administer G-CSF (Neupogen), Day 1 = 0
					if (input$GCSF2 == 2) {input.data$GCSF <- 1}	#Select input for when to administer G-CSF (Neupogen), Day 7 = 1
					input.data$DOSE <- input$DOSE2	#Slider input for Melphalan dose
					#Simulate
					if (input$PI952 == FALSE) {input.data <- input.data[input.data$ID == 1,]}
					sim.data2 <- ddply(input.data, .(ID,FFM,ETA2), simulate.conc.cmpf)
					sim.data2 <- as.data.frame(sim.data2)
				}
				sim.data2 <- Rsim.data2.function(Rinput.data())
				sim.data2 <- as.data.frame(sim.data2)
			}	#Brackets closing expression for "withProgress"
		)	#Brackets closing "withProgress"
	})	#Brackets closing "reactive" expression for "Rsim.data2"
	#Rsim.data3 simulates PK and PD profiles based on Rinput.data and DOSE2 and GCSF2
	Rsim.data3 <- reactive({
		withProgress(
			message = "Updating Regimen 3",
			value = 0,
			{
				Rsim.data3.function <- function(input.data) {
					if (input$GCSF3 == 1) {input.data$GCSF <- 0}	#Select input for when to administer G-CSF (Neupogen), Day 1 = 0
					if (input$GCSF3 == 2) {input.data$GCSF <- 1}	#Select input for when to administer G-CSF (Neupogen), Day 7 = 1
					input.data$DOSE <- input$DOSE3	#Slider input for Melphalan dose
					#Simulate
					if (input$PI953 == FALSE) {input.data <- input.data[input.data$ID == 1,]}
					sim.data3 <- ddply(input.data, .(ID,FFM,ETA2), simulate.conc.cmpf)
					sim.data3 <- as.data.frame(sim.data3)
				}
				sim.data3 <- Rsim.data3.function(Rinput.data())
				sim.data3 <- as.data.frame(sim.data3)
			}	#Brackets closing expression for "withProgress"
		)	#Brackets closing "withProgress"
	})	#Brackets closing "reactive" expression for "Rsim.data3"
	############
	##_OUTPUT_##
	############
	output$outputBMIText <- renderUI({	#Send calculated BMI to ui.R
		outputBMIText.function <- function(input.data) {
			withMathJax(
				helpText(
					"Body mass index = ",round(input.data$BMI[1],digits = 1)," \\(kg/m^2\\)"
				)	#Brackets closing "helpText"
			)	#Brackets closing "withMathJax"
		}
		outputBMIText.function(Rinput.data())
	})	#Brackets closing "renderUI" expression for "outputBMIText"
	output$outputBSAText <- renderUI({	#Send calculate BSA to ui.R
		outputBSAText.function <- function(input.data) {
			withMathJax(
				helpText(
					"Body surface area = ",round(input.data$BSA[1],digits = 1)," \\(m^2\\)"
				)	#Brackets closing "helpText"
			)	#Brackets closing "withMathJax"
		}
		outputBSAText.function(Rinput.data())
	})	#Brackets closing "renderUI" expression for "outputBSAText"
	output$outputFFMText <- renderUI({	#Send calculate FFM to ui.R
		outputFFMText.function <- function(input.data) {
			withMathJax(
				helpText(
					"Fat free mass = ",round(input.data$FFM[1],digits = 1)," \\(kg\\)"
				)	#Brackets closing "helpText"
			)	#Brackets closing "withMathJax"
		}
		outputFFMText.function(Rinput.data())
	})	#Brackets closing "renderUI" expression for "outputFFMText"
	output$outputCRCLText <- renderUI({	#Send calculated CRCL to ui.R
		outputCRCLText.function <- function(input.data) {
			withMathJax(
				helpText(
					"Creatinine clearance = ",round(input.data$CRCL[1],digits = 1)," \\(mL/min\\)"
				)	#Brackets closing "helpText"
			)	#Brackets closing "withMathJax"
		}
		outputCRCLText.function(Rinput.data())
	})	#Brackets closing "renderUI" expression for "outputCRCLText"
	#Plot simulated PD profile (ANC)
	output$outputPDPlot <- renderPlot({
		outputPDPlot.function <- function(sim.data1,sim.data2,sim.data3) {
			plotobj <- ggplot()
			#Show population predicted PD profile as a solid line
			#C4 is the "amount" in compartment 4, i.e., PD dependent variable
			plotobj <- plotobj + geom_hline(aes(yintercept = 0.5),linetype = "dashed") #Grade 4 neutropenia
			if (input$TUNITS == FALSE) {
				plotobj <- plotobj + geom_line(aes(x = time/24,y = C4),data = sim.data1[sim.data1$ID == 1,],colour = "red")	#Regimen 1
				if (input$NREG > 1) {	#Regimen 2
					plotobj <- plotobj + geom_line(aes(x = time/24,y = C4),data = sim.data2[sim.data2$ID == 1,],colour = "blue")
				}
				if (input$NREG > 2) {	#Regimen 3
					plotobj <- plotobj + geom_line(aes(x = time/24,y = C4),data = sim.data3[sim.data3$ID == 1,],colour = "darkgreen")
				}
				#Show 95% prediction intervals as shaded ribbons
				if (input$PI951 == TRUE) {	#Regimen 1
					plotobj <- plotobj + stat_summary(aes(x = time/24,y = C4),data = sim.data1[sim.data1$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
				}
				if (input$PI952 == TRUE & input$NREG > 1) {	#Regimen 2
					plotobj <- plotobj + stat_summary(aes(x = time/24,y = C4),data = sim.data2[sim.data2$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "blue",alpha = 0.3)
				}
				if (input$PI953 == TRUE & input$NREG > 2) {	#Regimen 3
					plotobj <- plotobj + stat_summary(aes(x = time/24,y = C4),data = sim.data3[sim.data3$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "darkgreen",alpha = 0.3)
				}
				plotobj <- plotobj + scale_x_continuous("\nTime Since Melphalan Dose (days)",breaks = seq(from = 0,to = max(PD.TIME)/24,by = 5))
				plotobj <- plotobj + annotate("text",x = 27,y = 0.8,label = "Grade 4 Neutropenia",size = 5)
			}
			if (input$TUNITS == TRUE) {	#Change time-axis units
				plotobj <- plotobj + geom_line(aes(x = time,y = C4),data = sim.data1[sim.data1$ID == 1,],colour = "red")	#Regimen 1
				if (input$NREG > 1) {	#Regimen 2
					plotobj <- plotobj + geom_line(aes(x = time,y = C4),data = sim.data2[sim.data2$ID == 1,],colour = "blue")
				}
				if (input$NREG > 2) {	#Regimen 3
					plotobj <- plotobj + geom_line(aes(x = time,y = C4),data = sim.data3[sim.data3$ID == 1,],colour = "darkgreen")
				}
				#Show 95% prediction intervals as shaded ribbons
				if (input$PI951 == TRUE) {	#Regimen 1
					plotobj <- plotobj + stat_summary(aes(x = time,y = C4),data = sim.data1[sim.data1$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
				}
				if (input$PI952 == TRUE & input$NREG > 1) {	#Regimen 2
				plotobj <- plotobj + stat_summary(aes(x = time,y = C4),data = sim.data2[sim.data2$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "blue",alpha = 0.3)
				}
				if (input$PI953 == TRUE & input$NREG > 2) {	#Regimen 3
					plotobj <- plotobj + stat_summary(aes(x = time,y = C4),data = sim.data3[sim.data3$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "darkgreen",alpha = 0.3)
				}
				plotobj <- plotobj + scale_x_continuous("\nTime Since Melphalan Dose (hours)",breaks = seq(from = 0,to = max(PD.TIME),by = 100))	
				plotobj <- plotobj + annotate("text",x = 648,y = 0.8,label = "Grade 4 Neutropenia",size = 5)				
			}
			plotobj <- plotobj + scale_y_log10("Absolute Neutrophil Count (K/µL)\n",breaks = c(0.1,1,10),lim = c(0.01,50))
			print(plotobj)
		}
		outputPDPlot.function(Rsim.data1(),Rsim.data2(),Rsim.data3())
	})	#Brackets closing "renderPlot" expression for "outputPDPlot"
	#Plot simulated PK profile (melphalan concentration)
	output$outputPKPlot <- renderPlot({
		outputPKPlot.function <- function(sim.data1,sim.data2,sim.data3) {
			plotobj <- ggplot()
			#Show population predicted PK profile as a solid line
			#A1 is the "amount" in compartment 1, i.e., PK dependent variable
			plotobj <- plotobj + geom_line(aes(x = time,y = A1/POPV1*(FFM/59.9)),data = sim.data1[sim.data1$ID == 1,],colour = "red")	#Regimen 1
			if (input$NREG > 1) {	#Regimen 2
			plotobj <- plotobj + geom_line(aes(x = time,y = A1/POPV1*(FFM/59.9)),data = sim.data2[sim.data2$ID == 1,],colour = "blue")
			}
			if (input$NREG > 2) {	#Regimen 3
			plotobj <- plotobj + geom_line(aes(x = time,y = A1/POPV1*(FFM/59.9)),data = sim.data3[sim.data3$ID == 1,],colour = "darkgreen")
			}
			#Show 95% prediction intervals as shaded ribbons
			if (input$PI951 == TRUE) {	#Regimen 1
			plotobj <- plotobj + stat_summary(aes(x = time,y = A1/POPV1*(FFM/59.9)*exp(ETA2)),data = sim.data1[sim.data1$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "red",alpha = 0.3)
			}
			if (input$PI952 == TRUE & input$NREG > 1) {	#Regimen 2
			plotobj <- plotobj + stat_summary(aes(x = time,y = A1/POPV1*(FFM/59.9)*exp(ETA2)),data = sim.data2[sim.data2$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "blue",alpha = 0.3)
			}
			if (input$PI953 == TRUE & input$NREG > 2) {	#Regmen 3
			plotobj <- plotobj + stat_summary(aes(x = time,y = A1/POPV1*(FFM/59.9)*exp(ETA2)),data = sim.data3[sim.data3$ID > 1,],geom = "ribbon",fun.ymin = "CI95lo",fun.ymax = "CI95hi",fill = "darkgreen",alpha = 0.3)
			}
			plotobj <- plotobj + scale_x_continuous("\nTime (hours)",lim = c(0,max(PK.TIME)))
			plotobj <- plotobj + scale_y_log10("Melphalan Concentration (mg/L)\n",breaks = c(0.1,1,10),lim = c(0.01,30))
			print(plotobj)
		}
		outputPKPlot.function(Rsim.data1(),Rsim.data2(),Rsim.data3())
	})	#Brackets closing "renderPlot" expression for "outputPKPlot"
	#Plot simulated duration in Grade 4 Neutropenia
	output$outputG4NPlot <- renderPlot({
		outputG4NPlot.function <- function(sim.data1,sim.data2,sim.data3) {
			sim.data1.tail <- ddply(sim.data1, .(ID), lastperID)	#Take the last line of each individual
			sim.data2.tail <- ddply(sim.data2, .(ID), lastperID)
			sim.data3.tail <- ddply(sim.data3, .(ID), lastperID)
			plotobj <- ggplot()
			#Show population predicted duration in Grade 4 Neutropenia as points
			#S9 is the overall duration spent in severe neutropenia
			plotobj <- plotobj + geom_point(aes(x = input$DOSE1,y = S9),data = sim.data1.tail[sim.data1.tail$ID == 1,],colour = "red",size = 3)	#Regimen 1
			if (input$NREG > 1) {	#Regimen 2
				plotobj <- plotobj + geom_point(aes(x = input$DOSE2,y = S9),data = sim.data2.tail[sim.data2.tail$ID == 1,],colour = "blue",size = 3)
			}
			if (input$NREG > 2) {	#Regimen 3
				plotobj <- plotobj + geom_point(aes(x = input$DOSE3,y = S9),data = sim.data3.tail[sim.data3.tail$ID == 1,],colour = "darkgreen",size = 3)
			}
			#Show error bars if "Show 95% Prediction Intervals" checkbox is ticked
			if (input$PI951 == TRUE) {	#Regimen 1
				plotobj <- plotobj + geom_errorbar(aes(x = input$DOSE1,ymax = CI95hi(S9),ymin = CI95lo(S9)),data = sim.data1.tail[sim.data1.tail$ID > 1,],colour = "red",width = 5)
			}
			if (input$PI952 == TRUE & input$NREG > 1) {	#Regimen 2
				plotobj <- plotobj + geom_errorbar(aes(x = input$DOSE2,ymax = CI95hi(S9),ymin = CI95lo(S9)),data = sim.data2.tail[sim.data2.tail$ID > 1,],colour = "blue",width = 5)
			}
			if (input$PI953 == TRUE & input$NREG > 2) {	#Regimen 3
				plotobj <- plotobj + geom_errorbar(aes(x = input$DOSE3,ymax = CI95hi(S9),ymin = CI95lo(S9)),data = sim.data3.tail[sim.data3.tail$ID > 1,],colour = "darkgreen", width = 5)
			}
			plotobj <- plotobj + xlab(expression(paste("Melphalan Dose (",mg/m^2,")",sep="")))
			plotobj <- plotobj + scale_y_continuous("Duration of Severe Neutropenia (hours)\n")
			print(plotobj)
		}
		outputG4NPlot.function(Rsim.data1(),Rsim.data2(),Rsim.data3())
	})	#Brackets closing "renderPlot" expression for "outputG4NPlot"
	#Send summary of duration in Grade 4 neutropenia to ui.R for regimen 1
	output$outputG4N1Text <- renderText({
		outputG4N1Text.function <- function(sim.data1) {
			sim.data1.tail <- ddply(sim.data1, .(ID), lastperID)	#Find the last row for each individual
			pop.G4N1 <- round(sim.data1.tail$S9[sim.data1.tail$ID == 1])	#Population predicted duration in Grade 4 Neutropenia
			outputG4N1Text <- paste("Duration of Severe Neutropenia = ",pop.G4N1," hours",sep="")	#Population predicted summary text
			if (input$PI951 == TRUE) {	#If "Show 95% Prediction Intervals" for Regimen 1 is selected, calculate 95% prediction intervals 
				lo.G4N1 <- round(mean(sim.data1.tail$S9[sim.data1.tail$ID > 1]))	#2.5th percentile
				hi.G4N1 <- round(CI95hi(sim.data1.tail$S9[sim.data1.tail$ID > 1]))	#97.5th percentile
				outputG4N1Text <- paste("Duration of Severe Neutropenia = ",pop.G4N1," hours (",lo.G4N1," - ",hi.G4N1,")",sep="")	#Population predicted and 95% prediction intervals summary text
			}
			outputG4N1Text	#Print text
		}
		outputG4N1Text.function(Rsim.data1())
	})	#Brackets closing "renderText" expression for "outputG4N1Text"
	#Send summary of duration in Grade 4 neutropenia to ui.R for regimen 2
	output$outputG4N2Text <- renderText({
		outputG4N2Text.function <- function(sim.data2) {
			sim.data2.tail <- ddply(sim.data2, .(ID), lastperID)	#Find the last row for each individual
			pop.G4N2 <- round(sim.data2.tail$S9[sim.data2.tail$ID == 1])	#Population predicted duration in Grade 4 Neutropenia
			outputG4N2Text <- paste("Duration of Severe Neutropenia = ",pop.G4N2," hours",sep="")	#Population predicted summary text
			if (input$PI952 == TRUE) {	#If "Show 95% Prediction Intervals" for Regimen 2 is selected, calculate 95% prediction intervals
				lo.G4N2 <- round(CI95lo(sim.data2.tail$S9[sim.data2.tail$ID > 1]))	#2.5th percentile
				hi.G4N2 <- round(CI95hi(sim.data2.tail$S9[sim.data2.tail$ID > 1]))	#97.5th percentile
				outputG4N2Text <- paste("Duration of Severe Neutropenia = ",pop.G4N2," hours (",lo.G4N2," - ",hi.G4N2,")",sep="")	#Population predicted and 95% prediction intervals summary text
			}
			outputG4N2Text	#Print text
		}
		outputG4N2Text.function(Rsim.data2())
	})	#Brackets closing "renderText" expression for "outputG4N2Text"
	#Send summary of duration in Grade 4 neutropenia to ui.R for regimen 3
	output$outputG4N3Text <- renderText({
		outputG4N3Text.function <- function(sim.data3) {
			sim.data3.tail <- ddply(sim.data3, .(ID), lastperID)	#Find the last row for each individual
			pop.G4N3 <- round(sim.data3.tail$S9[sim.data3.tail$ID == 1])	#Population predicted duration in Grade 4 Neutropenia
			outputG4N3Text <- paste("Duration of Severe Neutropenia = ",pop.G4N3," hours",sep="")	#Population predicted summary text
			if (input$PI953 == TRUE) {	#If "Show 95% Prediction Intervals" for Regimen 3 is selected, calculate 95% prediction intervals
				lo.G4N3 <- round(CI95lo(sim.data3.tail$S9[sim.data3.tail$ID > 1]))	#2.5th percentile
				hi.G4N3 <- round(CI95hi(sim.data3.tail$S9[sim.data3.tail$ID > 1]))	#97.5th percentile
				outputG4N3Text <- paste("Duration of Severe Neutropenia = ",pop.G4N3," hours (",lo.G4N3," - ",hi.G4N3,")",sep="")	#Population predicted and 95% prediction intervals summary text
			}
			outputG4N3Text	#Print text
		}
		outputG4N3Text.function(Rsim.data3())
	})	#Brackets closing "renderText" expression for "outputG4N3Text"
	#############
	##_SESSION_##
	#############
	session$onSessionEnded(function() {
		stopApp()	#Stops the R session when Google Chrome/R Studio window closes
	})	#Brackets closing "session"
})	#Brackets closing "shinyServer"
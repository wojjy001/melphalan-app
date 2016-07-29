# server.R script for MelphalanApp
# Reactive objects (i.e., those dependent on widget input) are written here
# ------------------------------------------------------------------------------
# Define the "server" part of the Shiny application
shinyServer(function(input,output,session) {
	###########
	##_INPUT_##
	###########
	# Create an input data frame that stores input patient characteristics
	# Will be used for the three different simulation scenarios
		Rinput.data <- reactive({
			# Call in user-defined widget values
				AGE <- input$AGE	# Numeric input for patient's age
				TBW <- input$TBW	# Numeric input for patient's total body weight
				HT <- input$HT	# Numeric input for patient's height
				SECR <- input$SECR	# Numeric input for patient's serum creatinine
				HCT <- input$HCT	# Numeric input for patient's haematocrit
				ANCBASE <- input$ANCBASE	# Numeric input for patient's baseline absolute neutrophil count
				if (input$SEX == 1) SEX <- 0	# Select input for patient's gender, female = 0
				if (input$SEX == 2) SEX <- 1	# Select input for patient's gender, male = 1
				if (input$RACE == 1 | input$RACE == 3) RACE <- 0	# Select input for patient's race, Caucasian and Unknown = 0
				if (input$RACE == 2) RACE <- 1	# Select input for patient's race, African-American = 1
				if (input$SLC7A5 == 1) SLC7A5 <- 0	# Select input for patient's SLC7A5 genotype, AA or AG = 0
				if (input$SLC7A5 == 2) SLC7A5 <- 1	# Select input for patient's SLC7A5 genotype, GG = 1

			# Calculate secondary parameters based on input
				# Body mass index (BMI)
					BMI <- TBW/(HT/100)^2	# Used to calculate fat free mass
				# Body surface area (BSA)
					BSA <- sqrt((HT*TBW)/3600)	# Use to calculate amount to be administered
				# Creatinine clearance (CRCL) and fat free mass (FFM) based on gender
					if (SEX == 0) {	# Females
						CRCL <- (((140-AGE)*TBW)/(SECR*0.815))*0.85
						FFM <- 9270*TBW/(8780+(244*BMI))
					} else {	# Males
						CRCL <- ((140-AGE)*TBW)/(SECR*0.815)
						FFM <- 9270*TBW/(6680+(216*BMI))
					}

			# Set up input.data
			# Only columns missing values will be the amount to be administered and when G-CSF rescue was administered as each "sim.data" data frame will have a different value for "amt" and "G-CSF"
				input.data <- expand.ev(
					ID = 1:(n+1), # n individuals (plus an additional because the first ID is PRED)
					time = 0,	# time that melphalan dose will be administered
					amt = NA,	# amt in mg/m^2, currently amount per m^2 is unknown
					evid = 1,	# dosing event
					cmt = 1, # dose into compartment 1, i.e., CENT
					rate = -2,	# infusion duration is specified in the model file
					BSA = BSA,	# Required to be stored to calculate "amt"
					BMI = BMI,	# Required to be stored for ui
					FFM = FFM,	# Fat free mass
					CRCL = CRCL,	# Creatinine clearance
					HCT = HCT,	# Haematocrit
					ANCBASE = ANCBASE,	# Baseline absolute neutrophil count
					SEX = SEX,	# Gender
					RACE = RACE,	# Race
					SLC7A5 = SLC7A5,	# Genotype
					GCSF = 0	# Time of administration of G-CSF (will be different for different sim.data), default is on Day 1 (0)
				)
		})	# Brackets closing "Rinput.data"

	###########
	##_DOSE1_##
	###########

	# Simulate a population based on input characteristics
	# Will have it's own specific dose and time of G-CSF administration
		Rsim.data1 <- reactive({
			# Read in reactive input.data
				input.data <- Rinput.data()
			# Read in simulation specific value for G-CSF
				if (input$GCSF1 == 2) input.data$GCSF <- 1	# Select input for when to administer G-CSF (Neupogen), Day 7 = 1
			# Calculate amt to be administered based on patient's BSA and DOSE1
				input.data$amt <- input.data$BSA*input$DOSE1
		  # Simulate
				sim.data1 <- mod %>% data_set(input.data) %>% mrgsim(add = time)
				sim.data1 <- as.data.frame(sim.data1)	#Convert to a data frame so that it is more useful for me!
		})	# Brackets closing "Rsim.data1"

	# Create a data frame that only contains the "PRED" data
		Rpred.data1 <- reactive({
			# Read in reactive expressions
				sim.data1 <- Rsim.data1()
			# Subset out only ID == 1 (PRED individual)
				pred.data1 <- sim.data1[sim.data1$ID == 1,]
		})	# Brackets closing "Rpred.data1"

	# Summarise simulated data as prediction intervals when option is selected
		Rsummary.data1 <- reactive({
			# Read in reactive expressions
				sim.data1 <- Rsim.data1()
			# Summarise data
				sim.data1 <- sim.data1[sim.data1$ID != 1,]	# Do not include ID == 1 - they are PRED
				summary.data1 <- ddply(sim.data1, .(time), summary.function)
		})	# Brackets closing "Rsummary.data1"

	###########
	##_DOSE2_##
	###########

	# Simulate a population based on input characteristics
	# Will have it's own specific dose and time of G-CSF administration
		Rsim.data2 <- reactive({
			# Read in reactive input.data
				input.data <- Rinput.data()
			# Read in simulation specific value for G-CSF
				if (input$GCSF2 == 2) input.data$GCSF <- 1	# Select input for when to administer G-CSF (Neupogen), Day 7 = 1
			# Calculate amt to be administered based on patient's BSA and DOSE1
				input.data$amt <- input.data$BSA*input$DOSE2
		  # Simulate
				sim.data2 <- mod %>% data_set(input.data) %>% mrgsim(add = time)
				sim.data2 <- as.data.frame(sim.data2)	#Convert to a data frame so that it is more useful for me!
		})	# Brackets closing "Rsim.data2"

	# Create a data frame that only contains the "PRED" data
		Rpred.data2 <- reactive({
			# Read in reactive expressions
				sim.data2 <- Rsim.data2()
			# Subset out only ID == 1 (PRED individual)
				pred.data2 <- sim.data2[sim.data2$ID == 1,]
		})	# Brackets closing "Rpred.data2"

	# Summarise simulated data as prediction intervals when option is selected
		Rsummary.data2 <- reactive({
			# Read in reactive expressions
				sim.data2 <- Rsim.data2()
			# Summarise data
				sim.data2 <- sim.data2[sim.data2$ID != 1,]	# Do not include ID == 1 - they are PRED
				summary.data2 <- ddply(sim.data2, .(time), summary.function)
		})	# Brackets closing "Rsummary.data2"

	############
	##_OUTPUT_##
	############
	output$BMI.text <- renderUI({
		input.data <- Rinput.data()
		withMathJax(helpText("Body mass index = ",round(input.data$BMI[1],digits = 1)," \\(kg/m^2\\)"))
	})	# Brackets closing "renderUI" expression

	output$BSA.text <- renderUI({
		input.data <- Rinput.data()
		withMathJax(helpText("Body surface area = ",round(input.data$BSA[1],digits = 1)," \\(m^2\\)"))
	})	# Brackets closing "renderUI" expression

	output$FFM.text <- renderUI({
		input.data <- Rinput.data()
		helpText("Fat free mass = ",round(input.data$FFM[1],digits = 1)," \\(kg\\)")
	})	# Brackets closing "renderUI" expression

	output$CRCL.text <- renderUI({
		input.data <- Rinput.data()
		helpText("Creatinine clearance = ",round(input.data$CRCL[1],digits = 1)," \\(mL/min\\)")
	})	# Brackets closing "renderUI" expression

	# Simulation results for ANC
		output$anc.plot <- renderPlot({
			# Read in reactive data
				pred.data1 <- Rpred.data1()
				pred.data2 <- Rpred.data2()
				summary.data1 <- Rsummary.data1()
				summary.data2 <- Rsummary.data2()

			# Plot ANC over time
				plotobj1 <- NULL
				plotobj1 <- ggplot()
			# Population predicted
				plotobj1 <- plotobj1 + geom_line(aes(x = time,y = ANC),data = pred.data1,colour = "red",size = 1)	# DOSE1
				if (input$NREG > 1) plotobj1 <- plotobj1 + geom_line(aes(x = time,y = ANC),data = pred.data2,colour = "blue",size = 1)	# DOSE2
			# 95% prediction intervals
				if (input$PI1 == TRUE) plotobj1 <- plotobj1 + geom_ribbon(aes(x = time,ymin = CIlo_ANC,ymax = CIhi_ANC),data = summary.data1,fill = "red",alpha = 0.3)	# DOSE1
				if (input$PI2 == TRUE & input$NREG > 1) plotobj1 <- plotobj1 + geom_ribbon(aes(x = time,ymin = CIlo_ANC,ymax = CIhi_ANC),data = summary.data2,fill = "blue",alpha = 0.3)	# DOSE2
			# Grade 4 neutropenia
				plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0.5),linetype = "dashed")
				plotobj1 <- plotobj1 + annotate("text",x = 648,y = 0.8,label = "Grade 4 Neutropenia",size = 5)
			# Axes
				plotobj1 <- plotobj1 + scale_x_continuous("\nTime since melphalan dose (days)",breaks = seq(from = 0,to = max(time.PD),by = 100))
				plotobj1 <- plotobj1 + scale_y_log10("Absolute Neutrophil Count (K/µL)\n",breaks = log.plot.breaks,labels = log.plot.breaks,lim = c(0.0001,NA))
			# Return plot
				print(plotobj1)
		})	# Brackets closing "renderPlot"

	# Simulation results for melphalan concentrations
		output$melph.plot <- renderPlot({
			# Read in reactive data
			pred.data1 <- Rpred.data1()
			pred.data2 <- Rpred.data2()
			summary.data1 <- Rsummary.data1()
			summary.data2 <- Rsummary.data2()

			# Plot ANC over time
				plotobj2 <- NULL
				plotobj2 <- ggplot()
			# Population predicted
				plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE),data = pred.data1,colour = "red",size = 1)	# DOSE1
				if (input$NREG > 1) plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE),data = pred.data2,colour = "blue",size = 1)	# DOSE2
			# 95% prediction intervals
				if (input$PI1 == TRUE) plotobj2 <- plotobj2 + geom_ribbon(aes(x = time,ymin = CIlo_IPRE,ymax = CIhi_IPRE),data = summary.data1,fill = "red",alpha = 0.3)	# DOSE1
				if (input$PI2 == TRUE & input$NREG > 1) plotobj2 <- plotobj2 + geom_ribbon(aes(x = time,ymin = CIlo_IPRE,ymax = CIhi_IPRE),data = summary.data2,fill = "blue",alpha = 0.3)	# DOSE2
			# Axes
				plotobj2 <- plotobj2 + scale_x_continuous("\nTime since melphalan dose (hours)",lim = c(0,12))
				plotobj2 <- plotobj2 + scale_y_log10("Melphalan Concentration (mg/L)\n",lim = c(0.001,NA),breaks = log.plot.breaks,labels = log.plot.breaks)
			# Return plot
				print(plotobj2)
		})	# Brackets closing "renderPlot"

	# Summary of Time spent in Grade 4 Neutropenia for DOSE1
		output$G4N1.text.DOSE1 <- renderText({
			# Read in reactive data
				pred.data1 <- Rpred.data1()
				summary.data1 <- Rsummary.data1()
			# Create a text object
				pred.G4N1 <- round(tail(pred.data1$G4N1,1))	# PRED
				G4N1.text.DOSE1 <- paste0("Duration on Severe Neutropenia = ",pred.G4N1," hours")
				if (input$PI1 == TRUE) {
					CIlo.G4N1 <- round(tail(summary.data1$CIlo_G4N1,1))	# 2.5th percentile
					CIhi.G4N1 <- round(tail(summary.data1$CIhi_G4N1,1))	# 97.5th percentile
					G4N1.text.DOSE1 <- paste0("Duration on Severe Neutropenia = ",pred.G4N1," hours (",CIlo.G4N1," - ",CIhi.G4N1,")")
				}
				G4N1.text.DOSE1
		})	# Brackets closing "renderText"

	# Summary of Time spent in Grade 4 Neutropenia for DOSE2
		output$G4N1.text.DOSE2 <- renderText({
			# Read in reactive data
				pred.data2 <- Rpred.data2()
				summary.data2 <- Rsummary.data2()
			# Create a text object
				pred.G4N1 <- round(tail(pred.data2$G4N1,1))	# PRED
				G4N1.text.DOSE2 <- paste0("Duration on Severe Neutropenia = ",pred.G4N1," hours")
				if (input$PI2 == TRUE) {
					CIlo.G4N1 <- round(tail(summary.data2$CIlo_G4N1,1))	# 2.5th percentile
					CIhi.G4N1 <- round(tail(summary.data2$CIhi_G4N1,1))	# 97.5th percentile
					G4N1.text.DOSE2 <- paste0("Duration on Severe Neutropenia = ",pred.G4N1," hours (",CIlo.G4N1," - ",CIhi.G4N1,")")
				}
				G4N1.text.DOSE2
		})	# Brackets closing "renderText"

  #############
  ##_SESSION_##
  #############
  # Close the R session when Chrome closes
  session$onSessionEnded(function() {
    stopApp()
  })
})  # Brackets closing "shinyServer" function

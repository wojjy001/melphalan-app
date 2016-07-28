#Calculate the time above MIC
#Uses rFIT.data for calculating time above MIC for fitted concentrations
#Uses rSIM.data for calculating time above MIC for simulated concentrations
rMIC.data <- reactive({
	mic.function <- function(FIT.data,SIM.data) {
		if (input$SIMNEW == FALSE) {
			EMIC.data <- subset(FIT.data, ICONC > 20)
			time.EMIC <- max(EMIC.data$TIME) - min(EMIC.data$TIME)
			time.SMIC <- NA
			time.MIC <- data.frame(EMIC = time.EMIC, SMIC = time.SMIC)
		}
		
		if (input$SIMNEW == TRUE) {
			EMIC.data <- subset(FIT.data, ICONC > 20)
			time.EMIC <- max(EMIC.data$TIME) - min(EMIC.data$TIME)			
			SMIC.data <- subset(SIM.data, SCONC > 20)
			time.SMIC <- max(SMIC.data$TIME) - min(SMIC.data$TIME)
			time.MIC <- data.frame(EMIC = time.EMIC, SMIC = time.SMIC)			
		}
		time.MIC
	}	
	MIC.data <- mic.function(rFIT.data(),rSIM.data())
})

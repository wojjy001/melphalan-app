#Define UI (user-interface) for MelphPrototypeApp1
#Sends user-defined input to server.R
#Calls created output from server.R
#-----------------------------------------------------------------------------------
fixedPage(
	fixedRow(
		h2("Melphalan", align = "center")	#Application Title
	),	#Brackets closing "fixedRow"
	hr(),	#Horizontal line - divider
	sidebarLayout(
		sidebarPanel(	#sidebarPanel containing user-defined widgets
			h4("Patient Information"),	#Header for the following widgets
			numericInput("AGE","Age (years):",min = 0,max = 100,value = 70, step = 1),	#Numeric input for patient's age
			numericInput("TBW","Total Body Weight (kg):",min = 0,max = 150,value = 70,step = 0.1),	#Numeric input for patient's total body weight
			numericInput("HT","Height (cm):",min = 0,max = 210,value = 170,step = 1),	#Numeric input for patient's height
			uiOutput("outputBMIText"),	#Patient's calculated body mass index (kg/m^2)
			uiOutput("outputBSAText"),	#Patient's calculated body surface area (m^2)
			uiOutput("outputFFMText"),	#Patient's calculated fat free mass (kg)
			br(),
			numericInput("SECR","Serum creatinine (µmol/L):",min = 0,max = 300,value = 70,step = 1),	#Numeric input for patient's serum creatinine
			uiOutput("outputCRCLText"),	#Patient's calculated creatinine clearance (mL/min)
			br(),
			numericInput("HCT","Haematocrit (%):",min = 0,max = 100,value = 32.5,step = 0.1),	#Numeric input for patient's haematocrit
			numericInput("BUN","Blood Urea Nitrogen (mg/dL):",min = 0,max = 100,value = 14,step = 0.1),	#Numeric input for patient's blood urea nitrogen
			numericInput("WBC","White Blood Cell Count (K/µL):",min = 0,max = 20,value = 4.9,step = 0.1), #Numeric input for patient's white blood cell count
			numericInput("ANC","Baseline Absolute Neutrophil Count (K/µL):",min = 0,max = 100,value = 3.5,step = 0.1), #Numeric input for patient's absolute neutrophil count
			selectInput("SEX","Gender:",choices = list("Female" = 1,"Male" = 2),selected = 1),	#Select input for patient's gender
			selectInput("SLC7A5","SLC7A5 Genotype:",choices = list("AA or AG" = 1,"GG" = 2),selected = 1),	#Select input for patient's SLC7A5
			withMathJax(
				numericInput("LNP53FOLD","Log-transformed p53 mRNA expression level change in PBMCs after treating with 75 µg/mL melphalan \\(ex \\, vivo\\) versus baseline level:",min = 0,max = 100,value = 2.62,step = 0.01)	#Numeric input for patient's LNP53FOLD
			)	#Brackets closing "withMathJax"
		),	#Brackets closing "sidebarPanel"
		mainPanel(	#mainPanel containing output (plots)
			tabsetPanel(
				tabPanel("ANC Profile",
					plotOutput("outputPDPlot"),
					checkboxInput("TUNITS","Change time-axis units from days to hours",value = FALSE)
				),	#Brackets closing "tabPanel"
				tabPanel("Duration of Severe Neutropenia",
					plotOutput("outputG4NPlot")
				),	#Brackets closing "tabPanel"
				tabPanel("Melphalan Profile",
					plotOutput("outputPKPlot")
				)	#Brackets closing "tabPanel"
			),	#Brackets closing "tabsetPanel"
			hr(),
			h4("Dosing Information"),	#Header for the following widgets
			selectInput("NREG","Select number of dosing regimens to compare:",choices = list("1" = 1,"2" = 2,"3" = 3),selected = 1,width = 325),	#Select input for number of dosing regimens to compare
			fixedRow(
				column(width = 6,
					HTML('<h5 style="color:#FF0000">Regimen 1:</h5>')
				),	#Brackets closing "column"
				column(width = 6,
					checkboxInput("PI951",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE)
				)	#Brackets closing "column"
			),	#Brackets closing "fixedRow"
			fixedRow(
				column(width = 6,
					withMathJax(
						sliderInput("DOSE1","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10)	#Slider input for Melphalan dose
					)	#Brackets closing "withMathJax"
				),	#Brackets closing "column"
				column(width = 6,
					selectInput("GCSF1","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1)	#Select input for when to administer G-CSF (Neupogen)	
				)	#Brackets closing "column"
			),	#Brackets closing "fixedRow"
			fixedRow(
				column(width = 6, 
					textOutput("outputG4N1Text")
				)	#Brackets closing "column"
			),	#Brackets closing "fixedRow"
			hr(),
			conditionalPanel(condition = "input.NREG > 1",
				fixedRow(
					column(width = 6,
						HTML('<h5 style="color:#0000FF">Regimen 2:</h5>')
					),	#Brackets closing "column"
					column(width = 6,
						checkboxInput("PI952",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE)	
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"		
				fixedRow(
					column(width = 6,
						withMathJax(
							sliderInput("DOSE2","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 200,step = 10)	#Slider input for Melphalan dose
						)	#Brackets closing "withMathJax"
					),	#Brackets closing "column"
					column(width = 6,
						selectInput("GCSF2","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1)	#Select input for when to administer G-CSF (Neupogen)	
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"
				fixedRow(
					column(width = 6,
						textOutput("outputG4N2Text")
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"
				hr()
			),	#Brackets closing "conditionalPanel"
			conditionalPanel(condition = "input.NREG > 2",
				fixedRow(
					column(width = 6,
						HTML('<h5 style="color:#006600">Regimen 3:</h5>')
					),	#Brackets closing "column"
					column(width = 6,
						checkboxInput("PI953",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE)	
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"		
				fixedRow(
					column(width = 6,
					withMathJax(
							sliderInput("DOSE3","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 300,step = 10)	#Slider input for Melphalan dose
						)	#Brackets closing "withMathJax"
					),	#Brackets closing "column"
					column(width = 6,
						selectInput("GCSF3","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1)	#Select input for when to administer G-CSF (Neupogen)	
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"
				fixedRow(
					column(width = 6,
						textOutput("outputG4N3Text")
					)	#Brackets closing "column"
				),	#Brackets closing "fixedRow"
				hr()
			)	#Brackets closing "conditionalPanel"
		)	#Brackets closing "mainPanel"
	)	#Brackets closing "sidebarLayout"
)	#Brackets closing "fixedPage"

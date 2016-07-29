# Define UI (user-interface) for MelphalanApp
# Sends user-defined input to server.R
# Calls created output from server.R
# ------------------------------------------------------------------------------
fixedPage(
	fixedRow(
		h2("Melphalan", align = "center")	# Application Title
	),	# Brackets closing "fixedRow"
	hr(),	# Horizontal line - divider
	sidebarLayout(
		sidebarPanel(	# sidebarPanel containing user-defined widgets
			h4("Patient Information"),	# Header for the following widgets
			numericInput("AGE","Age (years):",min = 0,max = 100,value = 70, step = 1),	# Numeric input for patient's age
			numericInput("TBW","Total Body Weight (kg):",min = 0,max = 150,value = 70,step = 0.1),	# Numeric input for patient's total body weight
			numericInput("HT","Height (cm):",min = 0,max = 210,value = 170,step = 1),	# Numeric input for patient's height
			uiOutput("BMI.text"),	# Patient's calculated body mass index (kg/m^2)
			uiOutput("BSA.text"),	# Patient's calculated body surface area (m^2)
			uiOutput("FFM.text"),	# Patient's calculated fat free mass (kg)
			br(),
			numericInput("SECR","Serum creatinine (µmol/L):",min = 0,max = 300,value = 70,step = 1),	# Numeric input for patient's serum creatinine
			uiOutput("CRCL.text"),	# Patient's calculated creatinine clearance (mL/min)
			br(),
			numericInput("HCT","Haematocrit (%):",min = 0,max = 100,value = 32.5,step = 0.1),	# Numeric input for patient's haematocrit
			numericInput("ANCBASE","Baseline Absolute Neutrophil Count (K/µL):",min = 0,max = 100,value = 3.5,step = 0.1), # Numeric input for patient's absolute neutrophil count
			selectInput("SEX","Gender:",choices = list("Female" = 1,"Male" = 2),selected = 1),	# Select input for patient's gender
			selectInput("RACE","Race/Ethnicity:",choices = list("Caucasian" = 1,"African-American" = 2,"Unknown" = 3),selected = 1),	# Select input for patient's race
			selectInput("SLC7A5","SLC7A5 Genotype:",choices = list("AA or AG" = 1,"GG" = 2),selected = 1)	# Select input for patient's SLC7A5
		),	# Brackets closing "sidebarPanel"
		mainPanel(	# mainPanel containing output (plots)
			tabsetPanel(
				tabPanel("ANC Profile",
					plotOutput("anc.plot")
				),	# Brackets closing "tabPanel"
				tabPanel("Duration of Severe Neutropenia"
				),	# Brackets closing "tabPanel"
				tabPanel("Melphalan Profile",
					plotOutput("melph.plot")
				)	# Brackets closing "tabPanel"
			),	# Brackets closing "tabsetPanel"
			hr(),
			h4("Dosing Information"),	# Header for the following widgets
			selectInput("NREG","Select number of dosing regimens to compare:",choices = list("1" = 1,"2" = 2,"3" = 3),selected = 1,width = 325),	# Select input for number of dosing regimens to compare
			fixedRow(
				column(6,
					HTML('<h5 style="color:#FF0000">Regimen 1:</h5>')
				),	# Brackets closing "column"
				column(6,
					checkboxInput("PI1",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE)
				)	# Brackets closing "column"
			),	# Brackets closing "fixedRow"
			fixedRow(
				column(6,
					withMathJax(
						sliderInput("DOSE1","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10)	# Slider input for Melphalan dose
					)	# Brackets closing "withMathJax"
				),	# Brackets closing "column"
				column(6,
					selectInput("GCSF1","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1)	# Select input for when to administer G-CSF (Neupogen)
				)	# Brackets closing "column"
			),	# Brackets closing "fixedRow"
			fixedRow(
				column(6,
					textOutput("G4N1.text.DOSE1")
				)	# Brackets closing "column"
			),	# Brackets closing "fixedRow"
			hr(),
			conditionalPanel(condition = "input.NREG > 1",
				fixedRow(
					column(6,
						HTML('<h5 style="color:#0000FF">Regimen 2:</h5>')
					),	# Brackets closing "column"
					column(6,
						checkboxInput("PI2",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE)
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(6,
						withMathJax(
							sliderInput("DOSE2","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10)	# Slider input for Melphalan dose
						)	# Brackets closing "withMathJax"
					),	# Brackets closing "column"
					column(6,
						selectInput("GCSF2","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1)	# Select input for when to administer G-CSF (Neupogen)
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(6,
						textOutput("G4N1.text.DOSE2")
					)	# Brackets closing "column"
				)	# Brackets closing "fixedRow"
			)	# Brackets closing "conditionalPanel"
		)	# Brackets closing "mainPanel"
	)	# Brackets closing "sidebarLayout"
)	# Brackets closing "fixedPage"

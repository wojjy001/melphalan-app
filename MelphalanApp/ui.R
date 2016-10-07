# Define UI (user-interface) for MelphalanApp
# Sends user-defined input to server.R
# Calls created output from server.R
# Now using shinydashboard for the user-interface
# ------------------------------------------------------------------------------
# Application's header
header <-
	dashboardHeader(
		title = "Melphalan",
		titleWidth = 400
	)	# Brackets closing "dashboardHeader"
# Application's sidebar
sidebar <-
	dashboardSidebar(
		width = 400,	# Width of sidebar the same as width of header
		sidebarMenu(id = "sidebarmenu",
			menuItem("About",tabName = "about",icon = icon("question-circle"),
				menuSubItem("Background and Objective",tabName = "objective",icon = icon("angle-double-right")),
				menuSubItem("ASCPT 2016 Abstract",tabName = "abstract",icon = icon("angle-double-right")),
				menuSubItem("Population PK/PD Model",tabName = "model",icon = icon("angle-double-right")),
				menuSubItem("Resources",tabName = "packages",icon = icon("angle-double-right")),
				menuSubItem("Acknowledgements",tabName = "acknowledgements",icon = icon("angle-double-right"))
			),
			menuItem("Application",tabName = "application",icon = icon("line-chart")),
			conditionalPanel(condition = "input.sidebarmenu === 'application'",
				fixedRow(
					column(11,offset = 1,
						numericInput("AGE","Age (years):",min = 0,max = 100,value = 70, step = 1),	# Numeric input for patient's age
						numericInput("TBW","Total Body Weight (kg):",min = 0,max = 150,value = 70,step = 0.1),	# Numeric input for patient's total body weight
						numericInput("HT","Height (cm):",min = 0,max = 210,value = 170,step = 1)	# Numeric input for patient's height
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(11,offset = 2,
						em(uiOutput("BMI.text")),	# Patient's calculated body mass index (kg/m^2)
						em(uiOutput("BSA.text")),	# Patient's calculated body surface area (m^2)
						em(uiOutput("FFM.text"))	# Patient's calculated fat free mass (kg)
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(11,offset = 1,
						numericInput("SECR","Serum creatinine (mg/dL):",min = 0,max = 10,value = 1,step = 0.01)	# Numeric input for patient's serum creatinine
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(11,offset = 2,
						em(uiOutput("CRCL.text"))	# Patient's calculated creatinine clearance (mL/min)
					)	# Brackets closing "column"
				),	# Brackets closing "fixedRow"
				fixedRow(
					column(11,offset = 1,
						numericInput("HCT","Haematocrit (%):",min = 0,max = 100,value = 32.5,step = 0.1),	# Numeric input for patient's haematocrit
						numericInput("ANCBASE","Baseline Absolute Neutrophil Count (K/µL):",min = 0,max = 100,value = 3.5,step = 0.1), # Numeric input for patient's absolute neutrophil count
						selectInput("SEX","Gender:",choices = list("Female" = 1,"Male" = 2),selected = 1),	# Select input for patient's gender
						selectInput("RACE","Race/Ethnicity:",choices = list("Caucasian" = 1,"African-American" = 2,"Unknown" = 3),selected = 1),	# Select input for patient's race
						selectInput("SLC7A5","SLC7A5 Genotype:",choices = list("AA or AG" = 1,"GG" = 2),selected = 1)	# Select input for patient's SLC7A5
					)	# Brackets closing "column"
				)	# Brackets closing "fixedRow"
			)	# Brackets closing "conditonalPanel"
		)	# Brackets closing "sidebarMenu"
	)	# Brackets closing "dashboardSidebar"
# Application's body
body <-
	dashboardBody(
		tabItems(
			tabItem(tabName = "abstract",
				includeMarkdown("ascpt_abstract.Rmd")
			),	# Brackets closing "tabItem" for "abstract"
			tabItem(tabName = "model",
				pre(includeText("model.R"))
			),	# Brackets closing "tabItem" for "model"
			tabItem(tabName = "packages",
				pre(htmlOutput("session.info"))
			),	# Brackets closing "tabItem" for "packages"
			tabItem(tabName = "application",
				tabBox(
					tabPanel("ANC Profile",
						plotOutput("anc.plot",height = 450)
					),	# Brackets closing "tabPanel"
					tabPanel("Duration in Grade 4 Neutropenia",
						plotOutput("g4n.plot",height = 450)
					),	# Brackets closing "tabPanel"
					tabPanel("Melphalan Profile",
						plotOutput("melph.plot",height = 450)
					),	# Brackets closing "tabPanel"
					checkboxInput("PI",paste("Show 95% Prediction Intervals (n = ",n,")",sep=""),value = FALSE),
					title = tagList(icon("line-chart"),"Simulated Output"),
					width = 12
				),	# Brackets closing "tabBox"
				box(
					selectInput("NREG","Select number of dosing regimens to compare:",choices = list("1" = 1,"2" = 2,"3" = 3),selected = 1,width = 325),	# Select input for number of dosing regimens to compare
					hr(),
					fixedRow(
						column(4,
							h5(strong("Regimen 1:"),style = "color:#F8766D"),
							withMathJax(
								sliderInput("DOSE1","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10,width = 325)	# Slider input for Melphalan dose
							),	# Brackets closing "withMathJax"
							selectInput("GCSF1","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1,width = 325),	# Select input for when to administer G-CSF (Neupogen)
							strong(textOutput("G4N1.text.DOSE1"),style = "color:#F8766D")
						),	# Brackets closing "column"
						conditionalPanel(condition = "input.NREG > 1",
							column(4,
								h5(strong("Regimen 2:"),style = "color:#619CFF"),
								withMathJax(
									sliderInput("DOSE2","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10,width = 325)	# Slider input for Melphalan dose
								),	# Brackets closing "withMathJax"
								selectInput("GCSF2","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1,width = 325),	# Select input for when to administer G-CSF (Neupogen)
								strong(textOutput("G4N1.text.DOSE2"),style = "color:#619CFF")
							)	# Brackets closing "column"
						),	# Brackets closing "conditionalPanel"
						conditionalPanel(condition = "input.NREG > 2",
							column(4,
								h5(strong("Regimen 3:"),style = "color:#00BA38"),
								withMathJax(
									sliderInput("DOSE3","Melphalan Dose (\\(mg/m^2\\)):",min = 0,max = 400,value = 100,step = 10,width = 325)	# Slider input for Melphalan dose
								),	# Brackets closing "withMathJax"
								selectInput("GCSF3","Administration of G-CSF (Neupogen):",choices = list("Day 1" = 1,"Day 7" = 2),selected = 1,width = 325),	# Select input for when to administer G-CSF (Neupogen)
								strong(textOutput("G4N1.text.DOSE3"),style = "color:#00BA38")
								# textOutput("G4N1.text.DOSE3")
							)	# Brackets closing "column"
						)	# Brackets closing "conditionalPanel"
					),	# Brackets closing "fixedRow"
					width = 12,
					status = "primary"
				)	# Brackets closing "box"
			)	# Brackets closing "tabItem" for "application"
		)	# Brackets closing "tabItems"
	)	# Brackets closing "dashboardBody"
# ------------------------------------------------------------------------------
# User-interface Object
  dashboardPage(header,sidebar,body,skin = "blue")

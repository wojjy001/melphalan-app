# ui.R script for DoxyApp
# The user-interface and widget input for the Shiny application is defined here
# Sends user-defined input to server.R, calls created output from server.R
# Now using shinydashboard for the user-interface
# ------------------------------------------------------------------------------
# Application's header
header <-
  dashboardHeader(
		title = "Doxycycline",
		titleWidth = 250
	)	# Brackets closing "dashboardHeader"
# Application's sidebar
sidebar <-
	dashboardSidebar(
		width = 250,	# Width of sidebar the same as width of header
		sidebarMenu(
		  menuItem("About",tabName = "about",icon = icon("question-circle"),
        menuSubItem("Objective",tabName = "objective",icon = icon("angle-double-right")),
        menuSubItem("Model (mrgsolve code)",tabName = "model",icon = icon("angle-double-right")),
        menuSubItem("Acknowledgements",tabName = "acknowledgements",icon = icon("angle-double-right"))
      ),  # Brackets closing "menuItem"
      menuItem("Simulations",tabName = "sim",icon = icon("line-chart"))
		)	# Brackets closing "sidebarMenu"
	) # Brackets closing "dashboardSidebar"
# Application's body
body <-
	dashboardBody(
    tags$head(
			tags$link(rel = "stylesheet",type = "text/css",href = "custom.css")
		),
		tabItems(
      tabItem(tabName = "objective",
        includeMarkdown("objective.Rmd")
			), # Brackets closing "tabItem" for "objective"
      tabItem(tabName = "model",
        pre(includeText("model.R"))
      ),  # Brackets closing "tabItem" for "model"
      tabItem(tabName = "acknowledgements",
        includeMarkdown("acknowledgements.Rmd"),
        img(src = "ACP_logo.png",width = 225,height = 75) # University of South Australia, Australian Centre for Pharmacometrics logo
      ),  # Brackets closing "tabItem" for "acknowledgements"
      tabItem(tabName = "sim",
        box(
          fixedRow(
            column(4,
              selectInput("DOSE_REG","Dose Regimen:",choices = list("Single-dose" = 1,"Multiple-dose (standard infection)" = 2,"Multiple-dose (severe infection)"= 3)),
              fixedRow(
                column(3,
                  p(strong("Description:"))
                ),  # Brackets closing "column"
                column(9,
                  conditionalPanel(condition = "input.DOSE_REG == 1",
                    p("120 mg Doryx MPC over 96 hours ",em(strong("or"))),
                    p("100 mg Doryx Tablet over 96 hours")
                  ),  # Brackets closing "conditionalPanel"
                  conditionalPanel(condition = "input.DOSE_REG == 2",
                    p("120 mg Doryx MPC every 12 hours on Day 1, then 120 mg every 24 hours on Days 2 to 7 ",em(strong("or"))),
                    p("100 mg Doryx Tablet every 12 hours on Day 1, then 100 mg every 24 hours on Days 2 to 7")
                  ),  # Brackets closing "conditionalPanel"
                  conditionalPanel(condition = "input.DOSE_REG == 3",
                    p("120 mg Doryx MPC every 12 hours for 7 days ",em(strong("or"))),
                    p("120 mg Doryx Tablet every 12 hours for 7 days")
                  )  # Brackets closing "conditionalPanel"
                ) # Brackets closingn "column"
              ) # Brackets closing "fixedRow"
            ),  # Brackets closing "column"
            column(4,
              selectInput("SIM_STUDY","Simulation Study:",choices = list("Fed versus Fasted" = 1,"Doryx MPC versus Doryx Tablet" = 2,"Male versus Female" = 3)),
              fixedRow(
                column(4,
                  p(strong("Plot Legend:"))
                ),  # Brackets closing "column"
                column(8,
                  conditionalPanel(condition = "input.SIM_STUDY == 1",
                    p(strong("Fasted"),style = "color:#B22222"),
                    p(strong("Fed"),style = "color:#3c8dbc")
                  ),  # Brackets closing "conditionalPanel"
                  conditionalPanel(condition = "input.SIM_STUDY == 2",
                    p(strong("Doryx MPC"),style = "color:#B22222"),
                    p(strong("Doryx Tablet"),style = "color:#3c8dbc")
                  ),  # Brackets closing "conditionalPanel"
                  conditionalPanel(condition = "input.SIM_STUDY == 3",
                    p(strong("Female"),style = "color:#B22222"),
                    p(strong("Male"),style = "color:#3c8dbc")
                  )  # Brackets closing "conditionalPanel"
                ) # Brackets closing "column"
              ) # Brackets closing "fixedRow"
            ),  # Brackets closing "column"
            column(4,
              selectInput("PI","Prediction intervals:",choices = list("No Prediction Intervals" = 1,"90% Prediction Intervals" = 2,"95% Prediction Intervals" = 3)),
              checkboxInput("SUMSTATS","Show summary statistics tables",value = FALSE), # Calculate Tmax, Cmax and AUC. Show prediction intervals if a "type" of prediction intervals is previously selected (as above). Show for each facet if "FACET" is selected above.
              checkboxInput("LOGS","Plot concentrations on a log-scale",value = FALSE)
            ) # Brackets closing "column"
          ),  # Brackets closing "fixedRow"
          title = strong("Simulation Options"),
          solidHeader = TRUE,
          status = "primary",
          width = 12
        ),  # Brackets closing "box"
        box(
          fixedRow(
            column(6,
              conditionalPanel(condition = "input.SIM_STUDY == 1",
                h4(strong("Doryx MPC"))
              ),  # Brackets closing "conditionalPanel"
              conditionalPanel(condition = "input.SIM_STUDY == 2",
                h4(strong("Fasted"))
              ),  # Brackets closing "conditionalPanel"
              conditionalPanel(condition = "input.SIM_STUDY == 3",
                h4(strong("Doryx MPC"))
              ),  # Brackets closing "conditionalPanel"
              plotOutput("Rplot1"),
              conditionalPanel(condition = "input.SUMSTATS",
                hr(),
                uiOutput("Rtable1")
              ) # Brackets closing "conditionalPanel"
            ),  # Brackets closing "column"
            column(6,
              conditionalPanel(condition = "input.SIM_STUDY == 1",
                h4(strong("Doryx Tablet"))
              ),  # Brackets closing "conditionalPanel"
              conditionalPanel(condition = "input.SIM_STUDY == 2",
                h4(strong("Fed"))
              ),  # Brackets closing "conditionalPanel"
              conditionalPanel(condition = "input.SIM_STUDY == 3",
                h4(strong("Doryx Tablet"))
              ),  # Brackets closing "conditionalPanel"
              plotOutput("Rplot2"),
              conditionalPanel(condition = "input.SUMSTATS",
                hr(),
                uiOutput("Rtable2")
              ) # Brackets closing "conditionalPanel"
            ),  # Brackets closing "column"
            align = "center"
          ), # Brackets closing "fixedRow"
          title = strong("Simulated Concentration-Time Profiles"),
          solidHeader = TRUE,
          status = "primary",
          width = 12
        ) # Brackets closing "box"
      )  # Brackets closing "tabItem" for "sim"
		)  # Brackets closing "tabItems"
	) # Brackets closing "dashboardBody"
# ------------------------------------------------------------------------------
# User-interface Object
  dashboardPage(header,sidebar,body,skin = "blue")

# global.R script for MelphalanApp
# Objects that are not reactive are written here
# This also a safe place for functions that are then used in server.R
# ------------------------------------------------------------------------------
# Load package libraries
  library(shiny)
  library(shinydashboard)  # Package for making cooler user-interface for Shiny applications
  library(ggplot2)  # Plotting
  library(grid)  # Plotting
  library(plyr)  # Split and rearrange data, ddply function
  library(dplyr)  # New plyr
  library(mrgsolve) # Metrum differential equation solver for pharmacometrics
  library(xtable) # Formatting tables

# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# Function for making column headings in "bold"
  bold <- function(x) {paste0('{\\textbf{',x,'}}')}

# Source model code
  source("model.R")

# ------------------------------------------------------------------------------
# Define time sequence - using mrgsolve's tgrid function
  # Simulate concentrations and ANC for 30 day time-period (720 hours)
    time.PK <- seq(from = 0,to = 12,by = 0.25) # Intense sampling for PK concentrations
    time.PD <- seq(from = 0,to = 720,by = 12)  # Daily sampling for ANC
    time <- sort(c(unique(c(time.PK,time.PD))))

# Plot breaks
  log.plot.breaks <- c(0.001,0.01,0.1,1,10,100,1000,10000)

# Set number of individuals that make up the 95% prediction intervals
	n <- 2000
# Set seed for reproducible numbers
	set.seed(123456)
# One per ID function
  oneperID <- function(x) tail(x,1)

# Specify percentile probabilities given input
	input.CIlo <- 0.025	# 2.5th percentile
	input.CIhi <- 0.975	# 97.5th percentile
# Summary function for median and prediction intervals
	summary.function <- function(x) {
		# Calculate and combine results
			summary <- c(
				# Melphalan concentrations
				"CIlo_IPRE" = quantile(x$IPRE,probs = input.CIlo,names = FALSE),
				"CIhi_IPRE" = quantile(x$IPRE,probs = input.CIhi,names = FALSE),
				# Absolute neutrophil counts
				"CIlo_ANC" = quantile(x$ANC,probs = input.CIlo,names = FALSE),
				"CIhi_ANC" = quantile(x$ANC,probs = input.CIhi,names = FALSE),
				# Time spent in Grade 4 Neutropenia
				"CIlo_G4N1" = quantile(x$G4N1,probs = input.CIlo,names = FALSE),
				"CIhi_G4N1" = quantile(x$G4N1,probs = input.CIhi,names = FALSE)
			)
	}

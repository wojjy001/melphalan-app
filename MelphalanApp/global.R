# global.R script for DoxyApp
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
  # Simulate concentrations for 96 hour time-period
    tgrid.fine <- tgrid(0,24,0.25) # Intense sampling
    tgrid.coarse <- tgrid(36,96,12)  # Less intense sampling later in the interval
    TIME.tgrid <- c(tgrid.fine,tgrid.coarse)
  # Time sequences for the multiple dosing scenarios (10 days)
    time.multiple <- 0:240  # Intense sampling (every hour)

# Plot breaks
  plot.breaks <- c(10,seq(from = 500,to = 6000,by = 500))
  log.plot.breaks <- c(10,100,1000,10000)

# Set number of individuals that make up the 95% prediction intervals
	n <- 1000
# Set seed for reproducible numbers
	set.seed(123456)
# One per ID function
  oneperID <- function(x) tail(x,1)

# Process NONMEM simulation output and compare with mrgsolve
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()

# Set working directory
	setwd("/Volumes/Prosecutor/PhD/melphalan-app/Simulation")
	master.dir <- "/Volumes/Prosecutor/PhD/melphalan-app/Simulation/"
	dir.extension <- ".nm7"

# Load libraries
	library(ggplot2)
	library(grid)

# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# ---------------------------------------------------------------------------
# Run name
	runname <- "PKPD_Final_raceSIM_JW"

# # Process the fit file
# 	processSIMdata <- function(model.path.file.ext)
# 	{ #begin processSIMdata
#
# 	#Work out the nonmem model dir and model file name
# 	model.path.file <- gsub(".ctl", "", model.path.file.ext)
# 	model.path <- dirname(model.path.file)
# 	model.file.ext <- basename(model.path.file.ext)
# 	model.file <- gsub(".ctl","",model.file.ext)
#
# 	model.dir <- paste(master.dir,"/",model.path, sep="")
#
# 	#Work out the nonmem output dir
# 	 output.dir <- paste(master.dir,"/",model.path.file,dir.extension, sep="")
# 	 setwd(output.dir)
#
# 	#Work out the SIM fit filename
# 	SIM.file.ext <- paste(model.file,".fit",sep="")
# 	#Work out the SIM fit filename
# 	SIM.file.ext.out <- paste(model.file,".fit.csv",sep="")
#
#
# 	#Need to remove lines with "TABLE NO. 1 in them" and the header lines for each new subject
# 	 #Strip the unwanted lines
# 	 indata <- readLines(SIM.file.ext)
# 	 tablelines <- grep("TABLE NO.  1",indata)  #May be installation specific
# 	 headerlines <- grep(" ID",indata)          #May be installation specific
#
# 	 #Extract the information in the header line for column names
# 	 header <- indata[headerlines[1]]
# 	 header <- scan(textConnection(header), what="char")
# 	 colnum <- length(header)
#
# 	 #Strip out the unwanted lines
# 	 badlines <- c(tablelines,headerlines)
# 	 indata <- indata[-badlines]
#
# 	 #replace white space with commas
# 	 for (i in 1:length(indata))
# 		{
# 		indata[i] <- gsub('[[:space:]]+',',',indata[i])
# 		}
#
# 	 #write to a file
# 	 writeLines(indata, "SIMtemp.txt")
#
# 	 #read again as csv file
# 	 SIMdata <- read.csv("SIMtemp.txt", header=F)
# 	 SIMdata <- SIMdata[,-1]     #delete first blank column
# 	 names(SIMdata) <- header
#
# 	 #The NONMEM output does not make a new ID number for each simulation
# 	 #Therefore make a list of SIM numbers
# 	 nsims <- length(tablelines)
# 	 numtimepoints <- length(SIMdata$ID)
# 	 numtimespersubject <- numtimepoints/nsims
# 	 SIMdata$SIM <- rep(1:nsims, each=numtimespersubject)
#
# 	 #Write the SIMdata to a file
# 	 write.csv(SIMdata,SIM.file.ext.out, row.names=F)
#
# 	 #tidy up
# 	 #consider deleting the original fit file?
# 	 file.remove("SIMtemp.txt")
# 	 setwd(master.dir)
#
# 	 #SIMdata
#
# 	} #end processSIMdata
# 	processSIMdata(paste0(runname,".ctl"))

#Read the simulated data
	SIM.data <- read.csv(paste0(runname,".nm7/",runname,".fit.csv"),stringsAsFactors=F)

#Change working directory
	setwd(paste0(master.dir,"/",runname,".nm7"))

# ------------------------------------------------------------------------------
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)

# Remove the first row for each individual
	SIM.data <- SIM.data[SIM.data$TIME > 0 | c(SIM.data$TIME == 0 & SIM.data$IPRE == 0),]
# ANC predictions were accidentally transformed
	SIM.data$IPRED[SIM.data$CMT == 4] <- (SIM.data$IPRED[SIM.data$CMT == 4]*0.2+1)^(1/0.2)
	SIM.data$PRED[SIM.data$CMT == 4] <- (SIM.data$PRED[SIM.data$CMT == 4]*0.2+1)^(1/0.2)

# Plot results - melphalan concentrations
	plotobj1 <- NULL
	plotobj1 <- ggplot(SIM.data[SIM.data$CMT == 1,])
	plotobj1 <- plotobj1 + stat_summary(aes(x = TIME,y = PRED),geom = "line",fun.y = median,colour = "red")
	plotobj1 <- plotobj1 + stat_summary(aes(x = TIME,y = IPRED),geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)",lim = c(0,10))
	plotobj1 <- plotobj1 + scale_y_log10("Melphalan Concentration (mg/L)\n",breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
	print(plotobj1)

# Plot results - ANC
	plotobj2 <- NULL
	plotobj2 <- ggplot(SIM.data[SIM.data$CMT == 4,])
	plotobj2 <- plotobj2 + stat_summary(aes(x = TIME,y = PRED),geom = "line",fun.y = median,colour = "red")
	plotobj2 <- plotobj2 + stat_summary(aes(x = TIME,y = IPRED),geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.5),linetype = "dashed")
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
	plotobj2 <- plotobj2 + scale_y_log10("Melphalan Concentration (mg/L)\n",breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
	print(plotobj2)

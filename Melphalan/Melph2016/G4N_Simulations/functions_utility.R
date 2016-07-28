#file paths------------------------------------------------------------------------------

#Master data directory
  master.dir <- "D:/Wojciechowski/Melphalan/G4N_Simulations/"		#Uni Computer


#Extension for model directories made by WfN

  dir.extension <- ".nm7"
 

#Opposite of %in% functions--------------------------------------------------------------- 
 "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y from match help
 
 
#ggplot2 specific functions---------------------------------------------------------------
count.unique <- function(x) length(unique(x))

#Function to count numbers in a boxplot
boxplot.give.n <- function(x)
{
   return(c(y = median(x), label = length(x)))
}



to.png <- function(plotobj.in,maintext.in)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=900, height=700, pointsize=20)
  print(plotobj.in)
  
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8))
		 
  dev.off()
 }
 
to.png.sqr <- function(plotobj.in,maintext.in)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=650, height=500, pointsize=14)
   print(plotobj.in)
   
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
  dev.off()
 }    


to.png.tiny <- function(plotobj.in,maintext.in)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=250, height=250, pointsize=12)
   print(plotobj.in)
   
  # #Stamp the plot
  # time.date <- Sys.time()
  # working.dir <- getwd()
  # grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            # x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            # just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
   dev.off()
 }      

 
 to.png.wide <- function(plotobj.in,maintext.in)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=960, height=500, pointsize=14)
   print(plotobj.in)
   
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
  dev.off()
 }      

 
to.png.long <- function(plotobj.in,maintext.in)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=900, height=1200, pointsize=14)
  print(plotobj.in)
  
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8))
		 
  dev.off()
 } 
 
 


 
#--------------------------------------------------------------- 
#Data summary functions

#Define a function for length without NA's
   lengthNA <- function(x) length(na.omit(x))
   
#90% confidence interval functions
 CI90lo <- function(x) quantile(x, probs=0.05)
 CI90hi <- function(x) quantile(x, probs=0.95)
     
   
#Define a function for geometric mean
  geomean <- function(x, na.rm=F)
  {  
  if (na.rm==T) x <- x[is.na(x)==F]
  exp(mean(log(x)))
  }
  #Note x cannot be negative or zero 

#Median and 90% tolerance intervals
  sumfunc90 <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat3 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    stat4 <-  lengthNA(x)
    result <- c("median"=stat1, "low90"=stat2, "hi90"=stat3, "n"=stat4)
    result
  }
  
 #Median and 95% tolerance intervals
  sumfunc95 <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.025, na.rm=T, names=F)  #95%CI
    stat3 <-  quantile(x, probs=0.975, na.rm=T, names=F)
    stat4 <-  lengthNA(x)
    result <- c("median"=stat1, "low95"=stat2, "hi95"=stat3, "n"=stat4)
    result
  }
  
  
#Mean, sd and CV
  sumfuncCV <- function(x)
  {
    stat1 <-  mean(x, na.rm=T)
    stat2 <-  sd(x, na.rm=T)  
    stat3 <-  stat2/stat1*100
    stat4 <-  lengthNA(x)
    result <- c("mean"=stat1, "sd"=stat2, "cv"=stat3, "n"=stat4)
    result
  }
  
  
#Geomean, mean, sd, CV, min, max, n - for Millenium MLN 8237
  sumfuncMLN <- function(x)
  {
    stat1 <-  geomean(x, na.rm=T)
	stat2 <-  mean(x, na.rm=T)
    stat3 <-  sd(x, na.rm=T)  
    stat4 <-  stat3/stat2*100
	stat5 <-  min(x, na.rm=T)  
	stat6 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat7 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    stat8 <-  max(x, na.rm=T) 
	stat9 <-  lengthNA(x)
    result <- c("gmean"=stat1, "mean"=stat2, "sd"=stat3, "cv"=stat4, "min"=stat5, "lo90"=stat6, "hi90"=stat7, "max"=stat8, "n"=stat9)
    result
  }
  
  
 #Summarize distribution by percentiles
  sumfuncPercentile <- function(x)
  {
    stat1 <-  quantile(x, probs=0.05, na.rm=T, names=F) 
	stat2 <-  quantile(x, probs=0.10, na.rm=T, names=F) 
    stat3 <-  quantile(x, probs=0.25, na.rm=T, names=F) 
    stat4 <-  quantile(x, probs=0.50, na.rm=T, names=F) 
	stat5 <-  quantile(x, probs=0.75, na.rm=T, names=F)  
	stat6 <-  quantile(x, probs=0.90, na.rm=T, names=F)  
    stat7 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    result <- c("05perct"=stat1, "10perct"=stat2, "25perct"=stat3, "50perct"=stat4, "75perct"=stat5, "90perct"=stat6, "95perct"=stat7)
    result
  } 

  
#Mean, sd, min and max & n
  sumfuncRange <- function(x)
  {
    stat1 <-  mean(x, na.rm=T)
	stat2 <-  sd(x, na.rm=T)  
    stat3 <-  min(x, na.rm=T)  
    stat4 <-  max(x, na.rm=T)  
    stat5 <-  lengthNA(x)
    result <- c("mean"=stat1,"sd"=stat2, "min"=stat3, "max"=stat4, "n"=stat5)
    result
  }
    
  

#Median etc for boxplot
  sumfuncBOX <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.025, na.rm=T, names=F) 
	stat3 <-  quantile(x, probs=0.25, na.rm=T, names=F)
	stat4 <-  quantile(x, probs=0.75, na.rm=T, names=F)
    stat5 <-  quantile(x, probs=0.975, na.rm=T, names=F)
    result <- c("median"=stat1, "q025"=stat2, "q25"=stat3, "q75"=stat4, "q975"=stat5)
    result
  }
    
 
#---------------------------------------------------------------------------------------------------------------------
#Functions for trimming extreme values  

trimtoNAci90 <- function(x)
{
    #set values outside 5th and 95th percentile to NA
    limits <- quantile(x, prob=c(.05, .95),na.rm=T)
    x[x <= limits[1]] <- NA
	x[x >= limits[2]] <- NA
	x
}	

#test1 <- rnorm(1000, mean=100, sd=20)
#test2 <- trimtoNA1(test1)
#hist(test1)
#hist(test2)


trimtoNAz3 <- function(x)
{
    #set values outside zscore=3 to NA
    limits <- quantile(x, prob=c(.0027, .9973),na.rm=T)
    x[x <= limits[1]] <- NA
	x[x >= limits[2]] <- NA
	x
}	




#Functions for trimming extreme values  
trimtolimitz2log <- function(x)
{
    #set log values outside zscore=3 to limits
	x <- log(x)
    limits <- quantile(x, prob=c(.00455, .9545),na.rm=T)
    x[x <= limits[1]] <- limits[1]
	x[x >= limits[2]] <- limits[2]
	x <- exp(x)
	x
}	




#Functions for trimming extreme values  
trimtoNAz2log <- function(x)
{
    #set log values outside zscore=3 to limits
	x <- log(x)
    limits <- quantile(x, prob=c(.00455, .9545),na.rm=T)
    x[x <= limits[1]] <- NA
	x[x >= limits[2]] <- NA
	x <- exp(x)
	x
}	



#Functions for trimming extreme values  
trimtoNAci80log <- function(x)
{
    #set log values outside 80% limits
	x <- log(x)
    limits <- quantile(x, prob=c(.10, .90),na.rm=T)
    x[x <= limits[1]] <- NA
	x[x >= limits[2]] <- NA
	x <- exp(x)
	x
}	




#Functions for trimming extreme values  
trimtolimitci80log <- function(x)
{
    #set log values outside 80% limits
	x <- log(x)
    limits <- quantile(x, prob=c(.10, .90),na.rm=T)
    x[x <= limits[1]] <- limits[1]
	x[x >= limits[2]] <- limits[2]
	x <- exp(x)
	x
}	




#Functions for trimming extreme values  
trimtolimitz3log <- function(x)
{
    #set log values outside zscore=3 to limits
	x <- log(x)
    limits <- quantile(x, prob=c(.0027, .9973),na.rm=T)
    x[x <= limits[1]] <- limits[1]
	x[x >= limits[2]] <- limits[2]
	x <- exp(x)
	x
}	


#test1 <- rlnorm(1000, mean=2, sd=0.2)
#test2 <- trimtolimitz3log(test1)
#hist(test1)
#hist(test2)


#---------------------------------------------------------------------------------------------------------------------
#Bioequivalence - geometric mean and 90% CI
BElimits2 <- function(ratiox)
#Takes ratio of AUC (unlogged)
	{
	     #debug
		 logratio <- log(ratiox)
		 meanlog <- mean(logratio,na.rm=T)
		 n <- lengthNA(logratio)
		 semlog <- sd(logratio,na.rm=T)/sqrt(n)
		 lo90log <- meanlog - qt(0.95,n-1)*semlog   #90% CI - qt is one sided
		 hi90log <- meanlog + qt(0.95,n-1)*semlog
		 mean <- exp(meanlog)
		 lo90 <- exp(lo90log)
		 hi90 <- exp(hi90log)
		 result <- c("mean"=mean,"CI90lo"=lo90,"CI90hi"=hi90)
		 result
	}	
	
#see BElimits_test4.R for testing

  
#--------------------------------------------------------------- 
#Utility function to bind a list of dataframes - must have matching columns
  bind.list <- function(x)
   {
    #Bind a list of smaller dataframes into one big dataframe
    #Access the 2nd level of the list with [[x]]
   alldata <- x[[1]]
   for (i in 2:length(x))
    {
    alldata <- rbind( alldata, x[[i]] )
    as.data.frame(alldata)
    }
   alldata 
   }
   
   
 
 
 
 #NCA related functions--------------------------------------------------------------
 #see AUCtest.R for an evaluation of the trapz function

oneperID <- function(x)
{ # returns a single value for each ID of dataframe
  ID <- head(x, n=1)
  ID
}


lastperID <- function(x)
{ # returns a single value for each ID of dataframe
  ID <- tail(x, n=1)
  ID
}

AUCtrapz <- function(x, y)
{ # computes the integral of y with respect to x using trapezoidal integration.
  idx <- 2:length(x)
  AUC <- (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
  return(AUC)
}



AUCinfinity <- function(x, y, ntail)
{ # computes an AUC to infinity - does not handle missing data
  #debug
  #x <- c(0,10,20,30,40,50,60)
  #y <- c(16,8.2,3.8,2.1,0.9,0.51,0.24)
  #ntail <- 3
 
  AUC0t <- AUCtrapz(x,y)
    
  taildatax <- tail(x,ntail)
  taildatay <- tail(y,ntail)
  tailfit <- lm(log(taildatay) ~ taildatax)
  
  k <- -1*tailfit$coefficients["taildatax"]
  R2 <- summary(tailfit)["r.squared"]
  
  AUCinf <- NA	
  if (k > 0 & R2 > 0.4) 
    {
      Clast <- tail(y,1)
      AUCexp <- Clast/k
      AUCinf <- AUC0t+AUCexp
    }
   
  names(AUCinf) <- "AUCinf"
  AUCinf
}

Cmax <- function(x, y)
{ # computes the Cmax
 cmax <- max(y, na.rm=T)
 tindex <- which(y==cmax)
 tmax <- x[tindex]
 unique(cmax)[1] #as there can be 2 or more equal Cmax's, chose the first
}


Cmin <- function(x, y)
{ # computes the Cmin
 cmin <- min(y, na.rm=T)
 tindex <- which(y==cmin)
 tmax <- x[tindex]
 unique(cmin)[1] #as there can be 2 or more equal Cmin's, chose the first
}


tmax <- function(x, y)
{ # computes the time of Cmax
 cmax <- max(y)
 tindex <- which(y==cmax)
 tmax <- x[tindex]
 head(tmax, n=1)   #as there can be 2 or more equal Cmax's, choose the first
}


#Function for imputing last observation carried forward 
  locf <- function (x) 
{
    #Last observation carried forward
    #Finds an NA and carries forward the previous value
    good <- !is.na(x)
    positions <- seq(length(x))
    good.positions <- good * positions
    last.good.position <- cummax(good.positions)
    last.good.position[last.good.position == 0] <- NA
    x[last.good.position]
}                                      


#Function for next observation carried backward   
  nocb <- function (x) 
  {
   #Next observation carried backward
   #Reverses the locf function
   rev(locf(rev(x)))
  }  
  
#Function for calculating percent missing - use apply to do column by column  
 percent.missing <- function(x)
    {
	length(x[is.na(x)==T])/length(x)*100
	} 

#Function for taking hydraulic rate constants and calculating half-lives (2 compartment)	
	HydraulicToHalflife <- function(CL, V1, Q, V2)
		{
		k10 <- CL/V1
		k12 <- Q/V1
		k21 <- k12*V1/V2

		beta <- 0.5*(k12+k21+k10-sqrt((k12+k21+k10)^2-4*k21*k10))
		alpha <- k21*k10/beta

		thalfAlpha <- log(2)/alpha
		thalfBeta <- log(2)/beta

		result <- c("k10"=k10,"k12"=k12,"k21"=k21,"alpha"=alpha,"beta"=beta,"thalfa"=thalfAlpha,"thalfb"=thalfBeta)
		result
		}
  
  
################################################################################
# Some functions used via R2html to comment an R script and produce html output
# needs library(R2HTML)


#Echos and command only in the HTML file!
HE <- function(num.commands)
#HTML ECHO of command - echos the last "n" R commands and writes them to an HTML file
  {
   savehistory()
   command.history <- readLines(".Rhistory")
   n <- num.commands
   l <- length(command.history)
   last.commands <- command.history[(l-n):(l)]
   for (i in 1:n)
   {
   HTML(paste(">",last.commands[i],sep=""))
   HTML("")
   }
  }

#Echos and command AND output in the HTML file - brilliant!
HEO <- function(num.commands)
#HTML ECHO of command - echos the last "n" R commands and writes them to an HTML file
  {
   savehistory()
   command.history <- readLines(".Rhistory")
   n <- num.commands
   l <- length(command.history)
   last.commands <- command.history[(l-n):(l)]
   for (i in 1:n)
   {
   HTML(paste(">",last.commands[i],sep=""))
   HTML((eval(parse(text=last.commands[i]))))
   HTML("")
   }
  }

HC <- function(Hlevel,text)
#HTML COMMENT - writes a series of comments to an HTML file
#The HTML heading level is specified by Hlevel, the text is a string
#html code can be inserted and will (probably!) be understood
  {
  if (Hlevel == 1)
  HTML(paste("<H1>",text,"</H1>"))
  if (Hlevel == 2)
  HTML(paste("<H2>",text,"</H2>"))
  if (Hlevel == 3)
  HTML(paste("<H3>",text,"</H3>"))
  if (Hlevel == 4)
  HTML(paste("<H4>",text,"</H4>"))
  if (Hlevel == 5)
  HTML(paste("<p>",text,"</p>"))
  }


HG <- function(file.path)
#HTML GRAPH - inserts a png graph file into the html file
#requires master.dir from wfnviaR settings file
 {
  #debug file.path <- "1_1comp_ka_lag.g77/conc_vs_time.jpg"
  HTMLInsertGraph(GraphFileName=file.path, WidthHTML=720)
 }


HGG2 <- function(plotobj,maintext)
#HTML GGPLOT2 GRAPH - creates a png file of a ggplot2 object, then inserts the png into the html file
 {
 setwd(dir.path)
 png.file.name <- paste(maintext,".png",sep="")
 to.png(plotobj,maintext)
 HG(png.file.name)
 setwd(master.dir)
 }
 


HT <- function(file.path)
#HTML TABLE - retrieves a *.csv file as a dataframe, then writes the dataframe as an html table
 {
  #debug file.path <- "1_1comp_ka_lag.g77/1_1comp_ka_lag.fit.param.csv"
  temp.data.frame <- read.csv(file.path, stringsAsFactors=F)
  print(temp.data.frame)                                      #output in console
  HTML(temp.data.frame)                                #output in html
  rm(temp.data.frame)
 }



HL <- function(file.path)
#HTML LINK - inserts a hyperlink to a file into the html file
#requires master.dir from wfnviaR settings file
 {
 #Debug file.path <- "concdata.csv"
 #link.file <- paste(master.dir,file.path, sep="/")
 #link.file <- paste(master.dir,file.path, sep="/")
 #paste("<a href=",file.path,">",link.file,"</a>", sep=" ")
 link.string <- paste("<a href=",file.path,">",file.path,"</a>", sep=" ")
 HTML(link.string)
 #<a href="../rawdata1.xls">../rawdata1.xls</a>
 
 #HL("concdata.csv")
 }
 
 
#-------------------------------------------------------------------------------------------- 
#Function for ordering factors - from gdata package but renames to avoid conflict with reorder from stats
# This function changes the order of the levels of a factor. It can do so via three different mechanisms, depending on whether, X and FUN, new.order or sort are provided.
# If X and Fun are provided: The data in X is grouped by the levels of x and FUN is applied. The groups are then sorted by this value, and the resulting order is used for the new factor level names.
# If new.order is provided: For a numeric vector, the new factor level names are constructed by reordering the factor levels according to the numeric values. For vectors, new.order gives the list of new factor level names. In either case levels omitted from new.order will become missing (NA) values.
# If sort is provided (as it is by default): The new factor level names are generated by applying the supplied function to the existing factor level names. With sort=mixedsort the factor levels are sorted so that combined numeric and character strings are sorted in according to character rules on the character sections (including ignoring case), and the numeric rules for the numeric sections. See mixedsort for details.

#needs gdata library for mixed sort
 
 reorder2 <- function (x, X, FUN, ..., order = is.ordered(x), new.order, sort = mixedsort) 
{
    constructor <- if (order) 
        ordered
    else factor
    if (!missing(new.order)) {
        if (is.numeric(new.order)) 
            new.order <- levels(x)[new.order]
        else new.order <- new.order
    }
    else if (!missing(FUN)) 
        new.order <- names(sort(tapply(X, x, FUN, ...)))
    else new.order <- sort(levels(x))
    constructor(x, levels = new.order)
}


#--------------------------------------------------------------------------------------------  
  processSIMdata <- function(model.path.file.ext)
{ #begin processSIMdata

#Debug model.path.file.ext <- "NM14/base_CRCLCL_min_VPC.ctl"

#Work out the nonmem model dir and model file name
  model.path.file <- gsub(".ctl", "", model.path.file.ext)
  model.path <- dirname(model.path.file)
  model.file.ext <- basename(model.path.file.ext)
  model.file <- gsub(".ctl","",model.file.ext)

  model.dir <- paste(master.dir,"/",model.path, sep="")

#Work out the nonmem output dir
   output.dir <- paste(master.dir,"/",model.path.file,dir.extension, sep="")
   setwd(output.dir)
  
#Work out the SIM fit filename
  SIM.file.ext <- paste(model.file,".fit",sep="")
#Work out the SIM fit filename
  SIM.file.ext.out <- paste(model.file,".fit.csv",sep="")


#Need to remove lines with "TABLE NO. 1 in them" and the header lines for each new subject
   #Strip the unwanted lines
   indata <- readLines(SIM.file.ext)
   tablelines <- grep("TABLE NO.  1",indata)  #May be installation specific
   headerlines <- grep(" ID",indata)          #May be installation specific 

   #Extract the information in the header line for column names
   header <- indata[headerlines[1]]
   header <- scan(textConnection(header), what="char")
   colnum <- length(header)

   #Strip out the unwanted lines
   badlines <- c(tablelines,headerlines)
   indata <- indata[-badlines]

   #replace white space with commas
   for (i in 1:length(indata))
    {
    indata[i] <- gsub('[[:space:]]+',',',indata[i])
    }

   #write to a file
   writeLines(indata, "SIMtemp.txt")

   #read again as csv file
   SIMdata <- read.csv("SIMtemp.txt", header=F)
   SIMdata <- SIMdata[,-1]     #delete first blank column
   names(SIMdata) <- header

   #The NONMEM output does not make a new ID number for each simulation
   #Therefore make a list of SIM numbers
   nsims <- length(tablelines)
   numtimepoints <- length(SIMdata$ID)
   numtimespersubject <- numtimepoints/nsims
   SIMdata$SIM <- rep(1:nsims, each=numtimespersubject)

   #Write the SIMdata to a file
   write.csv(SIMdata,SIM.file.ext.out, row.names=F)

   #tidy up
   #consider deleting the original fit file?
   file.remove("SIMtemp.txt")
   setwd(master.dir)
   
   #SIMdata

} #end processSIMdata



#-----------------------------------------------------------------------------------------------------------
#Define a function that performs the calculates TAFD from DATE & TIME for one subject only
#Will be applied later to all subjects using "lapplyBy" of the doBy package
calculate.TAFD <- function(nmdata)

{ #begin calculate.TAFD

  #Convert date and time to POSIXct format in seconds
    datetimestring <- paste(nmdata$DATE, nmdata$TIME)
    datetimestring <- strptime(datetimestring, format = "%m/%d/%Y %H:%M", tz="GMT")  #see Note 1, see Note 2
    #strptime converts a string of characters representing time and date to a standard format called POSIXlt.
    nmdata$DATETIME <- as.POSIXct(datetimestring)
    #POSIXct is a standard date & time format based on seconds since the start of 1970
	
  #Make sure the data are sorted by time
    nmdata <- orderBy( ~DATETIME, data=nmdata)  

  #Find the index of the first dose
    doseindex <- which(is.na(nmdata$AMT)==F)
    firstdoseindex <- doseindex[1]
  #Calculate time elapsed since the first dose
    nmdata$TAFD <- nmdata$DATETIME - nmdata$DATETIME[firstdoseindex]
    
  #Convert result from seconds to appropriate units
    conversion.factor = 60*60     #in this case from seconds to hours, see Note 3
    nmdata$TAFD <- nmdata$TAFD/conversion.factor
	nmdata$TAFD  <- as.numeric(nmdata$TAFD) 

   #Return the modified dataframe as output from the function
       nmdata

} #end calculate.TAFD



#Function that performs the cacluates TAD from TAFD for one subject only
#Not suitable for ADDL doses - an alternative version exists for this!
#Will be applied later to all subjects using "lapplyBy" of the doBy package
calculate.TAD <- function(nmdata)
   { #begin calculate.TAD
   
     #debug
	 #nmdata <- subset(dataall, UID==1302)
	 #nmdata <- subset(nmdata, select=c(UID,DATE,TIME,RTLD,TAFD,AMT,DV,DVID,BQL,MDV))
	
   
     #Calculate dose times attributable to AMT without ADDL item
      doseindex1 <- which(is.na(nmdata$AMT)==F)
	  
      dosetimes1 <- nmdata$TAFD[doseindex1]
     
      #Combine all the dose times for AMT's with and without ADDL items
      dosetimes <- c(dosetimes1,Inf)             #Append a last dose time of infinity
      dosetimes <- sort(unique(dosetimes))       #Remove duplicate dose times and sort

      #Calculate the total number of doses for the subject
       numdoses <- length(dosetimes)
	   
      #Calculate the total number of time points for the subject
       numtimes <- length(nmdata$TAFD)
	   

     #For each dose time, cycle through the sample times, calculating DNUM and TAD
	 if (numdoses > 1) 
	 {
      for (i in 1:(numdoses-1))
       {
         for (k in 1:numtimes)
          {
		   if(nmdata$TAFD[k] >= dosetimes[i]  & nmdata$TAFD[k] < dosetimes[i+1] ) #if time is between dose1 and dose2
            {
             nmdata$DNUM[k] <- i
             nmdata$TAD[k] <- nmdata$TAFD[k] - dosetimes[i]
            }
          }
       }
	  }
	  
	  if (numdoses ==1) 
	  {
	  nmdata$DNUM <- 1
	  nmdata$TAD <- nmdata$TAFD
	  }

      #Return the modified dataframe as output from the function
       nmdata
     
   } #end calculate.TAD
 

tauANCruns <- function(TAFDin,ANCGin,ANCGrun=4)
{
#Function to look at a sequence of ANC grades, and find the runs of 4 and calculate how long they were (tau)
#Debug
	#TAFDin <- testdata$TAFD[testdata$ID==1]   #Time after first dose
	#ANCGin <- testdata$ANCG[testdata$ID==1]  #Previousll calculated neutropenia grade (0-4)
	#ANCGrun <- 4  #Grade of neutropenia for which to count durations of run

#Trap errors
if (length(ANCGin) > 1) 
{	
	
#Calculate the runs
	rundata <- rle(ANCGin)

#Turn run count into flag
	rundata$values[rundata$lengths <= 1] <- 0
	rundata$values[rundata$lengths > 1] <- 1
	RunFlag <- inverse.rle(rundata)

#Resort to good old index based computing
	TauCounter <- 0
	output <- as.data.frame(cbind(TAFDin,ANCGin,RunFlag,"RunTau"=0))

#Look at each row and processing when runflag condition passes
	 for (i in 1:(length(ANCGin)-1))
	  {
	   if (output$RunFlag[i]==1 & output$ANCGin[i]==ANCGrun)
	   {
	   TauCounter <- TauCounter + output$TAFDin[i+1]-output$TAFDin[i]
	   output$RunTau[i+1] <- TauCounter
	   }
	   else
	   {
	   TauCounter <- 0
	   output$RunTau[i] <- 0
	   }
	  } 
	#print(output)
    max(output$RunTau)  #The maximum duration of a run on ANCG > 4
}

}


  prob01234 <- function(x)
	 {
	 #probability for grade data coded as 0, 1, 2, 3 & 4
	  nsub = length(x)
      P0 <- length(x[x==0])/nsub  #probability of 0	  
	  P1 <- length(x[x>=1])/nsub  #probability of 1 or worse
	  P2 <- length(x[x>=2])/nsub  #probability of 2 or worse
	  P3 <- length(x[x>=3])/nsub  #probability of 3 or worse
	  P4 <- length(x[x>=4])/nsub  #probability of 4 or worse
	  
	  PR0 = 1-P1      # Probability of Score=0
	  PR1 = P1-P2     # Probability of Score=1
	  PR2 = P2-P3     # Probability of Score=2
      PR3 = P3-P4     # Probability of Score=3
      PR4 = P4
	  
	  c("P0"=P0,"P1"=P1,"P2"=P2,"P3"=P3,"P4"=P4,"PR0"=PR0,"PR1"=PR1,"PR2"=PR2,"PR3"=PR3,"PR4"=PR4,"nsub"=nsub)
	 }
	


#http://ggorjan.blogspot.com/2008/11/trimmed-standard-deviation.html	
sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
     if(is.complex(x)) stop("trimmed sd are not defined for complex data")
     if(trim >= 0.5) return(0)
     lo <- floor(n * trim) + 1
     hi <- n + 1 - lo
     x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}


 

# Create a data frame for input into NONMEM
  ID <- 1
  PK.TIME <- seq(from = 0,to = 24,by = 0.25)
  PD.TIME <- seq(from = 0,to = 720,by = 24)
  TIME <- c(0,PK.TIME,PD.TIME)

# Default covariate values for simulation
  CRCL <- 91.94  # Creatinine clearance (mL/min)
  HCT <- 32.5  # Haematocrit (%)
  SEX <- 0 # Male (1) or female (0)
  FFM <- 59.9  # Fat free mass (kg)
  SLC7A5 <- 0  # AA or AG (0) versus GG (1)
  GCSF <- 0  # Neupogen administered on Day 1 (0) or Day 7 (1)
  RACE <- 0  # Caucasian or unknown (0) versus African-American (1)
  ANCBASE <- 3.5 # Baseline ANC (K/microL)

# Dosing information
  INFD <- 0.5
  AMT <- 100*1.8
  RATE <- AMT/INFD

# Create input data frame
  NONMEM.data <- data.frame(
    CID = ID,
    TIME,
    AMT = c(AMT,rep(".",times = length(c(PK.TIME,PD.TIME)))),
    RATE = c(RATE,rep(".",times = length(c(PK.TIME,PD.TIME)))),
    EVID = c(1,rep(0,times = length(PK.TIME)),rep(0,times = length(PD.TIME))),
    CMT = c(1,rep(1,times = length(PK.TIME)),rep(4,times = length(PD.TIME))),
    DV = ".",
    MDV = c(1,rep(0,times = length(c(PK.TIME,PD.TIME)))),
    CRCL,
    HCT,
    SEX,
    FFM,
    SLC7A5,
    GCSF,
    RACE,
    ANCBASE
  )

  NONMEM.data <- NONMEM.data[with(NONMEM.data, order(NONMEM.data$TIME)),]
  NONMEM.data <- as.data.frame(NONMEM.data)

  write.csv(NONMEM.data,"E:/Wojciechowski/melphalan-app/Simulation/input.NONMEM.data.csv",row.names = FALSE,quote = FALSE)

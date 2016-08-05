$PROB MELPHALAN-OSU11055 TWO COMPARTMENT INFUSION PKPD
$INPUT ID TIME AMT RATE EVID CMT DV MDV CRCL HCT SEX FFM SLC7A5 GCSF RACE ANCBASE
$DATA ..\input.NONMEM.data.csv IGNORE=C

$SUBR ADVAN6 TOL=9

$MODEL
      COMP=(PKCENTR)
      COMP=(PKPERI)
       COMP=(STEM)
       COMP=(ANC,DEFOBS)
       COMP=(TRANSIT1)
       COMP=(TRANSIT2)
       COMP=(TRANSIT3)
    COMP=(INPUT) ;CAPACITY DIRECTLY INDUCE ANC, NOT THROUGH DELAYED TRANSIT


$PK

   TVCL=THETA(1)*((CRCL/91.94)**THETA(6))*((FFM/59.90)**(0.75))*((HCT/32.50)**THETA(8))
   TVV1=THETA(2)*(FFM/59.90)
   TVQ=THETA(3)*((FFM/59.90)**(0.75))
   TVV2=THETA(4)*(1+SLC7A5*THETA(7))*(FFM/59.90)

   CL=TVCL*EXP(ETA(1))
   V1=TVV1*EXP(ETA(2))
   Q=TVQ*EXP(ETA(3))
   V2=TVV2*EXP(ETA(1)*THETA(5))

   S1=V1

   K10=CL/V1
   K12=Q/V1
   K21=Q/V2


   BASE = THETA(10)*EXP(ETA(4))*((ANCBASE/3.5)**THETA(20))*((HCT/32.50)**THETA(26))
   SLOPE = THETA(11)*EXP(ETA(5))*(1+GCSF*THETA(21))*(1+SEX*THETA(24))*((HCT/32.50)**THETA(25))
   MTT  = THETA(12)*EXP(ETA(6))*(1+GCSF*THETA(19))*(1+Race*THETA(22))*((CRCL/91.94)**THETA(23))
   POWER1EN = THETA(13)
   POWER1EX = THETA(14)


   K    = 4/MTT

   KE = LOG(2)/THETA(15)

   ON=1
   IF (GCSF.EQ.1) ON=0
   TVIP0=ON*THETA(16)+(1-ON)*THETA(17)
   IP0=TVIP0*EXP(ETA(7))
;  IS PRECURSOR UNEXPRESSED BASELINE IN A3, AS A CAPACITY TO INDUCE ANC RAPIDLY BY GCSF

   IPT = THETA(18)*EXP(ETA(8))

   A_0(3)  = KE*BASE/K
   A_0(4)  = BASE
   A_0(5)  = KE*BASE/K
   A_0(6)  = KE*BASE/K
   A_0(7)  = KE*BASE/K
   A_0(8)  = IP0


$DES
CP      =  A(1)/V1
DRUG    =  SLOPE*CP
;****************KINETICS****************************
DADT(1)=A(2)*K21-A(1)*K10-A(1)*K12
DADT(2)=A(1)*K12-A(2)*K21
;****************DYNAMICS****************************
QN=0
IF(GCSF.EQ.0.AND.TIME.GT.72) QN=1
IF(GCSF.EQ.1.AND.TIME.GT.216) QN=1
POWER1 = POWER1EN + QN*POWER1EX

FN = (BASE/0.001)**POWER1
IF(A(4).GT.0.001) FN = (BASE/A(4))**POWER1

KIN=0
IF(GCSF.EQ.0.AND.TIME.GT.72) KIN=1/IPT
IF(GCSF.EQ.1.AND.TIME.GT.216) KIN=1/IPT

DADT(3) =  K*A(3)*(1-DRUG)*FN - K*A(3)    ;STEM
DADT(4) =  K*A(7) - KE*A(4) + KIN*A(8)   ;ANC
DADT(5) =  K*A(3) - K*A(5)              ;TRANSIT 1
DADT(6) =  K*A(5) - K*A(6)             ;TRANSIT 2
DADT(7) =  K*A(6) - K*A(7)               ;TRANSIT 3
DADT(8) =  -KIN*A(8)             ; INPUT1


$ERROR
   IPRED=F
IF (CMT .EQ. 4) IPRED=((F)**0.2-1)/0.2 ; BOXCOX
   IRES   = DV-IPRED

  W=THETA(9)

   IF(W.EQ.0) W=1
   SD=SQRT(W**2)
   IWRES  = IRES/SD

QK=0
QD=0
IF (CMT .EQ. 1) QK=1
IF (CMT .EQ. 4) QD=1


PKY = IPRED * EXP(ERR(1))
PDY = IPRED + SD*ERR(2)

 Y=QK*PKY+QD*PDY

SET=IREP

$THETA
28.9 FIX ; POPCL
19.5 FIX  ; POPV1
28.5 FIX  ; POPQ
20.9 FIX   ; POPV2
1.31 FIX ; THETA5
0.193 FIX; THETA6
-0.130 FIX; THETA7
0.216 FIX; THETA8

0.693 FIX ; THETA9 W
5.66 FIX   ;THETA10 BASE
7.41 FIX         ;THETA11 SLOPE
96 FIX	;THETA12 MTT
0.188 FIX     ;THETA13 EndoGCSFPOWER1
0.0313 FIX ; THETA14 NeupogenPOWER1
7 FIX ; THETA15 KE halflife
109 FIX; THETA16 IP01
0.0564 FIX; THETA17 IP02
14.7 FIX; THETA18 IPT
0.157 FIX; THETA19 FAC GCSF_MTT
0.226 FIX; THETA20 FAC2 ANC_BASE
0.421 FIX; THETA21 FAC3 GCSF_SLOPE
-0.0894 FIX; THETA22 FAC4 Race_MTT
-0.056 FIX; THETA23 FAC5 CrCL_MTT
0.214 FIX; THETA24 FAC6 SEX_SLOPE
0.646 FIX; THETA25 FAC7 HCT_SLOPE
0.582 FIX; THETA26 FAC8 HCT_BASE

$OMEGA
0.0792 FIX ;PPVCL
0.0740 FIX  ;PPVV1
0.297 FIX ;PPVQ
0.1 FIX  ;ETA4 BASE
0.0631 FIX  ;ETA5 SLOPE
0.00561 FIX  ;ETA6 MTT
0.23 FIX ;ETA7 IP0
0 FIX ;ETA8 IPT


$SIGMA
0.0584 FIX ; EPS1
1 FIX ; EPS2

$SIMULATION (9215690) ONLYSIM SUBPROBLEMS=1000

$TABLE SET GCSF ID TIME CMT Y IPRED
ONEHEADER NOPRINT FILE=*.fit

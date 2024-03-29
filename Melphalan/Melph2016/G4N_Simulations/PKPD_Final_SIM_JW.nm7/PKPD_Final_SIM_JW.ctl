$PROBLEM MELPHALAN-OSU11055 TWO COMPARTMENT INFUSION PKPD
$INPUT ID TIME ANC DV AMT RATE MDV DVID CMT WT CRCL FFM HCT GCSF SLC7A5 BUN ANCBASE SEX WBC LNP53FOLD
$DATA ..\test_data_1.csv IGNORE=#

$SUBR ADVAN6 TOL = 5

$MODEL 
 COMP=(PKCENTR)
 COMP=(PKPERI)
 COMP=(STEM)
 COMP=(ANC,DEFOBS)
 COMP=(TRANSIT1)
 COMP=(TRANSIT2)
 COMP=(TRANSIT3)    
 COMP=(INPUT) ;CAPACITY DIRECTLY INDUCE ANC, NOT THROUGH DELAYED TRANSIT  
 COMP=(G4NDUR) ;COMPARTMENT FOR COLLECTING TIME SPENT IN GRADE 4 NEUTROPENIA

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
  
 BASE=THETA(10)*EXP(ETA(4))*((ANCBASE/3.5)**THETA(20))*((HCT/32.50)**THETA(26))
 SLOPE=THETA(11)*EXP(ETA(5))*(1+GCSF*THETA(21))*(1+SEX*THETA(23))*((HCT/32.50)**THETA(24))
 MTT=THETA(12)*EXP(ETA(6))*(1+GCSF*THETA(19))*((BUN/14)**THETA(22))
 POWER1EN=THETA(13)
 POWER1EX=THETA(14)
 
 K=4/MTT
 
 KE=LOG(2)/THETA(15)

 IF(LNP53FOLD.EQ.-99) THEN
  FCOV0=1
 ELSE 
  FCOV0=(LNP53FOLD/2.62)**THETA(25)
 ENDIF

 ON=1
 IF(GCSF.EQ.1) ON=0
 TVIP0=ON*THETA(16)+(1-ON)*THETA(17)
 IP0=TVIP0*EXP(ETA(7))*FCOV0*((WBC/4.90)**THETA(27))  
 IPT=THETA(18)*EXP(ETA(8))
 
 A_0(3)=KE*BASE/K
 A_0(4)=BASE
 A_0(5)=KE*BASE/K
 A_0(6)=KE*BASE/K
 A_0(7)=KE*BASE/K      
 A_0(8)=IP0 

$DES
 CP=A(1)/V1
 DRUG=SLOPE*CP
 
 ;****************KINETICS****************************
 DADT(1)=A(2)*K21-A(1)*K10-A(1)*K12
 DADT(2)=A(1)*K12-A(2)*K21
 
 ;****************DYNAMICS****************************
 QN=0
 IF(GCSF.EQ.0.AND.TIME.GT.72) QN=1
 IF(GCSF.EQ.1.AND.TIME.GT.216) QN=1
 POWER1=POWER1EN + QN*POWER1EX

 FN=(BASE/0.001)**POWER1
 IF(A(4).GT.0.001) FN=(BASE/A(4))**POWER1

 KIN=0
 IF(GCSF.EQ.0.AND.TIME.GT.72) KIN=1/IPT
 IF(GCSF.EQ.1.AND.TIME.GT.216) KIN=1/IPT

 DADT(3)=K*A(3)*(1-DRUG)*FN -K*A(3) ;STEM
 DADT(4)=K*A(7) -KE*A(4) +KIN*A(8) ;ANC 
 DADT(5)=K*A(3) -K*A(5) ;TRANSIT 1
 DADT(6)=K*A(5) -K*A(6) ;TRANSIT 2
 DADT(7)=K*A(6) -K*A(7) ;TRANSIT 3 
 DADT(8)=-KIN*A(8) ; INPUT1

 DADT(9)=0
 IF(A(4).LT.0.5) DADT(9)=1

$ERROR 
 G4ND=A(9)

 IPRED=F
 IF (CMT.EQ.4) IPRED=((F)**0.2-1)/0.2 ;BOXCOX
 IRES=DV-IPRED
 
 W=THETA(9) 
 IF(W.EQ.0) W=1
 SD=SQRT(W**2)
 IWRES=IRES/SD

 QK=0
 QD=0
 IF(CMT.EQ.1) QK=1
 IF(CMT.EQ.4) QD=1

 PKY=IPRED * EXP(ERR(1))
 PDY=IPRED + SD*ERR(2)
  
 Y=QK*PKY+QD*PDY    
  
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
5.65 FIX    ;THETA10 BASE
7.4 FIX   ;THETA11 SLOPE
95.8 FIX  ;THETA12 MTT
0.182 FIX     ;THETA13 ENDOGCSFPOWER1
0.0364 FIX ; THETA14 NEUPOGENPOWER1
7 FIX ; THETA15 KE HALFLIFE
110 FIX; THETA16 IP01
0.035 FIX; THETA17 IP02
14.7 FIX; THETA18 IPT
0.153 FIX; THETA19 FAC GCSF_MTT
0.225 FIX; THETA20 FAC2 ANC_BASE
0.405 FIX; THETA21 FAC3 GCSF_SLOPE
0.0629 FIX; THETA22 FAC4 BUN_MTT
0.208 FIX; THETA23 FAC5 SEX_SLOPE
0.653 FIX; THETA24 FAC6 HCT_SLOPE
-0.777 FIX; THETA25 FAC7 LNP53_IP0
0.584 FIX; THETA26 FAC8 HCT_BASE
0.411 FIX; THETA27 FAC9


$OMEGA
0.0792 FIX ;PPVCL
0.0740 FIX  ;PPVV1
0.297 FIX ;PPVQ
0.101 FIX  ;ETA4 BASE
0.0608 FIX  ;ETA5 SLOPE
0.00674 FIX  ;ETA6 MTT
0.159 FIX ;ETA7 IP0
0 FIX ;ETA8 IPT

$SIGMA
0.0584 FIX ; EPS1
1 FIX ; EPS2

$SIMULATION (9215690) ONLYSIM SUBPROBLEMS=100

$TABLE ID TIME AMT MDV DVID CMT Y IPRED PRED G4ND
WT CRCL FFM HCT GCSF SLC7A5 BUN ANCBASE SEX WBC LNP53FOLD
NOPRINT ONEHEADER FILE=pkpd_final_sim_jw.fit 

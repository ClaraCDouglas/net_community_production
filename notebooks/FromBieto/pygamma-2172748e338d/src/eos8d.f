C APRIL 12 1984
C  EOS80 DERIVATIVES TEMP. & SALT
      double precision FUNCTION EOS8D(S,T,P0,DRV)
C  MODIFIED RCM
C ******************************************************
C SPECIFIC VOLUME ANOMALY (STERIC ANOMALY) BASED ON 1980 EQUATION
C OF STATE FOR SEAWATER AND 1978 PRACTICAL SALINITY SCALE.
C REFERENCES
C MILLERO, ET AL (1980) DEEP-SEA RES.,27A,255-264
C MILLERO AND POISSON 1981,DEEP-SEA RES.,28A PP 625-629.
C BOTH ABOVE REFERENCES ARE ALSO FOUND IN UNESCO REPORT 38 (1981)
C UNITS:
C       PRESSURE        P0       DECIBARS
C       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
C       SALINITY        S        (IPSS-78)
C       SPEC. VOL. ANA. EOS8D    M**3/KG *1.0E-8
C       DENSITY ANA.    SIGMA    KG/M**3
C DRV MATRIX FORMAT
C     1    2     3
C  1   V  ,VT  ,VTT    TEMP DERIV. S,T,P
C  2   V0 ,VOT ,V0TT   FOR S,T,0
C  3   RO ,ROT ,ROTT   FOR S,T,P  DENSITY DERIV
C  4   K0 ,K0T ,K0TT   FOR S,T,0 SEC BULK MOD
C  5   A  ,AT  ,ATT
C  6   B  ,BT  ,BTT    BULK MOD PRESS COEFFS
C  7 DRDP ,K   ,DVDP   PRESSURE DERIVATIVE
C  8   R0S,    ,VS      SALINITY DERIVATIVES
C
CHECK VALUE: FOR S = 40 (IPSS-78) , T = 40 DEG C, P0= 10000 DECIBARS.
C        DR/DP                  DR/DT                 DR/DS
C       DRV(1,7)              DRV(2,3)             DRV(1,8)
C
C FINITE DIFFERENCE WITH 3RD ORDER CORRECTION DONE IN DOUBLE PRECSION
C
C       3.46969238E-3       -.43311722           .705110777
C
C EXPLICIT DIFFERENTIATION SINGLE PRECISION FORMULATION EOS80
C
C       3.4696929E-3        -.4331173            .7051107
C
C *******************************************************
C
      IMPLICIT double precision (A-Z)
      double precision K,K0,KW,K35
C      REAL P,T,S,SIG,SR,R1,R2,R3,R4
C      REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW,K35
      double precision DRV(3,8)
C ********************
C DATA
      DATA R3500,R4/1028.1063,4.8314D-4/
      DATA DR350/28.106331/
C   R4 IS REFERED TO AS  C  IN MILLERO AND POISSON 1981
C CONVERT PRESSURE TO BARS AND TAKE SQUARE ROOT SALINITY.
      P=P0/10.
      R3500=1028.1063
      SAL=S
      SR = SQRT(ABS(S))
C *********************************************************
C PURE WATER DENSITY AT ATMOSPHERIC PRESSURE
C   BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
C
      R1 = ((((6.536332d-9*T-1.120083d-6)*T+1.001685d-4)*T
     X-9.095290d-3)*T+6.793952d-2)*T-28.263737
C SEAWATER DENSITY ATM PRESS.
C  COEFFICIENTS INVOLVING SALINITY
C  R2 = A   IN NOTATION OF MILLERO AND POISSON 1981
      R2 = (((5.3875d-9*T-8.2467d-7)*T+7.6438d-5)*T-4.0899d-3)*T
     X+8.24493d-1
C  R3 = B  IN NOTATION OF MILLERO AND POISSON 1981
      R3 = (-1.6546d-6*T+1.0227d-4)*T-5.72466d-3
C  INTERNATIONAL ONE-ATMOSPHERE EQUATION OF STATE OF SEAWATER
      SIG = (R4*S + R3*SR + R2)*S + R1
C SPECIFIC VOLUME AT ATMOSPHERIC PRESSURE
      V350P = 1.0/R3500
      SVA = -SIG*V350P/(R3500+SIG)
      SIGMA=SIG+DR350
      DRV(1,3) = SIGMA
      V0 = 1.0/(1000.0 + SIGMA)
      DRV(1,2) = V0
C COMPUTE DERIV WRT SALT OF RHO
      R4S=9.6628d-4
      RHOS=R4S*SAL+1.5*R3*SR+R2
C*************************************
C COMPUTE DERIV WRT TEMP OF RHO
      R1 =(((3.268166d-8*T-4.480332d-6)*T+3.005055d-4)*T
     X -1.819058d-2)*T+6.793952d-2
      R2 = ((2.155d-8*T-2.47401d-6)*T+1.52876d-4)*T-4.0899d-3
      R3 = -3.3092d-6*T+1.0227d-4
C
      RHOT = (R3*SR + R2)*SAL + R1
      DRDT=RHOT
      DRV(2,3) = RHOT
      RHO1 = 1000.0 + SIGMA
      RHO2 = RHO1*RHO1
      V0T = -RHOT/(RHO2)
C*********** SPECIFIC VOL. DERIV WRT S ***********
      V0S=-RHOS/RHO2
C****************************
      DRV(1,8)=RHOS
      DRV(2,2) = V0T
C COMPUTE SECOND DERIVATIVE OF RHO
      R1 = ((1.3072664d-7*T-1.3440996d-5)*T+6.01011d-4)*T-1.819058d-2
      R2 = (6.465d-8*T-4.94802d-6)*T+1.52876d-4
      R3 = -3.3092d-6
C
      RHOTT = (R3*SR + R2)*SAL + R1
      DRV(3,3) = RHOTT
      V0TT = (2.0*RHOT*RHOT/RHO1 - RHOTT)/(RHO2)
      DRV(3,2) = V0TT
C  SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8
      EOS8D=SVAN
C ******************************************************************
C ******  NEW HIGH PRESSURE EQUATION OF STATE FOR SEAWATER ********
C ******************************************************************
C        MILLERO, ET AL , 1980 DSR 27A, PP 255-264
C               CONSTANT NOTATION FOLLOWS ARTICLE
C********************************************************
C COMPUTE COMPRESSION TERMS
      E = (9.1697d-10*T+2.0816d-8)*T-9.9348d-7
      BW = (5.2787d-8*T-6.12293d-6)*T+3.47718d-5
      B = BW + E*S
C
C******* DERIV B WRT SALT
      DBDS=E
C************************
C  CORRECT B FOR ANAMOLY BIAS CHANGE
      DRV(1,6) = B + 5.03217d-5
C DERIV OF B
      BW = 1.05574d-7*T-6.12293d-6
      E = 1.83394d-9*T +2.0816d-8
      BT = BW + E*SAL
      DRV(2,6) = BT
C   COEFFICIENTS OF A
C SECOND DERIV OF B
      E = 1.83394d-9
      BW = 1.05574d-7
      BTT = BW + E*SAL
      DRV(3,6) = BTT
      D = 1.91075d-4
      C = (-1.6078d-6*T-1.0981d-5)*T+2.2838d-3
      AW = ((-5.77905d-7*T+1.16092d-4)*T+1.43713d-3)*T
     X-0.1194975
      A = (D*SR + C)*S + AW
C
C  CORRECT A FOR ANAMOLY BIAS CHANGE
      DRV(1,5) = A + 3.3594055
C*****   DERIV A WRT SALT ************
      DADS=2.866125d-4*SR+C
C************************************
C DERIV OF A
      C = -3.2156d-6*T -1.0981d-5
      AW = (-1.733715d-6*T+2.32184d-4)*T+1.43713d-3
C
      AT = C*SAL + AW
      DRV(2,5) = AT
C SECOND DERIV OF A
      C = -3.2156d-6
      AW = -3.46743d-6*T + 2.32184d-4
C
      ATT = C*SAL + AW
      DRV(3,5) = ATT
C COEFFICIENT K0
      B1 = (-5.3009d-4*T+1.6483d-2)*T+7.944d-2
      A1 = ((-6.1670d-5*T+1.09987d-2)*T-0.603459)*T+54.6746
      KW = (((-5.155288d-5*T+1.360477d-2)*T-2.327105)*T
     X+148.4206)*T-1930.06
      K0 = (B1*SR + A1)*S + KW
C     ADD BIAS TO OUTPUT K0 VALUE
      DRV(1,4) = K0+21582.27
C****** DERIV K0 WRT SALT ************
      K0S=1.5*B1*SR+A1
C*************************************
C     DERIV K WRT SALT   *************
       KS=(DBDS*P+DADS)*P+K0S
C***********************************
C DERIV OF K0
      B1 = -1.06018d-3*T+1.6483d-2
C APRIL 9 1984 CORRECT A1 BIAS FROM -.603457 !!!
      A1 = (-1.8501d-4*T+2.19974d-2)*T-0.603459
      KW = ((-2.0621152d-4*T+4.081431d-2)*T-4.65421)*T+148.4206
      K0T = (B1*SR+A1)*SAL + KW
      DRV(2,4) = K0T
C SECOND DERIV OF K0
      B1 = -1.06018d-3
      A1 = -3.7002d-4*T + 2.19974d-2
      KW = (-6.1863456d-4*T+8.162862d-2)*T-4.65421
      K0TT = (B1*SR + A1)*SAL + KW
      DRV(3,4) = K0TT
C
C
C EVALUATE PRESSURE POLYNOMIAL
C ***********************************************
C   K EQUALS THE SECANT BULK MODULUS OF SEAWATER
C   DK=K(S,T,P)-K(35,0,P)
C  K35=K(35,0,P)
C ***********************************************
      DK = (B*P + A)*P + K0
      K35  = (5.03217d-5*P+3.359406)*P+21582.27
      GAM=P/K35
      PK = 1.0 - GAM
      SVA = SVA*PK + (V350P+SVA)*P*DK/(K35*(K35+DK))
C  SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8
      EOS8D=SVAN
      V350P = V350P*PK
C  ****************************************************
C COMPUTE DENSITY ANAMOLY WITH RESPECT TO 1000.0 KG/M**3
C  1) DR350: DENSITY ANAMOLY AT 35 (IPSS-78), 0 DEG. C AND 0 DECIBARS
C  2) DR35P: DENSITY ANAMOLY 35 (IPSS-78), 0 DEG. C ,  PRES. VARIATION
C  3) DVAN : DENSITY ANAMOLY VARIATIONS INVOLVING SPECFIC VOL. ANAMOLY
C ********************************************************************
C CHECK VALUE: SIGMA = 59.82037  KG/M**3 FOR S = 40 (IPSS-78),
C T = 40 DEG C, P0= 10000 DECIBARS.
C *******************************************************
      DR35P=GAM/V350P
      DVAN=SVA/(V350P*(V350P+SVA))
      SIGMA=DR350+DR35P-DVAN
      DRV(1,3)=SIGMA
      K=K35+DK
      VP=1.0-P/K
      KT = (BT*P + AT)*P + K0T
      KTT = (BTT*P + ATT)*P + K0TT
C
      V = 1.0/(SIGMA+1000.0D0)
      DRV(1,1) = V
      V2=V*V
C  DERIV SPECIFIC VOL. WRT SALT **********
      VS=V0S*VP+V0*P*KS/(K*K)
      RHOS=-VS/V2
C***************************************
      DRV(3,8)=VS
      DRV(1,8)=RHOS
C
      VT = V0T*VP + V0*P*KT/(K*K)
      VTT = V0TT*VP+P*(2.0*V0T*KT+KTT*V0-2.0*KT*KT*V0/K)/(K*K)
      R0TT=(2.0*VT*VT/V-VTT)/V2
      DRV(3,3)=R0TT
      DRV(2,1) = VT
      DRV(3,1) = VTT
      RHOT=-VT/V2
      DRDT=RHOT
      DRV(2,3)=RHOT
C PRESSURE DERIVATIVE DVDP
C   SET A & B TO UNBIASED VALUES
      A=DRV(1,5)
      B=DRV(1,6)
      DKDP = 2.0*B*P + A
C CORRECT DVDP TO PER DECIBAR BY MULTIPLE *.1
      DVDP = -.1*V0*(1.0 - P*DKDP/K)/K
      DRV(1,7) = -DVDP/V2
      DRV(2,7) = K
      DRV(3,7) = DVDP
      RETURN
      END

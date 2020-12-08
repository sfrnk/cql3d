      SUBROUTINE ADCECI(PTE, PNE)
C
C*********************************************************************
C     CALCULATE THE ELECTRON COLLISIONAL IONIZATION RATE USING
C     THE FORMULAS SELECTED BY THE LECI SWITCH.  LECIOK RETURNS
C     THE RATE FORMULAS ACTUALLY USED, WHICH MAY DIFFER FROM
C     LECI IF LECI WAS SET INAPPROPRIATELY FOR THE GIVEN ELEMENT.
C*********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON / ADSDAT / NUCZ, NSPC, APN(100,10), AQN(100,10),
     &                  NVALNC(100), SIGMA(10,10), QN(100,10),
     &                  EN(100,10), EIN(100,5), ENM(100,5,10),
     &                  FNM(100,5,10), FNN(100), ENN(100)
C
      COMMON / ADELEC / RCLION(100), RRAREC(100), RDIREC(100),
     &                  CIZLOS(100), RADRRC(100), RADDRC(100),
     &                  RADCLX(100), RADBRM(100), TOTREC(100),
     &                  TOTRAD(100)
C
      COMMON / ADPARM / LADEN, LADTIP, LTIPOK, LECI, LECIOK
C
      EXTERNAL ADCXECI,ADCBFI,ADCYGR
C
C     ILMAX = CURRENT MAXIMUM NUMBER OF RATE OPTIONS
C
      ILMAX = 3
C
C     CHECK FOR VALID LECI; RETURN LECI = 1 OPTION IF NOT VALID
C
      IF(LECI .LT. 1 .OR. LECI .GT. ILMAX) GO TO 10
C
      GO TO (10, 20, 30) LECI
C
C     LECI = 1; USE XSNQ DERIVED IONIZATION RATE CALCULATION
C
   10 CALL ADCXECI(PTE, PNE)
      LECIOK = 1
      RETURN
C
C     LECI = 2; USE BELFAST GROUP LIGHT ELEMENT RATE FITS
C
   20 CALL ADCBFI(NUCZ, PTE, PNE, RCLION, CIZLOS, IFLAG)
      IF(IFLAG .NE. 0) GO TO 10
      LECIOK = 2
      RETURN
C
C     LECI = 3; USE YOUNGER'S RATES FOR AVAILABLE ELEMENTS AND SHELLS,
C               XSNQ RATES OTHERWISE ON A CHARGE STATE-BY-CHARGE STATE
C               BASIS
C
   30 CALL ADCXECI(PTE, PNE)
      CALL ADCYGR(NUCZ, PTE, PNE, RCLION, CIZLOS, IFLAG)
      LECIOK = 3
      IF(IFLAG .NE. 0) LECIOK = 1
      RETURN

      END



      SUBROUTINE ADCXECI(PTE, PNE)
C
C
C***********************************************************************
C     REV: 6/15/78
C
C     COMPUTES THE XSNQ-DERIVED COLLISIONAL IONIZATION RATES RCLION
C     AT THE ELECTRON TEMPERATURE PTE (KEV), WHERE JQ=(IONIC CHARGE+1).
C     THESE RATES HAVE THE UNITS (SEC-1).
C
C     NOTES:
C
C     1) ELECTRONS ARE ASSUMED TO FILL IONIC SHELLS IN SIMPLE ASCENDING
C        ORDER.  THIS IS ASSUMED BOTH EXPLICITLY (THE APN ARRAY)
C        AS WELL AS IMPLICITLY THROUGH THE USE OF NVALNC VALUES
C        IN LOOP INDICES, ETC.
C     2) THROUGH USE OF THE FUNCTION EXPUND, EXP(X) IS SET EQUAL TO
C        ZERO IN THE RATE EQUATIONS FOR VALUES OF X < -45.
C     3) RATES FOR THE SPECIES JQ REFER TO TRANSFORMATIONS OF
C        THAT SPECIES INTO ADJACENT IONIZATION STATES, NOT THE
C        RATES AT WHICH THAT SPECIES IS FORMED.
C
C
C     PRINCIPAL PHYSICS REFERENCES:
C
C     1) POST, JENSEN, TARTER, GRASBERGER, LOKKE, 'STEADY STATE
C        RADIATIVE COOLING RATES FOR LOW-DENSITY HIGH-TEMPERATURE
C        PLASMAS', PPPL-1352 (SUB. TO N.F.) (WITH CORRECTIONS)
C     2) LOKKE AND GRASBERGER, 'XSNQ-U A NON-LTE EMISSION AND
C        ABSORPTION COEFFICIENT SUBROUTINE', LLL REPORT UCRL-52276
C        (WITH CORRECTIONS).
C     3) XSNQ FORTRAN LISTING
C
C*************************************************************************
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
CMPIINSERT_INCLUDE     
C
      COMMON / ADSDAT / NUCZ, NSPC, APN(100,10), AQN(100,10),
     &                  NVALNC(100), SIGMA(10,10), QN(100,10),
     &                  EN(100,10), EIN(100,5), ENM(100,5,10),
     &                  FNM(100,5,10), FNN(100), ENN(100)
C
      COMMON / ADELEC / RCLION(100), RRAREC(100), RDIREC(100),
     &                  CIZLOS(100), RADRRC(100), RADDRC(100),
     &                  RADCLX(100), RADBRM(100), TOTREC(100),
     &                  TOTRAD(100)
C
      EXTERNAL ADCEUND
C
      ZJLKEV = 1.6021E-16
C
C
      ZTE = PTE
C
      DO 100 JQ = 1 , NUCZ
C
      ZCHG = DBLE(JQ - 1)
      IVALNC = NVALNC(JQ)
C
C
C
C     ******************** IONIZATION RATE ********************
C
C
C     CALCULATE ELECTRON COLLISIONAL IONIZATION RATES PER SHELL
C     AND SUM FOR ALL POPULATED SHELLS. ALSO COMPUTE
C     IONIZATION ENERGY LOSS RATE CIZLOS.
C
      ZISUM = 0.d0 !0.
      ZESUM = 0.d0 !0.0
C
      DO 20 JN = 1 , IVALNC
      ZJN = JN
      ZXN = EIN(JQ,JN)/ZTE
C
      ZY = DLOG10(ZXN/4.)
      ZFY = .23 + .046 * ZY + .1074 * ZY**2
     &          - .0459 * ZY**3 - .01505 * ZY**4
      ZGAUNT = 12.18 * DEXP(-ZJN / (ZJN+5.)) * (1. + .0335*ZJN)
     &         * (1. - .622/(ZCHG+1.) - .0745/(ZCHG+1.)**2)
     &         * DABS(ZFY)
      ZITEMP = APN(JQ,JN) * ADCEUND(-ZXN,-45.d0) *
     &         (1.-ADCEUND(-ZXN,-45.d0)) * ZGAUNT / EIN(JQ,JN)**2
      !YuP[2020-11-16] Fixed BUG: Cannot use 45. in ADCEUND(-ZXN,-45.)
      !The argument should be real*8 !!!
      ZISUM = ZISUM + ZITEMP
      ZESUM = ZESUM + EIN(JQ,JN) * ZITEMP
   20 CONTINUE
C
      ZISUM = ZISUM * DSQRT(ZTE)
      ZESUM = ZESUM * DSQRT(ZTE)
      IF(ZISUM .LT. 1.d-26) ZISUM = 0.d0 !0.
      IF(ZESUM .LT. 1.d-26) ZESUM = 0.d0 !0.0
      RCLION(JQ) = (PNE*3.44E-11) * ZISUM
      CIZLOS(JQ) = (PNE*3.44E-11) * ZJLKEV * ZESUM
C
C
  100 CONTINUE
!--CMPIINSERT_IF_RANK_EQ_0      
!      WRITE(*,*)' ADCXECI: RCLION(1:3)=',RCLION(1:3)
!      WRITE(*,*)' ADCXECI: PNE,ZISUM=',PNE,ZISUM
!--CMPIINSERT_ENDIF_RANK
C
C
      RCLION(NUCZ + 1) = 0.
      CIZLOS(NUCZ + 1) = 0.0
C
C
      RETURN
      END SUBROUTINE ADCXECI



      SUBROUTINE ADCBFI(KNUCZ, PTE, PNE, PBFI, PBFIL, KFLAG)
C
C*********************************************************************
C     CALCULATE THE ELECTRON COLLISIONAL IONIZATION RATE COEFFICIENTS
C     USING FITTING FORMULAS FROM THE BELFAST GROUP, CULHAM
C     REPORT CLM-R216 (DECEMBER 1981):
C
C          ATOMIC AND MOLECULAR DATA FOR FUSION, PART I
C
C          RECOMMENDED CROSS SECTIONS AND RATES FOR ELECTRON
C          IONISATION OF LIGHT ATOMS AND IONS
C
C          K L BELL, H B GILBODY, J G HUGHES, A E KINGSTON, F J SMITH
C          THE QUEEN'S UNIVERSITY OF BELFAST, BELFAST, N IRELAND
C
C     KNUCZ = ELEMENT ATOMIC NUMBER ( Z=1 TO Z=8 ONLY)
C     PTE = ELECTRON TEMPERATURE (KEV)
C     PNE = ELECTRON DENSITY (CM-3)
C
C     PBFI = IONIZATION RATES FOR THE Z+1 CHARGE STATES (SEC-1)
C     PBFIL = ARRAY OF IONIZATION ENERGY LOSS RATES (WATTS/ION)
C     KFLAG = 0 FOR O.K., = 1 FOR INVALID ELEMENT
C
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      COMMON / COMBFI / EH(1),  AH(6,1),  AAH(1),  BTH(3,1),
     &                  EHE(2), AHE(6,2), AAHE(2), BTHE(3,2),
     &                  ELI(3), ALI(6,3), AALI(3), BTLI(3,3),
     &                  EBE(4), ABE(6,4), AABE(4), BTBE(3,4),
     &                  EB(5),  AB(6,5),  AAB(5),  BTB(3,5),
     &                  EC(6),  AC(6,6),  AAC(6),  BTC(3,6),
     &                  EN(7),  AN(6,7),  AAN(7),  BTN(3,7),
     &                  EO(8),  AO(6,8),  AAO(8),  BTO(3,8)
C
C
      DIMENSION PBFI(*), PBFIL(*)
C
C
C     HYDROGEN
C
      DATA EH / 13.60 /
      DATA AH / 2.3742E-8, -3.6866E-9, -1.0366E-8,
     &         -3.8010E-9,  3.4159E-9,  1.6834E-9 /
      DATA AAH / 2.4617E-8 /
      DATA BTH / 9.5986E-8, -9.2463E-7, 3.9973E-6 /
C
C     HELIUM
C
      DATA EHE / 24.6, 54.42 /
      DATA AHE / 1.4999E-8, 5.6656E-10, -6.0821E-9,
     1          -3.5894E-9, 1.5529E-9,   1.3207E-9,
     2           3.4356E-9, -1.6865E-9, -6.9236E-10,
     2           9.7863E-11, 1.5591E-10, 6.2236E-11 /
      DATA AAHE / 3.1373E-8, 3.0772E-9 /
      DATA BTHE / 4.7893E-8, -7.7359E-7, 3.7366E-6,
     2            1.1902E-8, -1.1514E-7, 5.0489E-7 /
C
C     LITHIUM
C
      DATA ELI / 5.39, 75.64, 122.45 /
      DATA ALI / 9.9655E-8, -5.5941E-8, -5.5228E-8,
     1           4.0589E-8, 1.4800E-8, -1.3120E-8,
     2           3.4023E-9, -7.6588E-10, -8.6078E-10,
     2           -8.9748E-10, 4.1661E-10, 3.3188E-10,
     3           1.1786E-9, -8.7637E-10, -9.3373E-11,
     3           2.1137E-10, 1.9017E-11, -4.0679E-11 /
      DATA AALI / 4.5456E-8, 7.3504E-9, 1.9767E-9 /
      DATA BTLI / 2.7800E-7, -1.5830E-6, 5.4650E-6,
     2            5.3459E-10, -5.6387E-8, 2.9577E-7,
     3           -1.0926E-9, 8.8700E-10, 6.0764E-9 /
C
C     BERYLLIUM
C
      DATA EBE / 9.32, 18.21, 153.89, 217.71 /
      DATA ABE / 7.4206E-8, -1.5520E-8, -3.9403E-8,
     1           7.2155E-9, 1.1098E-8, -2.5501E-9,
     2           1.7136E-8, -1.3997E-8, -2.9656E-11,
     2           3.0777E-9, -1.1676E-10, -5.1930E-10,
     3           1.6039E-9, -6.4336E-10, -7.7804E-10,
     3           3.3527E-10, 2.1889E-10, -1.0600E-10,
     4           4.9587E-10, -3.6870E-10, -3.9284E-11,
     4           8.8928E-11, 8.0007E-12, -1.7114E-11 /
      DATA AABE / 2.1732E-7, 2.7494E-8, 2.7873E-9, 8.3163E-10 /
      DATA BTBE / -2.1649E-7, 2.8120E-7, 5.3031E-7,
     2            -1.5900E-8, 2.5911E-8, 2.4061E-8,
     3            -2.4333E-10, -9.7973E-9, 5.2094E-8,
     4            -4.5966E-10, 3.7317E-10, 2.5564E-9 /
C
C     BORON
C
      DATA EB / 8.3, 25.15, 37.93, 259.37, 340.22 /
      DATA AB / 5.8365E-8, 1.0047E-8, -3.6230E-8,
     1         -7.3448E-9, 1.0220E-8, 1.6951E-9,
     2          2.0590E-8, -9.8899E-9, -6.0949E-9,
     2          2.7762E-9, 1.6499E-9, -6.7692E-10,
     3          5.7023E-9, -4.6578E-9, -9.8688E-12,
     3          1.0242E-9, -3.8855E-11, -1.7281E-10,
     4          7.3539E-10, -2.9498E-10, -3.5672E-10,
     4          1.5372E-10, 1.0036E-10, -4.8602E-11,
     5          2.5458E-10, -1.8930E-10, -2.0169E-11,
     5          4.5657E-11, 4.1076E-12, -8.7867E-12 /
      DATA AAB / 3.0952E-7, 4.7980E-8, 9.1493E-9, 1.2780E-9,
     &             4.2697E-10 /
      DATA BTB / -4.9240E-7, 1.3750E-6, -2.5387E-6,
     2           -4.1361E-8, 5.5259E-8, 9.9841E-8,
     3           -5.2910E-9, 8.6226E-9, 8.0067E-9,
     4           -1.1156E-10, -4.4920E-9, 2.3885E-8,
     5           -2.3600E-10, 1.9159E-10, 1.3125E-9 /
C
C     CARBON
C
      DATA EC / 11.26, 24.38, 47.89, 64.49, 392.08, 489.98 /
      DATA AC / 5.9848E-8, 1.1903E-8, -3.0140E-8,
     1          -1.3693E-8, 8.3748E-9, 4.0150E-9,
     2          2.8395E-8, -1.6698E-8, -2.3557E-9,
     2          3.2161E-10, 9.6016E-10, 5.2713E-10,
     3          9.0555E-9, -6.3206E-9, -1.3256E-9,
     3          1.7441E-9, 3.2680E-10, -3.8303E-10,
     4          2.7464E-9, -2.0070E-9, -2.3595E-11,
     4          4.2011E-10, -8.1600E-11, -3.9729E-11,
     5          3.9495E-10, -1.5842E-10, -1.9158E-10,
     5          8.2555E-11, 5.3899E-11, -2.6102E-11,
     6          1.4715E-10, -1.0941E-10, -1.1657E-11,
     6          2.6389E-11, 2.3742E-12, -5.0786E-12 /
      DATA AAC / 3.7442E-7, 6.0150E-8, 1.4501E-8,
     &           6.0330E-9, 6.8634E-10, 2.4679E-10 /
      DATA BTC / -6.5826E-7, 2.0521E-6, -4.4694E-6,
     2           -4.0217E-8, -2.7908E-8, 5.5499E-7,
     3           -5.3596E-9, -1.0473E-8, 1.0629E-7,
     4           -5.7710E-9, 9.9302E-9, 7.4462E-9,
     5           -5.9916E-11, -2.4124E-9, 1.2828E-8,
     6           -1.3640E-10, 1.1074E-10, 7.5862E-10 /
C
C     NITROGEN
C
      DATA EN / 14.53, 29.60, 47.45, 77.47, 97.89, 552.06, 667.03 /
      DATA AN / 4.6209E-8, 9.2264E-9, -1.2092E-8,
     1          -2.4852E-8, 5.1361E-9, 8.3068E-9,
     2          2.4369E-8, -2.2155E-9, -1.4805E-8,
     2          -4.4218E-10, 4.5211E-9, 1.7874E-10,
     3          1.2964E-8, -8.3408E-9, -2.3684E-9,
     3          2.2485E-9, 2.6234E-10, -2.6333E-10,
     4          4.6322E-9, -3.4645E-9, -3.1014E-10,
     4          8.0576E-10, 5.7791E-11, -1.4907E-10,
     5          1.5862E-9, -9.8633E-10, 7.5130E-11,
     5          3.1005E-11, -8.1970E-12, 2.6759E-11,
     6          2.3635E-10, -9.4805E-11, -1.1465E-10,
     6          4.9404E-11, 3.2255E-11, -1.5620E-11,
     7          9.2653E-11, -6.8892E-11, -7.3402E-12,
     7          1.6616E-11, 1.4949E-12, -3.1978E-12 /
      DATA AAN / 2.7367E-7, 4.4690E-8, 1.0243E-8,
     &           7.9743E-9, 5.7812E-9, 4.1073E-10,
     &           1.5539E-10 /
      DATA BTN / -4.2976E-7, 9.8352E-7, -9.5745E-7,
     2           3.0430E-8, -5.2696E-7, 2.4863E-6,
     3           3.4218E-8, -2.9155E-7, 1.2324E-6,
     4           -4.9184E-9, 6.3763E-9, 1.4917E-8,
     5           -8.5163E-9, 1.8527E-8, -7.8181E-9,
     6           -3.5856E-11, -1.4437E-9, 7.6765E-9,
     7           -8.5888E-11, 6.9728E-11, 4.7767E-10 /
C
C     OXYGEN
C
      DATA EO / 13.62, 35.12, 54.93, 77.41, 113.90, 138.12, 739.32,
     &          871.39 /
      DATA AO / 3.3559E-8, 1.3449E-8, -6.7112E-9,
     1          -1.9976E-8, 1.6214E-9, 6.5852E-9,
     2          2.4476E-8, -5.3141E-9, -7.3316E-9,
     2          -4.4515E-9, 2.4257E-9, 1.9791E-9,
     3          1.4741E-8, -8.7905E-9, -8.6099E-10,
     3          -2.4143E-10, 1.2598E-10, 6.4901E-10,
     4          6.2130E-9, -2.5047E-9, -3.0813E-9,
     4          1.3559E-9, 8.6816E-10, -4.3189E-10,
     5           2.6145E-9, -2.0276E-9, -1.6569E-10,
     5          5.0245E-10, 3.0067E-11, -9.8231E-11,
     6          1.0099E-9, -6.5165E-10, 2.8863E-12,
     6          3.0336E-11, -1.4065E-11, 4.5106E-12,
     7          1.5258E-10, -6.1203E-11, -7.4014E-11,
     7          3.1894E-11, 2.0823E-11, -1.0084E-11,
     8          6.2090E-11, -4.6167E-11, -4.9189E-12,
     8          1.1135E-11, 1.0018E-12, -2.1430E-12 /
      DATA AAO / 3.2721E-7, 4.9003E-8, 1.7523E-8,
     &           1.0270E-8, 4.0023E-9, 2.6228E-9,
     &           2.6516E-10, 1.0413E-10 /
      DATA BTO / -6.7484E-7, 2.3938E-6, -6.04388E-6,
     2           2.1747E-8, -5.4612E-7, 2.7213E-6,
     3           2.8129E-8, -3.1692E-7, 1.4381E-6,
     4           4.8134E-10, -4.3936E-8, 2.1924E-7,
     5           -1.5957E-9, -7.2485E-10, 1.9724E-8,
     6           -2.6662E-9, 2.8968E-9, 1.7246E-8,
     7           -2.3148E-11, -9.3201E-10, 4.9557E-9,
     8           -5.7557E-11, 4.6727E-11, 3.2010E-10 /
C
C
      INSPC = KNUCZ + 1
C
      DO 3 JSPC = 1, INSPC
      PBFI(JSPC) = 0.0
      PBFIL(JSPC) = 0.0
    3 CONTINUE
C
      KFLAG = 0
C
C     ERROR EXIT IF ELEMENT NOT AVAILABLE
C
      IF(KNUCZ .GE. 1 .AND. KNUCZ .LE. 8) GO TO 5
      KFLAG = 1
      RETURN
C
C     BRANCH ON ELEMENT TO CALL FIT EVALUATION ROUTINE WITH PROPER
C     CONSTANTS
C
    5 GO TO (10, 20, 30, 40, 50, 60, 70, 80) KNUCZ
C
   10 CALL ADCBFIT(KNUCZ, PTE, PNE, EH, AH, AAH, BTH, PBFI, PBFIL)
      RETURN
   20 CALL ADCBFIT(KNUCZ, PTE, PNE, EHE, AHE, AAHE, BTHE, PBFI,PBFIL)
      RETURN
   30 CALL ADCBFIT(KNUCZ, PTE, PNE, ELI, ALI, AALI, BTLI, PBFI,PBFIL)
      RETURN
   40 CALL ADCBFIT(KNUCZ, PTE, PNE, EBE, ABE, AABE, BTBE, PBFI,PBFIL)
      RETURN
   50 CALL ADCBFIT(KNUCZ, PTE, PNE, EB, AB, AAB, BTB, PBFI, PBFIL)
      RETURN
   60 CALL ADCBFIT(KNUCZ, PTE, PNE, EC, AC, AAC, BTC, PBFI, PBFIL)
      RETURN
   70 CALL ADCBFIT(KNUCZ, PTE, PNE, EN, AN, AAN, BTN, PBFI, PBFIL)
      RETURN
   80 CALL ADCBFIT(KNUCZ, PTE, PNE, EO, AO, AAO, BTO, PBFI, PBFIL)
      RETURN
      END




      SUBROUTINE ADCBFIT(KNUCZ, PTE, PNE, PE, PA, PAA, 
     1                   PBT, PBFI, PBFIL)
C
C*******************************************************************
C     EVALUATE THE BELFAST GROUP RATE COEFFICIENT FITTING FORMULAS.
C     NOTE THAT 0 IS RETURNED FOR TE < (IONIZATION POTENTIAL / 10).
C*******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      DIMENSION PE(KNUCZ), PA(6,KNUCZ), PAA(KNUCZ), PBT(3,KNUCZ)
      DIMENSION PBFI(1), PBFIL(1)
C
      ZJLEV = 1.6021E-19
C
C     CONVERT TE TO EV AND ZERO RETURN ARRAYS
C
      ZTE = 1000.0 * PTE
C
      INSPC = KNUCZ + 1
      DO 10 JSPC = 1, INSPC
      PBFI(JSPC) = 0.0
      PBFIL(JSPC) = 0.0
   10 CONTINUE
C
C     LOOP OVER ALL SPECIES EXCEPT FULLY STRIPPED
C
      DO 100 JSPC = 1, KNUCZ
C
C     RATE IS ZERO FOR KT < (IP/10)
C
      IF(ZTE .LT. PE(JSPC)/10.0) GO TO 100
C
C     BRANCH IF HIGH TEMPERATURE FIT IS REQUIRED
C
      IF(ZTE .GT. 10.0 * PE(JSPC)) GO TO 80
C
C     (IP/10) < TE < (10*IP) FITTING FORM
C
      ZX = PE(JSPC) / ZTE
      ZLOG = DLOG10(1.0 / ZX)
C
      ZSUM = PA(1,JSPC)
      DO 20 J = 2, 6
      ZSUM = ZSUM + PA(J,JSPC) * (ZLOG**(J-1))
   20 CONTINUE
C
      ZIRC = DEXP(-ZX) * DSQRT(1.0 / ZX) * ZSUM
      PBFI(JSPC) = PNE * ZIRC
      PBFIL(JSPC) = ZJLEV * PE(JSPC) * PBFI(JSPC)
      GO TO 100
C
C     HIGH TEMPERATURE FIT FOR TE > 10*IP
C
   80 ZX = PE(JSPC) / ZTE
C
      ZIRC = DSQRT(ZX) * (PAA(JSPC) * DLOG(1.0/ZX)
     &     + PBT(1,JSPC) + PBT(2,JSPC) * ZX + PBT(3,JSPC) * ZX * ZX)
C
      PBFI(JSPC) = PNE * ZIRC
      PBFIL(JSPC) = ZJLEV * PE(JSPC) * PBFI(JSPC)
C
  100 CONTINUE
C
      RETURN
      END



      SUBROUTINE ADCYGR(KNUCZ, PTE, PNE, PYGRI, PYGRL, KFLAG)
C
C*****************************************************************
C     COMPUTE DIRECT ELECTRON IONIZATION RATES FOR HYDROGEN-LIKE
C     THROUGH NEON-LIKE CHARGE STATES OF CERTAIN HEAVY ELEMENTS
C     (PRESENTLY SCANDIUM AND IRON).  USES FORMULAE DUE TO
C     S. YOUNGER, NBS.  (JQSRT 29, 61 (1983) AND PRIVATE
C     COMMUNICATIONS).
C
C     KNUCZ= ELEMENT (Z=21 OR 26 FOR NOW)
C     PTE = ELECTRON TEMPERATURE(KEV)
C     PNE = ELECTRON DENSITY (CM-3)
C
C     PYGRI = ARRAY OF IONIZATION RATES (SEC-1)
C             (ONLY THE 10 STATES FROM H-LIKE TO NE-LIKE WILL BE
C              OVERWRITTEN BY THE SUBROUTINE)
C     PYGRL = ARRAY OF IONIZATION ENERGY LOSS RATES (WATTS/ION)
C     KFLAG = 0 FOR O.K.;
C           = 1 FOR INVALID ELEMENT (NO VALUES ARE RETURNED)
C******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
      COMMON / COMYGR / EION(10,3), EIONFE(10,3), EIONSC(10,3),
     &          NSS(10,3), RATE(10,3),
     &          CA(3), CB(3), CC(3), CD(3)
C
      DIMENSION PYGRI(1), PYGRL(1)
C
      EXTERNAL ADCFCHI
C
      DATA ((NSS(JBE,JSS),JBE=1,10),JSS=1,3)/ 1,9*2,0,0,1,7*2,
     &      4*0,1,2,3,4,5,6/
C
      DATA (CA(J),J=1,3)/7.5E-14, 5.9E-14, 8.4E-14/
      DATA (CB(J), J=1,3)/ -2.605E-14, -1.635E-14, -2.68E-14/
      DATA (CC(J), J=1,3)/ 1.195E-14, 0.82E-14, 0.672E-14/
      DATA (CD(J), J=1,3)/ -6.15E-14, -3.79E-14, -5.08E-14/
C
      DATA ((EIONSC(JBE,JSS),JBE=1,10),JSS=1,3)/ 6034., 5675.,
     &     5566., 5491., 5318., 5204., 5095., 4990., 4890., 4794.,
     &     0., 0., 1288., 1213., 1135., 1068., 993.3,
     &     914.0, 846.3, 769.3, 4*0., 1090., 1004., 922.4,
     &     838.6, 761.8, 687.4 /
C
      DATA ((EIONFE(JBE,JSS),JBE=1,10),JSS=1,3)/ 9277., 8829.,
     &     8694., 8521., 8400., 8263., 8118., 7949., 7818., 7824.,
     &     0., 0., 2046., 1950., 1852., 1755., 1679., 1566., 1482.,
     &     1394., 4*0., 1789., 1678., 1583., 1463., 1366., 1266./
C
C
C
      TEV = 1000. * PTE
C
      EION(1,1) = 0.0
      DO 20 JSS = 1, 3
      DO 20 JBE = 1, 10
      IF(KNUCZ .EQ. 21) EION(JBE,JSS) = EIONSC(JBE,JSS)
      IF(KNUCZ .EQ. 26) EION(JBE,JSS) = EIONFE(JBE,JSS)
   20 CONTINUE
C
      KFLAG = 0
      IF(EION(1,1) .NE. 0.0)GO TO 50
      KFLAG = 1
      RETURN
C
   50 DO 100 JBE = 1, 10
      ICS = KNUCZ - JBE +1
      PYGRI(ICS) = 0.0
      PYGRL(ICS) = 0.0
C
      DO 90 JSS = 1, 3
      RATE(JBE,JSS) = 0.0
      IF(NSS(JBE,JSS) .EQ. 0) GO TO 90
C
      CHI = TEV / EION(JBE,JSS)
      A = CA(JSS) * DBLE(NSS(JBE,JSS))
      B = CB(JSS) * DBLE(NSS(JBE,JSS))
      C = CC(JSS) * DBLE(NSS(JBE,JSS))
      D = CD(JSS) * DBLE(NSS(JBE,JSS))
      RATE(JBE,JSS) = ADCFCHI(CHI,A,B,C,D) * (EION(JBE,JSS)**(-1.5)) *
     &                2.2E-06 * DSQRT(CHI) * DEXP(-1.0/CHI)
      PYGRI(ICS) = PYGRI(ICS) + PNE * RATE(JBE,JSS)
      PYGRL(ICS) = PYGRL(ICS) + 1.6021E-19 * PNE * EION(JBE,JSS)
     &              * RATE(JBE,JSS)
   90 CONTINUE
  100 CONTINUE
C
      RETURN
      END

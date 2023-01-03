SUBROUTINE PRODUC(NLAYER,DENSTY,SRR,PH,DAYLEN,LATHIK,TEMP,PTORSP,HR)
#	***********************************************************
# Subroutine PRODUC Ver. 2.14
#		Calculates productivity of a seagrass canopy
#		using both spectral and scalar equations.  
#		Usage: CALL PRODUC(NLAYER,DENSTY,PH,DAYLEN,LATHIK)
#
# Definitions of passed variables:
# ----------------------------------------------------
# 		NLAYER - Number of layers in the canopy (default = 100)
#  		DENSTY - Shoot density (shoots/m2)
# 		SRR - Shoot:Root ratio
#   	PH - Seawater pH
# 		DAYLEN - Length of the daylight period (hours)
# 		LATHICK - Width (thickness) of each canopy layer (mm)
#     TEMP - Water temperature, ∞C
#
#	COMMOM variables are array variables shared with other subroutines and the main program
#
# Definition of common variables:
# ----------------------------------------------------
# 		BIOMAS(I) - Relative biomass present in canopy layer (I).  All
#         				biomass elemenct should sum to 1.  Calculated in
#           			subroutine BIODIS
# 	  LAI(I)    - Leaf area index (m2 leaf/m2 area) in each canopy layer (I).
#         				Determined in subroutine BIODIS
# 		LAP(I)    - Leaf area projected toward the zenith.  Calculated as a
#         				function of bending angle in subroutine BIODIS 
# 		WAVLED(J) - Wavelength of spectral element j, determined in 
#       					subroutine EDCALC from irradiance inupt file
# 		PIZ(I)    - Spectral production coefficient Pi(z) calculated in
#           			subroutine EDCALC
# 		EO(I)     - Scalar irradiance (Ed + Eu) in canopy layer I calculated in
#         				subroutine EDCALC
#
# Files required by PROD_210
# ----------------------------------------------------
#     None
#
# Files created by PROD_201
# ----------------------------------------------------
#     None
#
# Files modified by PROD_201
# ----------------------------------------------------
#     GLRESULT.DAT (or user-specified name of optupt file)
#
# CHECK PAST FILE FOR UPDATE NOTES
# ********************************************************************************

	IMPLICIT REAL*8 (A-Z)
	INTEGER I,J,K,L,M,N,NLAYER,NUMINT,NWAVED,NWAVRB,NWAVLF,
     1	NCOL,MROW,ITER,DONE,FINL
	CHARACTER ANS,RUN,HDR*64
      CHARACTER*64 COPYRITE,VER
	COMMON /MASS/ASYMP,INFLEC,SHAPE,BETA0,BIOMAS(100)
     1	/AREA/LAI(100)
     2    /PROJEC/LAP(100)
     3	/WAVE/WAVLED(100)
     4    /PROD/PIZ(100),ALPHA(100)
     5    /IRR/EO(100)
     6    /CONVRG/ITER,DONE,FINL,RUN
	DIMENSION SCALRP(100),SPECP(100)
	DATA COPYRITE/'Copyright (c) 2016 by Richard C. Zimmerman'/
	DATA VER/'2.14'/
      PI=355./113.    !Define the constant pi
C
C				*** INITIAL PARAMETERIZATION ***
C			*** NOTE: Pm AND R ARE RELATIVE (UNITLESS) ***
C
c	WRITE(*,2)VER
c2	FORMAT(' Subroutine PRODUC Ver. ',A64)
	SUMBIO=0.
	SUMLAI=0.
	DO 5 I=1,NLAYER
		SUMBIO=SUMBIO+BIOMAS(I)
		SUMLAI=SUMLAI+LAI(I)
5	CONTINUE
	PMAX=81.988*EXP(-0.5316*PH) !Pmax is unitless (per hour) for easy scaling (= 1 at pH 8.2)
      T_INT = 727.76*EXP(-0.5504*PH)  !pH control of intercept for leaf P:R vs Temperature
      T_EXP_SLOPE = 0.0491-0.0093*PH  !pH control of exponential slope for leaf P:R vs Temp
      PTOR_LEAF = T_INT*EXP(T_EXP_SLOPE*TEMP) !Leaf P:R now depends on Temp and pH
      LEAFR = PMAX/PTOR_LEAF  !%Resulting temperature dependancy of leaf R
c      write(*,100)ph,pmax,temp,ptor_leaf,leafr
c100   format('pH =',f10.3,' pmax=',f10.3,'temp= ',f10.3,'p:R=',f10.3,
c     1    'Leaf R=',f10.3)           
	ROOTR=LEAFR/2.      !Root respiration, half of leaf R
	RHIZR=ROOTR/2.      !Rhizome respiration, half of root R
	BIOTOT=SRR+1
	LEAFBIO=SRR/BIOTOT
	ROOTBIO=0.5*(1.-LEAFBIO)
	RHIZBIO=ROOTBIO
C
C		*** CALCULATE RELATIVE RESPIRATION OF LEAVES ***
C					ROOTS AND RHIZOME
C			 note: RDAY also equals the daily
C				Hsat requirement of the shoots
C	
      RNOW=LEAFR+ROOTBIO*ROOTR+RHIZBIO*RHIZR  !Instantaneous respiration demand
	RDAY=(LEAFR*24.)+ROOTBIO*(DAYLEN*ROOTR+(24.-DAYLEN)*
     1	0.65*ROOTR)+RHIZBIO*(DAYLEN*RHIZR+(24.-DAYLEN)*0.65*RHIZR)

C
C			*** CALCULATE INSTANTANEOUS SPECTRAL PRODUCTION ***
C				FOR EACH LAYER OF THE CANOPY
C
	PNOW=0.
	DO 300 I=1,NLAYER
		SPECP(I)=BIOMAS(I)*PMAX*(1.-EXP(-PIZ(I)/PMAX))
		PNOW=PNOW+SPECP(I)
c	write(*,3010)i,piz(i),biomas(i),specp(i),pnoon
c3010	format(' Layer ',i3,'Pi(z)=',f7.3,' Biomass=',f7.3,' Spec P='
c     1	f7.3,' Integrated P=',f7.3)
c	if(mod(i,10).eq.0)pause
300	CONTINUE
C
C			*** Old Way to Integrate DAILY SPECTRAL PRODUCTION ***
C                             before Ver 2.10
C	OLDPSP=0.
C	DO 320 I=1,NLAYER
C		OLDPSP=OLDPSP+(BIOMAS(I)*PMAX*(1.-EXP(-0.67*PIZ(I)
C    1		/PMAX)))*DAYLEN
C320   CONTINUE
C
C			*** New Way Integrate DAILY SPECTRAL PRODUCTION ***
C                             Ver 2.10 and on
	IF(HR.EQ.12)THEN        !Only if input irradiance is at noon
          DAYPSP=0.
          TODINC=DAYLEN*6.   !# 10 min increments per day for daily integration
          TIMINC=1./TODINC        !Time increment, fractional days
          IF(MOD(TODINC,DAYLEN).NE.0)TODINC=TODINC+1.    !add 1 to counter if remainder
	    DO 400 I=1,NLAYER
              TOD=0.
              DO 350 J=1,INT(TODINC)
                  TOD=TOD+TIMINC
		        DAYPSP=DAYPSP+(BIOMAS(I)*PMAX*
     1                (1.-EXP(-SIN(PI*TOD)*PIZ(I)/PMAX)))/6.  !divide by 6 to properly integrate 10 min segments
C      write(*,3010)J,K,TOD,DAYLEN,DAYPSP,OLDPSP
C3010  format(' J=',I5,' K=',I5,' TOD='F10.3,' DAYLEN=',F10.3,
C     1' DAYPSP=',F10.3,' OLDPSP=',F10.3)
350           CONTINUE
C          PAUSE
400       CONTINUE
      ENDIF

C
C			*** CALCULATE DAILY P:R ***
C
	IF(SUMBIO.GT.0.)THEN
C		PTORSC=DAYPSC/RDAY
		PTORSP=DAYPSP/RDAY
	ELSE
C		PTORSC=-99.
		PTORSP=-99.
      ENDIF
      PTORTST=ABS(PTORSP-1.0)
C
C		*** WRITE OUTPUT DATA FILE ***
C
c	write(*,401)run,iter
c401   format('RUN=',A1,' ITER=',I3)
      IF(RUN.EQ.'o'   .OR.    RUN.EQ.'O')THEN
          WRITE(*,402)ITER,PTORTST,DENSTY
402       FORMAT('Iter ',I4,5x,' P:R test =',F8.5,5x,' Density =',F6.0)
          IF(PTORTST.GT.0.0001    .AND.   DENSTY.GT.1.0)GO TO 550
          DONE=1
          IF(PTORTST.LT.0.0001)THEN
              FINL=1
              WRITE(*,403)
              WRITE(3,403)
403           FORMAT(//'Daily P:R = 1.0, condition satisfied')
              GO TO 410
          ELSE IF(DENSTY.LT.1.)THEN
              FINL=1
              WRITE(*,404)
              WRITE(3,404)
404           FORMAT(//'Light environment cannot sustain',
     1            ' positive P:R')
          END IF
      END IF
410   WRITE(*,500)DENSTY,SUMLAI,SUMBIO
	WRITE(*,501)PH,DAYLEN,SRR,RDAY,PMAX,PNOW,RNOW
	WRITE(*,502)DAYPSP,PTORSP
	WRITE(3,500)DENSTY,SUMLAI,SUMBIO
	WRITE(3,501)PH,DAYLEN,SRR,RDAY,PMAX,PNOW,RNOW
	WRITE(3,502)DAYPSP,PTORSP
500	FORMAT(///' "DENSITY (shoots/m2)"',2X,F10.0/
     1    ' "LAI (m2/m2)"',13X,F10.3/
     2    ' "BIOMASS (rel)"',12X,F10.4)
501   FORMAT(' "pH"',21X,F10.2/' "Photoperiod (h)"',8X,F10.2
     1    /' "Shoot:Root"',14X,F10.3
     2    /' "Hsat requirement (h)"',4X,F10.3
     3    /' "Pmax (rel)"',14X,F10.3
     4    /' "Canopy Integrated P"',3x,F10.3
     5    /' "Canopy Integrated R"',3x,F10.3)
502	FORMAT(' "Daily P (Hsat equiv)"',4X,F10.3/
     1    ' "Whole plant P:R (per d)"',F10.2/)
      WRITE(3,510)
510	FORMAT(///' "Parameter values below are for',
     1    ' specified height layers; they are NOT',
     2    ' canopy integrals."'/' Height_(m)',9X,
     3    'Biomass_rel)',7X,'LAI_(m2/m2)',9X,
     4	'LAP_(m2/m2)',9x,'Noon_P(rel)')
      DO 540 I=NLAYER,1,-10
	    HT=(LATHIK*FLOAT(I))/1000.
	    WRITE(3,530)HT,BIOMAS(I),LAI(I),LAP(I),SPECP(I)
530	    FORMAT(1X,F7.2,6(5X,G15.3))
540   CONTINUE
      WRITE(3,545)
545   FORMAT(1X,'------------------------------------'///)
550	RETURN
	END
